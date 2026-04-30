# test/util/thermodynamic_identities.jl — model-agnostic self-validation harness.
#
# Issue #117: every QAtlas implementation that exposes (Energy, FreeEnergy,
# ThermalEntropy, SpecificHeat, ...) at a given BC should obey universal
# thermodynamic identities like
#
#   ε(β) = f(β) + T·s(β)               (Gibbs)
#   c_v(β) = -β² ∂ε/∂β                 (specific heat from energy)
#   m_α(h) = -∂f/∂h_α                  (magnetisation from free energy)
#   χ_αα   = ∂m_α/∂h_α = β·Var(M_α)/N  (linear response from variance)
#
# These hold *between the model's own outputs* — we are not comparing to
# literature here; we are checking internal consistency of the dispatch
# layer.  Self-validation catches per-site/total drift, sign errors,
# missing field-dependence, and ForwardDiff-incompatible kwargs.
#
# This file deliberately lives in `test/util/`: the harness is a
# debugging/verification tool, not part of QAtlas's public surface.  If
# downstream packages need it later, lift it to `src/verification/` as a
# weakdep extension on ForwardDiff.

using ForwardDiff
using QAtlas:
    fetch,
    AbstractQAtlasModel,
    AbstractQuantity,
    BoundaryCondition,
    Energy,
    FreeEnergy,
    ThermalEntropy,
    SpecificHeat,
    OBC,
    PBC,
    Infinite

"""
    ThermoIdentity(name, requires, check; model_filter = _ -> true)

A single self-validation rule.

- `name::String`: human-readable label, surfaced in `IdentityCheckResult`.
- `requires::Vector{Type}`: quantity types the identity needs to be
  evaluable on the target `(model, bc)`.  The harness checks dispatch
  existence (via `which(fetch, ...)`) and skips the identity if any
  required type lacks a non-catch-all method.
- `check::Function`: `(model, bc, params::NamedTuple) -> (lhs, rhs)`.
  The two sides should be equal up to the harness's `(rtol, atol)`.
- `model_filter::Function`: predicate `model -> Bool` for restricting
  the identity to specific model classes (e.g. only models that carry an
  `h` field for field-perturbation identities).  Default `_ -> true`
  applies the identity to every model that satisfies `requires`.
"""
struct ThermoIdentity
    name::String
    requires::Vector{Type}
    check::Function
    model_filter::Function
end

ThermoIdentity(name, requires, check; model_filter=_ -> true) =
    ThermoIdentity(name, requires, check, model_filter)

"""
    IdentityCheckResult

One row per `(identity, params)` evaluated by
[`verify_thermodynamic_identities`](@ref).

`status` is `:pass`, `:fail`, or `:skipped` (the latter when one of
`identity.requires` is not dispatchable for `(model, bc)`).  For
`:skipped` rows, the numerical fields are `NaN`.
"""
struct IdentityCheckResult
    model::Any
    bc::Any
    identity::String
    params::NamedTuple
    lhs::Float64
    rhs::Float64
    abs_err::Float64
    rel_err::Float64
    status::Symbol
end

# ──────────────────────────────────────────────────────────────────────
# Default identity set
# ──────────────────────────────────────────────────────────────────────

"""
    GIBBS_RELATION

Self-validation of `ε = f + T·s` (no AutoDiff — three independent
fetches whose values must reconcile).  Catches per-site/total mix-ups,
sign errors in entropy, and missing temperature factors.
"""
const GIBBS_RELATION = ThermoIdentity(
    "Gibbs ε = f + T·s",
    Type[Energy{:per_site}, FreeEnergy, ThermalEntropy],
    function (model, bc, params)
        β = params.β
        T = 1 / β
        ε = fetch(model, Energy(:per_site), bc; beta=β)
        f = fetch(model, FreeEnergy(), bc; beta=β)
        s = fetch(model, ThermalEntropy(), bc; beta=β)
        return Float64(ε), Float64(f + T * s)
    end,
)

"""
    SPECIFIC_HEAT_FROM_ENERGY

Self-validation of `c_v = -β² ∂ε/∂β` via `ForwardDiff.derivative` on
`Energy(:per_site)`.  Requires the model's `fetch` methods to accept a
`Real` `beta` kwarg (TFIM was relaxed in PR #115; new models inherit
this requirement).
"""
const SPECIFIC_HEAT_FROM_ENERGY = ThermoIdentity(
    "c_v = -β² ∂ε/∂β  (ForwardDiff)",
    Type[Energy{:per_site}, SpecificHeat],
    function (model, bc, params)
        β = params.β
        dε_dβ = ForwardDiff.derivative(b -> fetch(model, Energy(:per_site), bc; beta=b), β)
        c_v = fetch(model, SpecificHeat(), bc; beta=β)
        return -β^2 * Float64(dε_dβ), Float64(c_v)
    end,
)

"""
    SPECIFIC_HEAT_FROM_ENTROPY

Cross-method check `c_v = T · ∂s/∂T = -β · ∂s/∂β` via `ForwardDiff` on
`ThermalEntropy`.  Equivalent to `SPECIFIC_HEAT_FROM_ENERGY` only modulo
the Gibbs relation, so a discrepancy here pinpoints which of (s, ε)
disagrees with the c_v implementation.
"""
const SPECIFIC_HEAT_FROM_ENTROPY = ThermoIdentity(
    "c_v = -β · ∂s/∂β  (ForwardDiff on ThermalEntropy)",
    Type[ThermalEntropy, SpecificHeat],
    function (model, bc, params)
        β = params.β
        ds_dβ = ForwardDiff.derivative(b -> fetch(model, ThermalEntropy(), bc; beta=b), β)
        c_v = fetch(model, SpecificHeat(), bc; beta=β)
        return -β * Float64(ds_dβ), Float64(c_v)
    end,
)

# ──────────────────────────────────────────────────────────────────────
# Field-perturbation identities (require the model to carry an `h` field
# that can be reconstructed by positional argument)
# ──────────────────────────────────────────────────────────────────────

"""
    _perturb_field(model, field::Symbol, val) -> Union{typeof(model),Nothing}

Reconstruct `model` with the named `field` replaced by `val`, using the
positional constructor `typeof(model)(getfield.(model, propertynames)...)`.
Returns `nothing` if the field is absent or the constructor signature
does not match — the harness treats either as a skip signal.

ForwardDiff-friendly: the returned model carries `Dual` numbers in the
perturbed field if `val` is `Dual`.
"""
function _perturb_field(model, field::Symbol, val)
    field in propertynames(model) || return nothing
    args = map(propertynames(model)) do f
        f == field ? val : getfield(model, f)
    end
    try
        return typeof(model).name.wrapper(args...)
    catch
        return nothing
    end
end

_has_h_field(model) = :h in propertynames(model)

"""
    _central_diff(f, x; δ) -> Float64

Symmetric central finite difference `(f(x + δ) − f(x − δ)) / (2δ)`.
Used in place of `ForwardDiff` for identities that perturb a model's
*physical* field (e.g. `h`), since concrete-typed model structs
(`TFIM`, …) cannot store the `Dual` numbers that ForwardDiff would push.
The `O(δ²)` truncation error sets the achievable atol; default `δ`
balances against `O(eps/δ)` round-off (`eps^{1/3} ≈ 6e-6` is optimal).
"""
function _central_diff(f, x; δ::Real=1e-5)
    return (f(x + δ) - f(x - δ)) / (2δ)
end

"""
    MAGNETIZATION_X_FROM_FREE_ENERGY

Cross-method check `m_x = -∂f/∂h` (Helmholtz identity) via central
finite difference on `FreeEnergy` w.r.t. the model's transverse field.
Skipped on models without an `h` field (XXZ1D, S1Heisenberg1D, …).
For TFIM specifically this cross-validates the Pfeuty closed-form `m_x`
against a derivative of the free-fermion free energy.

Tolerance: limited by the `O(δ²) ~ 1e-10` central-difference truncation
error, so set the harness atol no tighter than `1e-7`.
"""
const MAGNETIZATION_X_FROM_FREE_ENERGY = ThermoIdentity(
    "m_x = -∂f/∂h  (central diff on FreeEnergy w.r.t. h)",
    Type[FreeEnergy, MagnetizationX],
    function (model, bc, params)
        β = params.β
        h0 = getfield(model, :h)
        df_dh = _central_diff(h0) do h
            perturbed = _perturb_field(model, :h, h)
            return fetch(perturbed, FreeEnergy(), bc; beta=β)
        end
        m_x = fetch(model, MagnetizationX(), bc; beta=β)
        return -Float64(df_dh), Float64(m_x)
    end;
    model_filter=_has_h_field,
)

"""
    SUSCEPTIBILITY_XX_KUBO_FROM_MAGNETIZATION

Cross-method check `χ_xx = ∂m_x/∂h` (Kubo / static linear response)
via central finite difference on `MagnetizationX` w.r.t. the model's
transverse field.

**Convention warning.**  This identity tests the *Kubo* static
susceptibility, which equals `(β/N) ∫₀^β dτ ⟨M_x(τ)M_x(0)⟩_c` — the
imaginary-time-integrated correlator at ω=0.  For a quantum
Hamiltonian where `[M_x, H] ≠ 0` the Kubo χ differs from the
*equal-time* variance `β·Var(M_x)/N` by operator-ordering corrections.

QAtlas's OBC dense-ED implementations of `SusceptibilityXX` (TFIM,
XXZ1D, S1Heisenberg1D) return the **variance form**, so this identity
is **not** in `DEFAULT_IDENTITIES` to avoid spurious failures on those
backends.  Only the closed-form `Infinite()` Kubo paths (e.g. TFIM
Infinite via the Calabrese-Mussardo integral) pass.  Use this identity
explicitly via the `identities=[…]` kwarg of
`verify_thermodynamic_identities` when validating a Kubo-convention
implementation.
"""
const SUSCEPTIBILITY_XX_KUBO_FROM_MAGNETIZATION = ThermoIdentity(
    "χ_xx = ∂m_x/∂h  (central diff, Kubo convention)",
    Type[MagnetizationX, SusceptibilityXX],
    function (model, bc, params)
        β = params.β
        h0 = getfield(model, :h)
        dmx_dh = _central_diff(h0) do h
            perturbed = _perturb_field(model, :h, h)
            return fetch(perturbed, MagnetizationX(), bc; beta=β)
        end
        χ_xx = fetch(model, SusceptibilityXX(), bc; beta=β)
        return Float64(dmx_dh), Float64(χ_xx)
    end;
    model_filter=_has_h_field,
)

"""
    DEFAULT_IDENTITIES

The identity set evaluated by `verify_thermodynamic_identities` when
the caller does not pass an explicit `identities` kwarg.  Includes the
universal Gibbs / specific-heat-from-energy / specific-heat-from-entropy
checks plus the Helmholtz `m_x = -∂f/∂h` check (skipped on models
without an `h` field).

The Kubo-vs-variance susceptibility identity
[`SUSCEPTIBILITY_XX_KUBO_FROM_MAGNETIZATION`](@ref) is *not* in this
set because the QAtlas OBC dense-ED backends use the equal-time
variance convention; see that identity's docstring for details.
"""
const DEFAULT_IDENTITIES = ThermoIdentity[
    GIBBS_RELATION,
    SPECIFIC_HEAT_FROM_ENERGY,
    SPECIFIC_HEAT_FROM_ENTROPY,
    MAGNETIZATION_X_FROM_FREE_ENERGY,
]

# ──────────────────────────────────────────────────────────────────────
# Dispatch-existence helper
# ──────────────────────────────────────────────────────────────────────

# Capture the catch-all `fetch(::AbstractQAtlasModel, ::AbstractQuantity,
# ::BoundaryCondition; ...)` once at file load; any (model, quantity, bc)
# triple whose `which(fetch, ...)` returns this same Method object is
# *not* implemented (the catch-all just throws an informative error).
const _CATCH_ALL_FETCH_METHOD = which(
    fetch, Tuple{AbstractQAtlasModel,AbstractQuantity,BoundaryCondition}
)

function _has_dispatch(model, ::Type{Q}, bc) where {Q}
    return which(fetch, Tuple{typeof(model),Q,typeof(bc)}) !== _CATCH_ALL_FETCH_METHOD
end

function _can_run(identity::ThermoIdentity, model, bc)
    identity.model_filter(model) || return false
    return all(Q -> _has_dispatch(model, Q, bc), identity.requires)
end

# ──────────────────────────────────────────────────────────────────────
# Harness
# ──────────────────────────────────────────────────────────────────────

"""
    verify_thermodynamic_identities(model, bc;
                                     βs,
                                     identities=DEFAULT_IDENTITIES,
                                     rtol=1e-8, atol=1e-10)
        -> Vector{IdentityCheckResult}

For every `(identity, β)` pair, evaluate `identity.check(model, bc, (;β))`
and record whether `lhs ≈ rhs` within `(rtol, atol)`.  Identities that
require a quantity not dispatchable on `(model, bc)` are recorded as
`:skipped` (numeric fields `NaN`).

The skip semantics let the same call be made on any model — including
ones missing some of the registry — without spurious test failures.
This is what makes the harness genuinely model-agnostic.

Use from `@testset`s with `@test all(r.status === :pass for r in results)`
when every required quantity is implemented (TFIM today), or
`@test all(r.status !== :fail for r in results)` for partially-populated
models.
"""
function verify_thermodynamic_identities(
    model::AbstractQAtlasModel,
    bc::BoundaryCondition;
    βs::AbstractVector{<:Real},
    identities::AbstractVector{ThermoIdentity}=DEFAULT_IDENTITIES,
    rtol::Real=1e-8,
    atol::Real=1e-10,
)
    results = IdentityCheckResult[]
    for identity in identities
        runnable = _can_run(identity, model, bc)
        for β in βs
            params = (; β=β)
            if !runnable
                push!(
                    results,
                    IdentityCheckResult(
                        model, bc, identity.name, params, NaN, NaN, NaN, NaN, :skipped
                    ),
                )
                continue
            end
            lhs, rhs = identity.check(model, bc, params)
            abs_err = abs(lhs - rhs)
            rel_err = abs_err / max(abs(lhs), abs(rhs), eps())
            status = (abs_err ≤ atol || rel_err ≤ rtol) ? :pass : :fail
            push!(
                results,
                IdentityCheckResult(
                    model, bc, identity.name, params, lhs, rhs, abs_err, rel_err, status
                ),
            )
        end
    end
    return results
end
