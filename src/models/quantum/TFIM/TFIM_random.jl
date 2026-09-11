# ─────────────────────────────────────────────────────────────────────────────
# Random transverse-field Ising chain: the Griffiths dynamical exponent.
#
# `TFIM` is the clean chain and carries two scalars, which is why the exact
# Griffiths result could not live there: it is a condition on the DISTRIBUTIONS
# of the couplings. That condition, [(J/h)^{1/z}]_av = 1, is itself
# distribution-free and model-independent in shape; the vocabulary it needs
# (`DisorderFamily`, `Disordered`) is in `core/disorder.jl`, and only its
# evaluation for THIS chain is here.
#
# Reference: Iglói-Monthus, [IgloiMonthus2005](@cite). Equation and section
# numbers are those of the arXiv version, `cond-mat/0502448`.
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Pauli σ, as in `TFIM.jl`. H = -Σ Jᵢ σᶻσᶻ - Σ hᵢ σˣ
#   Disorder:    Jᵢ = J·λᵢ and hᵢ = h·μᵢ, so the TFIM's J and h are the SCALES.

"""
    RandomTFIM

`Disordered{TFIM}`: the 1D random transverse-field Ising chain,

    H = -Σ_i J_i σᶻ_i σᶻ_{i+1} - Σ_i h_i σˣ_i

with `J_i = J·λ_i` and `h_i = h·μ_i`, the TFIM's own `J` and `h` being the
SCALES and `λ`, `μ` drawn from the families named for `:J` and `:h`.

Criticality is `[ln J]_av = [ln h]_av`, which for equal families is `J == h` and
in general is not.  Distance from it is [`rtfim_delta`](@ref).

| Quantity | BC | Coverage |
| --- | --- | --- |
| [`DynamicalExponent`](@ref) | `Infinite` | exact off criticality; throws at `δ = 0` |
| [`ActivatedExponent`](@ref) | `Infinite` | `1/2`, at criticality only |
| [`UniversalityClass`](@ref) | `Infinite` | `:IsingSDRG`, at criticality only |

Anything else is refused: a disordered chain does not inherit the clean one's
answers.
"""
const RandomTFIM = Disordered{TFIM}

"""
    RandomTFIM(; J = 1.0, h = 1.0, D = 1.0)

The symmetric [`PowerLawDisorder`](@ref) case, which is the default in the
strong-disorder RG literature.  Any other combination is built with
[`Disordered`](@ref) directly.
"""
function RandomTFIM(; J::Real=1.0, h::Real=1.0, D::Real=1.0)
    return Disordered(TFIM(; J=J, h=h); J=PowerLawDisorder(D), h=PowerLawDisorder(D))
end
export RandomTFIM

"""
    rtfim_delta(m::RandomTFIM) -> Float64

`δ = ([ln h]_av − [ln J]_av) / (var[ln h] + var[ln J])`.  Zero exactly at the
infinite-randomness critical point; positive in the disordered phase.
"""
function rtfim_delta(m::RandomTFIM)
    num =
        (log(clean_model(m).h) + mean_log(disorder(m, :h))) -
        (log(clean_model(m).J) + mean_log(disorder(m, :J)))
    return num / (var_log(disorder(m, :h)) + var_log(disorder(m, :J)))
end
export rtfim_delta

"""
    _rtfim_griffiths_residual(u, m)

`ln [(J/h)^{1/z}]_av` at `u = 1/z` ([IgloiMonthus2005](@cite) Eq. (4.15),
§4.1.3), which the condition sets to zero.

Solved in `u`, not `z`, for two reasons.  The admissible interval is bounded,
since `E[μ^{−u}]` needs `u < −moment_floor(fields)`, so a root arbitrarily close to
criticality is reached by bisection without an arbitrary upper cut-off on `z`.
And every term is `O(u)` near criticality, where the linear form
`r^{1/z} − (1 − D²/z²)` would difference two numbers both `1 − O(10⁻⁸)`.

`G(0) = 0` identically, and `G'(0)` is minus the numerator of
[`rtfim_delta`](@ref), so on the disordered side `G` leaves zero downwards and a
root is where it returns.
"""
function _rtfim_griffiths_residual(u, m::RandomTFIM)
    return u * log(clean_model(m).J / clean_model(m).h) +
           log_moment(disorder(m, :J), u) +
           log_moment(disorder(m, :h), -u)
end

# Duality interchanges bonds with fields: the WHOLE problem, not just the two
# scales ([IgloiMonthus2005](@cite), below Eq. (4.15)).  With unequal families
# `min(J,h)/max(J,h)` is not that swap and gets the ordered side wrong.
function _rtfim_dual(m::RandomTFIM)
    c = clean_model(m)
    return Disordered(TFIM(; J=c.h, h=c.J); J=disorder(m, :h), h=disorder(m, :J))
end

"""
    _rtfim_solve_u(m::RandomTFIM) -> Union{Float64,Nothing}

Bisect [`_rtfim_griffiths_residual`](@ref) on `0 < u < −moment_floor(fields)`,
returning `nothing` when the residual never turns positive there, which means
the condition has NO root, not that the search gave up.

That case is real and reachable.  For a family bounded away from zero the
residual's slope at large `u` is `ln(J/h) − ln(min μ)`, so a root exists only
while `J·max λ > h·min μ`: the strongest bond must beat the weakest field, or no
region can be locally ordered and there is no Griffiths phase to have an
exponent.  The bracket is therefore checked rather than assumed: assuming it
let the bisection collapse onto its own starting point and return that as an
answer.
"""
function _rtfim_solve_u(m::RandomTFIM)
    fl = moment_floor(disorder(m, :h))
    u_cap = isfinite(fl) ? -fl : Inf
    G(u) = _rtfim_griffiths_residual(u, m)

    hi = isfinite(u_cap) ? u_cap * (1 - 1e-12) : 1.0
    if isfinite(u_cap)
        G(hi) > 0 || return nothing
    else
        steps = 0
        while !(G(hi) > 0)
            hi *= 2
            (steps += 1) > 300 && return nothing
        end
    end

    lo = 0.0
    for _ in 1:300
        mid = 0.5 * (lo + hi)
        (hi - lo) <= 1e-14 * hi && break
        G(mid) > 0 ? (hi = mid) : (lo = mid)
    end
    return 0.5 * (lo + hi)
end

"""
    fetch(m::RandomTFIM, ::DynamicalExponent, ::Infinite) -> Float64

The Griffiths-phase dynamical exponent, exactly: the positive root of
`[(J/h)^{1/z}]_av = 1` ([IgloiMonthus2005](@cite) Eq. (4.15), §4.1.3).  It
varies continuously with the distance from criticality and diverges as `δ → 0`;
near criticality it reduces to the `1/z = 2|δ|` the review quotes with Eq. (4.51).

Throws at `δ = 0`: at the infinite-randomness fixed point the equation has no
root rather than a large one.  Ask for [`ActivatedExponent`](@ref) instead.
"""
function fetch(m::RandomTFIM, ::DynamicalExponent, ::Infinite; kwargs...)
    δ = rtfim_delta(m)
    isfinite(δ) || return error(
        "RandomTFIM: the distance from criticality is $δ, so there is nothing to " *
        "solve. It is not finite when both families have zero spread in ln λ, and " *
        "`BinaryDisorder(1.0)` is a deterministic coupling. A chain with no " *
        "disorder has no Griffiths phase.",
    )
    iszero(δ) && return error(
        "RandomTFIM at [ln J]_av == [ln h]_av is the infinite-randomness critical " *
        "point, where no finite dynamical exponent exists: the gap closes as " *
        "ln(1/Δ) ~ ξ^ψ, so [(J/h)^{1/z}]_av = 1 has no root. Ask for " *
        "`ActivatedExponent()` (= 1/2). A finite z exists on either side, and " *
        "diverges as criticality is approached.",
    )
    u = _rtfim_solve_u(δ > 0 ? m : _rtfim_dual(m))
    u === nothing && return error(
        "RandomTFIM: [(J/h)^{1/z}]_av = 1 has no root for these distributions, so " *
        "there is no Griffiths dynamical exponent to return. This happens when the " *
        "strongest bond cannot beat the weakest field (for two-valued couplings, " *
        "when J·max λ ≤ h·min μ): no region can be locally ordered, so there are no " *
        "rare regions and no Griffiths phase. δ = $δ says how far from criticality " *
        "the chain is, and says nothing about this.",
    )
    return 1 / u
end

"""
    fetch(::RandomTFIM, ::ActivatedExponent, ::Infinite) -> Rational{Int}

`ψ = 1/2` at the critical point, the exponent of `ln t_r ∼ ξ^ψ`
([IgloiMonthus2005](@cite) Eq. (4.13), §4.1.3).

Refused away from criticality, where the chain is in a Griffiths phase with a
finite [`DynamicalExponent`](@ref) instead. Returning `1/2` there would name the
exponent of a fixed point the model is not at.
"""
function fetch(m::RandomTFIM, ::ActivatedExponent, ::Infinite; kwargs...)
    iszero(rtfim_delta(m)) || return error(
        "RandomTFIM: ψ is the exponent of the infinite-randomness fixed point, " *
        "reached only at [ln J]_av == [ln h]_av; this chain has δ = " *
        "$(rtfim_delta(m)) and is in a Griffiths phase, where `DynamicalExponent()` " *
        "is finite and is the operative exponent.",
    )
    return 1 // 2
end

"""
    fetch(::RandomTFIM, ::UniversalityClass, ::Infinite) -> Universality{:IsingSDRG}

At criticality the chain flows to the infinite-randomness fixed point, whose
exact exponents are `fetch(Universality(:IsingSDRG), CriticalExponents(); d=2)`.
Refused off criticality, where the chain is not critical at all.
"""
function fetch(m::RandomTFIM, ::UniversalityClass, ::Infinite; kwargs...)
    iszero(rtfim_delta(m)) || return error(
        "RandomTFIM: only the critical chain has a universality class; this one " *
        "has δ = $(rtfim_delta(m)) and is off criticality.",
    )
    return Universality(:IsingSDRG)
end
