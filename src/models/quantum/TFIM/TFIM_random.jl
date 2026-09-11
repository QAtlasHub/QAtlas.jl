# ─────────────────────────────────────────────────────────────────────────────
# Random transverse-field Ising chain — the Griffiths dynamical exponent.
#
# `TFIM` is the clean chain and carries two scalars, which is why the exact
# Griffiths result could not live there: it is a condition on the DISTRIBUTIONS
# of the couplings.  That condition, [(J/h)^{1/z}]_av = 1, is itself
# distribution-free; only its evaluation is per-family, which is what
# `DisorderFamily` separates.
#
# Reference: Iglói–Monthus, [IgloiMonthus2005](@cite).  Equation and section
# numbers are those of the arXiv version, `cond-mat/0502448`.
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Pauli σ, as in `TFIM.jl` — H = -Σ Jᵢ σᶻσᶻ - Σ hᵢ σˣ
#   Disorder:    Jᵢ = J·λᵢ and hᵢ = h·μᵢ with λ, μ i.i.d. dimensionless and
#                J, h the SCALES, not the couplings; see `DisorderFamily`.

"""
    DisorderFamily

A distribution of dimensionless couplings `λ > 0`, carrying the three things the
Griffiths condition needs: [`log_moment`](@ref), [`mean_log`](@ref) and
[`var_log`](@ref).  A concrete family also gives [`moment_floor`](@ref), the
exponent below which `E[λ^s]` stops existing.
"""
abstract type DisorderFamily end
export DisorderFamily

"""
    log_moment(f::DisorderFamily, s::Real) -> Float64

`ln E[λ^s]`, computed in logs so the Griffiths residual stays cancellation-free
near criticality, where every term is `O(1/z)`.
"""
function log_moment end
export log_moment

"""    mean_log(f::DisorderFamily) -> Float64 — `E[ln λ]`."""
function mean_log end
export mean_log

"""    var_log(f::DisorderFamily) -> Float64 — `var[ln λ]`."""
function var_log end
export var_log

"""
    moment_floor(f::DisorderFamily) -> Float64

`inf{s : E[λ^s] < ∞}`, or `-Inf` where every moment exists.  The field family's
floor is what bounds `z` from below: the condition needs `E[μ^{−1/z}]`.
"""
function moment_floor end
export moment_floor

# A family missing one of the four is otherwise a raw MethodError at whichever
# call path happens to reach it first — and three of the four `fetch` routes do
# not touch `moment_floor` at all, so a partial family can pass them. Same shape
# as `derivative(::AbstractDiffBackend, ...)` in `core/derivative.jl`: name the
# contract rather than let dispatch report the symptom.
for _f in (:log_moment, :mean_log, :var_log, :moment_floor)
    @eval function $(_f)(f::DisorderFamily, args...; kwargs...)
        return error(
            "QAtlas.$($(QuoteNode(_f))): not defined for $(nameof(typeof(f))). A " *
            "`DisorderFamily` must implement all four of `log_moment`, `mean_log`, " *
            "`var_log` and `moment_floor` — `mean_log` and `var_log` are the first " *
            "two derivatives of `log_moment` at s = 0, so they must agree with it.",
        )
    end
end

"""
    PowerLawDisorder(D) <: DisorderFamily

`P(λ) = D⁻¹ λ^{−1+1/D}` on `(0, 1]`, the review's generic disorder-strength
family ([IgloiMonthus2005](@cite) §A.1): `D² = var(ln λ)`, and `D = 1` is
uniform on `[0, 1]`.  `E[λ^s] = 1/(1 + D s)`, so moments below `s = −1/D` do not
exist and `z > D` for a chain whose fields are drawn from it.
"""
struct PowerLawDisorder <: DisorderFamily
    D::Float64
    function PowerLawDisorder(D::Real)
        D > 0 || throw(ArgumentError("PowerLawDisorder: D must be > 0; got $D"))
        return new(Float64(D))
    end
end
export PowerLawDisorder

log_moment(f::PowerLawDisorder, s::Real) = -log1p(f.D * s)
mean_log(f::PowerLawDisorder) = -f.D
var_log(f::PowerLawDisorder) = f.D^2
moment_floor(f::PowerLawDisorder) = -1 / f.D

"""
    BinaryDisorder(κ) <: DisorderFamily

Two couplings, `λ ∈ {1, κ}` with equal probability, `0 < κ ≤ 1` — the other
standard RTFIM choice.  Bounded away from zero, so `E[λ^s]` exists for EVERY `s`
and `z` has no lower bound: the `z > D` of [`PowerLawDisorder`](@ref) is a
property of that family, not of the Griffiths condition.
"""
struct BinaryDisorder <: DisorderFamily
    κ::Float64
    function BinaryDisorder(κ::Real)
        0 < κ <= 1 || throw(ArgumentError("BinaryDisorder: need 0 < κ ≤ 1; got $κ"))
        return new(Float64(κ))
    end
end
export BinaryDisorder

# log((1 + κ^s)/2), guarded for the s → ±∞ tails the bisection walks through.
function log_moment(f::BinaryDisorder, s::Real)
    x = s * log(f.κ)
    return (x > 700 ? x : log1p(exp(x))) - log(2)
end
mean_log(f::BinaryDisorder) = log(f.κ) / 2
var_log(f::BinaryDisorder) = log(f.κ)^2 / 4
moment_floor(::BinaryDisorder) = -Inf

"""
    RandomTFIM(; J = 1.0, h = 1.0, D = 1.0)
    RandomTFIM(J, h, bonds::DisorderFamily, fields::DisorderFamily)

The 1D random transverse-field Ising chain,

    H = -Σ_i J_i σᶻ_i σᶻ_{i+1} - Σ_i h_i σˣ_i

with `J_i = J·λ_i` and `h_i = h·μ_i` for `λ ~ bonds`, `μ ~ fields`.  The keyword
form is the symmetric [`PowerLawDisorder`](@ref) case on both.

Criticality is `[ln J]_av = [ln h]_av`, which for equal families is `J == h` but
in general is not.  Distance from it is the standard control parameter

    δ = ([ln h]_av − [ln J]_av) / (var[ln h] + var[ln J]),

positive in the disordered (field-dominated) phase; see [`rtfim_delta`](@ref).

| Quantity | BC | Coverage |
| --- | --- | --- |
| [`DynamicalExponent`](@ref) | `Infinite` | exact, off criticality; throws at `δ = 0` |
| [`ActivatedExponent`](@ref) | `Infinite` | `1/2`, at criticality only |
| [`UniversalityClass`](@ref) | `Infinite` | `:IsingSDRG`, at criticality only |
"""
struct RandomTFIM{B<:DisorderFamily,F<:DisorderFamily} <: AbstractQAtlasModel
    J::Float64
    h::Float64
    bonds::B
    fields::F
    function RandomTFIM(J::Real, h::Real, bonds::B, fields::F) where {B,F}
        J > 0 || throw(ArgumentError("RandomTFIM: J must be > 0; got $J"))
        h > 0 || throw(ArgumentError("RandomTFIM: h must be > 0; got $h"))
        return new{B,F}(Float64(J), Float64(h), bonds, fields)
    end
end
function RandomTFIM(; J::Real=1.0, h::Real=1.0, D::Real=1.0)
    return RandomTFIM(J, h, PowerLawDisorder(D), PowerLawDisorder(D))
end
export RandomTFIM

"""
    rtfim_delta(m::RandomTFIM) -> Float64

`δ = ([ln h]_av − [ln J]_av) / (var[ln h] + var[ln J])`.  Zero exactly at the
infinite-randomness critical point; positive in the disordered phase.
"""
function rtfim_delta(m::RandomTFIM)
    num = (log(m.h) + mean_log(m.fields)) - (log(m.J) + mean_log(m.bonds))
    return num / (var_log(m.fields) + var_log(m.bonds))
end
export rtfim_delta

"""
    _rtfim_griffiths_residual(u, m)

`ln [(J/h)^{1/z}]_av` at `u = 1/z` ([`IgloiMonthus2005`](@cite) Eq. (4.15),
§4.1.3), which the condition sets to zero.

Solved in `u`, not `z`, for two reasons.  The admissible interval is bounded —
`E[μ^{−u}]` needs `u < −moment_floor(fields)` — so a root arbitrarily close to
criticality is reached by bisection without an arbitrary upper cut-off on `z`.
And every term is `O(u)` near criticality, where the linear form
`r^{1/z} − (1 − D²/z²)` would difference two numbers both `1 − O(10⁻⁸)`.

`G(0) = 0` identically, and `G'(0)` is minus the numerator of
[`rtfim_delta`](@ref), so on the disordered side `G` leaves zero downwards and a
root is where it returns.
"""
function _rtfim_griffiths_residual(u, m::RandomTFIM)
    return u * log(m.J / m.h) + log_moment(m.bonds, u) + log_moment(m.fields, -u)
end

# Duality interchanges bonds with fields — the WHOLE problem, not just the two
# scales ([IgloiMonthus2005](@cite), below Eq. (4.15)).  With unequal families
# `min(J,h)/max(J,h)` is not that swap and gets the ordered side wrong.
_rtfim_dual(m::RandomTFIM) = RandomTFIM(m.h, m.J, m.fields, m.bonds)

"""
    _rtfim_solve_u(m::RandomTFIM) -> Union{Float64,Nothing}

Bisect [`_rtfim_griffiths_residual`](@ref) on `0 < u < −moment_floor(fields)`,
returning `nothing` when the residual never turns positive there — which means
the condition has NO root, not that the search gave up.

That case is real and reachable.  For a family bounded away from zero the
residual's slope at large `u` is `ln(J/h) − ln(min μ)`, so a root exists only
while `J·max λ > h·min μ`: the strongest bond must beat the weakest field, or no
region can be locally ordered and there is no Griffiths phase to have an
exponent.  The bracket is therefore checked rather than assumed — assuming it
let the bisection collapse onto its own starting point and return that as an
answer.
"""
function _rtfim_solve_u(m::RandomTFIM)
    fl = moment_floor(m.fields)
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
        "solve. It is not finite when both families have zero spread in ln λ — " *
        "`BinaryDisorder(1.0)` is a deterministic coupling — and a chain with no " *
        "disorder has no Griffiths phase.",
    )
    iszero(δ) && return error(
        "RandomTFIM at [ln J]_av == [ln h]_av is the infinite-randomness critical " *
        "point, where no finite dynamical exponent exists — the gap closes as " *
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
finite [`DynamicalExponent`](@ref) instead — returning `1/2` there would name the
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
