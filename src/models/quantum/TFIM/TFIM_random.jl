# ─────────────────────────────────────────────────────────────────────────────
# Random transverse-field Ising chain — the Griffiths-phase dynamical exponent.
#
# `TFIM` is the clean chain and carries two scalars, which is why the exact
# Griffiths result could not live there: it is a condition on the DISTRIBUTIONS
# of the couplings, not on their values.
#
# Reference: Iglói–Monthus, [IgloiMonthus2005](@cite).  Equation and section
# numbers are those of the arXiv version, `cond-mat/0502448`.
# ─────────────────────────────────────────────────────────────────────────────

"""
    RandomTFIM(; J = 1.0, h = 1.0, D = 1.0) <: AbstractQAtlasModel

The 1D random transverse-field Ising chain,

    H = -Σ_i J_i σᶻ_i σᶻ_{i+1} - Σ_i h_i σˣ_i

with `J_i = J·λ_i` and `h_i = h·μ_i`, where `λ` and `μ` are i.i.d. from the
disorder-strength family `P(λ) = D⁻¹ λ^{−1+1/D}` on `[0, 1]`
([`IgloiMonthus2005`](@cite) Eq. (A.1), §A.1).  `D` is the strength: `D² =
var(ln λ)`, and `D = 1` is the uniform distribution on `[0, 1]`.

`J` and `h` are the upper cut-offs, not the couplings — the mean log couplings
are `ln J − D` and `ln h − D`, so the chain is critical exactly at `J == h`,
whatever `D`.

Distance from criticality is the standard control parameter

    δ = ([ln h]_av − [ln J]_av) / (var[ln h] + var[ln J]) = ln(h/J) / (2D²),

positive in the disordered (field-dominated) phase.  Both signs give the same
`z` by duality ([`IgloiMonthus2005`](@cite), below Eq. (4.15)).

| Quantity | BC | Coverage |
| --- | --- | --- |
| [`DynamicalExponent`](@ref) | `Infinite` | exact, off criticality; throws at `δ = 0` |
| [`ActivatedExponent`](@ref) | `Infinite` | `1/2`, at criticality only |
| [`UniversalityClass`](@ref) | `Infinite` | `:IsingSDRG`, at criticality only |
"""
struct RandomTFIM <: AbstractQAtlasModel
    J::Float64
    h::Float64
    D::Float64
    function RandomTFIM(J::Real, h::Real, D::Real)
        J > 0 || throw(ArgumentError("RandomTFIM: J must be > 0; got $J"))
        h > 0 || throw(ArgumentError("RandomTFIM: h must be > 0; got $h"))
        D > 0 || throw(ArgumentError("RandomTFIM: D must be > 0; got $D"))
        return new(Float64(J), Float64(h), Float64(D))
    end
end
RandomTFIM(; J::Real=1.0, h::Real=1.0, D::Real=1.0) = RandomTFIM(J, h, D)
export RandomTFIM

"""
    rtfim_delta(m::RandomTFIM) -> Float64

Distance from criticality, `δ = ln(h/J) / (2D²)`.  Zero exactly at the
infinite-randomness critical point; positive in the disordered phase.
"""
rtfim_delta(m::RandomTFIM) = log(m.h / m.J) / (2 * m.D^2)
export rtfim_delta

# The Griffiths condition [(J/h)^{1/z}]_av = 1 ([IgloiMonthus2005](@cite)
# Eq. (4.15), §4.1.3) on the family above.  Both factors are elementary:
#
#     E[λ^{ 1/z}] = 1 / (1 + D/z)          E[μ^{−1/z}] = 1 / (1 − D/z),  z > D
#
# so with r = J/h the condition collapses to  r^{1/z} = 1 − D²/z².  The second
# expectation is what fixes the domain: below z = D the field distribution has
# no 1/z-th inverse moment and the condition has no meaning, not merely no root.
#
# Written with `expm1` rather than as `r^(1/z) - (1 - D^2/z^2)`.  Near criticality
# both of those are 1 − O(10⁻⁸) and their difference loses eight digits to
# cancellation — in the regime where z diverges, which is the regime this oracle
# exists for.  `expm1(ln(r)/z) + D²/z²` is the same root with both terms small.
_rtfim_griffiths_residual(z, r, D) = expm1(log(r) / z) + D^2 / z^2

# Bisection rather than a solver dependency: the residual is monotone in z on
# (D, ∞) for r < 1, and the bracket is known analytically — it is positive as
# z → D⁺ (where the right-hand side vanishes) and negative as z → ∞ (where
# 2|δ|D²/z beats D²/z²).
function _rtfim_solve_z(r::Float64, D::Float64)
    lo, hi = D * (1 + 1e-12), D
    while _rtfim_griffiths_residual(hi, r, D) > 0
        hi *= 2
        hi > 1e12 && return Inf                     # indistinguishable from critical
    end
    for _ in 1:200
        mid = 0.5 * (lo + hi)
        (hi - lo) <= 1e-12 * max(1.0, mid) && return mid
        _rtfim_griffiths_residual(mid, r, D) > 0 ? (lo = mid) : (hi = mid)
    end
    return 0.5 * (lo + hi)
end

"""
    fetch(m::RandomTFIM, ::DynamicalExponent, ::Infinite) -> Float64

The Griffiths-phase dynamical exponent, exactly.

`z` is the positive root of `[(J/h)^{1/z}]_av = 1`
([`IgloiMonthus2005`](@cite) Eq. (4.15), §4.1.3), which on this model's disorder
family is

    (J/h)^{1/z} = 1 − D²/z²,   z > D.

It varies continuously with the distance from criticality and diverges as
`δ → 0`; near criticality it reduces to the closed form `1/z = 2|δ|` the review
quotes with Eq. (4.51).

Throws at `δ = 0`: at the infinite-randomness fixed point no finite `z` exists,
and the equation has no root there rather than a large one.  Ask for
[`ActivatedExponent`](@ref) instead.
"""
function fetch(m::RandomTFIM, ::DynamicalExponent, ::Infinite; kwargs...)
    δ = rtfim_delta(m)
    iszero(δ) && return error(
        "RandomTFIM at J == h is the infinite-randomness critical point, where no " *
        "finite dynamical exponent exists — the gap closes as ln(1/Δ) ~ ξ^ψ, so " *
        "[(J/h)^{1/z}]_av = 1 has no root. Ask for `ActivatedExponent()` (= 1/2). " *
        "A finite z exists on either side of J == h, and diverges as it is approached.",
    )
    # Duality: interchanging h and J gives the ordered-phase value, so the ratio
    # is taken smaller-over-larger and |δ| is what z depends on.
    r = min(m.J, m.h) / max(m.J, m.h)
    return _rtfim_solve_z(r, m.D)
end

"""
    fetch(::RandomTFIM, ::ActivatedExponent, ::Infinite) -> Rational{Int}

`ψ = 1/2` at the critical point, the exponent of `ln t_r ∼ ξ^ψ`
([`IgloiMonthus2005`](@cite) Eq. (4.13), §4.1.3).

Refused away from `J == h`: `ψ` describes the infinite-randomness fixed point,
and off criticality the chain is in a Griffiths phase with a finite
[`DynamicalExponent`](@ref) instead.  Returning `1/2` there would name the
exponent of a fixed point the model is not at.
"""
function fetch(m::RandomTFIM, ::ActivatedExponent, ::Infinite; kwargs...)
    iszero(rtfim_delta(m)) || return error(
        "RandomTFIM: ψ is the exponent of the infinite-randomness fixed point, " *
        "reached only at J == h; this chain has δ = $(rtfim_delta(m)) and is in a " *
        "Griffiths phase, where `DynamicalExponent()` is finite and is the " *
        "operative exponent.",
    )
    return 1 // 2
end

"""
    fetch(::RandomTFIM, ::UniversalityClass, ::Infinite) -> Universality{:IsingSDRG}

At `J == h` the chain flows to the infinite-randomness fixed point, whose exact
exponents are `fetch(Universality(:IsingSDRG), CriticalExponents(); d=2)`.
Refused off criticality, where the chain is not critical at all.
"""
function fetch(m::RandomTFIM, ::UniversalityClass, ::Infinite; kwargs...)
    iszero(rtfim_delta(m)) || return error(
        "RandomTFIM: only the critical chain (J == h) has a universality class; " *
        "this one has δ = $(rtfim_delta(m)) and is off criticality.",
    )
    return Universality(:IsingSDRG)
end
