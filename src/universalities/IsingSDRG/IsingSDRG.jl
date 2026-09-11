# ─────────────────────────────────────────────────────────────────────────────
# IsingSDRG — the exact exponent table of the infinite-randomness fixed point
# of the 1D random transverse-field Ising chain (RTFIC).
#
# The effective central charge and the refusal of the Calabrese–Cardy closed
# forms live in `universalities/behaviour/CardyEntanglement.jl`, next to the
# formulas they are about.  Only the exponents are here.
#
# References: Fisher, [FisherDS1995](@cite); the collected set is Table 1d of
# Iglói–Monthus, [IgloiMonthus2005](@cite).
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Universality{:IsingSDRG}, ::CriticalExponents; d::Int=2) -> NamedTuple

Exact exponents of the infinite-randomness fixed point of the 1D random
transverse-field Ising chain.

| field | value | |
| --- | --- | --- |
| `β` | `(3−√5)/2` | bulk magnetisation |
| `x_m` | `(3−√5)/4` | its scaling dimension, `β/ν` |
| `β_s` | `1` | surface magnetisation |
| `x_m_s` | `1/2` | its scaling dimension |
| `ν` | `2` | average correlation length, `ξ ∼ |δ|^{−ν}` |
| `ν_typ` | `1` | typical correlation length — a different, smaller exponent |
| `ψ` | `1/2` | activated dynamic scaling, `ln(1/Δ) ∼ ξ^ψ` |
| `φ` | `(1+√5)/2` | cluster moment, `μ ∼ |ln Ω|^φ` — the golden mean |

Rational where the value is rational; `β`, `x_m` and `φ` are irrational.

`α`, `γ`, `δ` and `η` are **not** here. Table 1d does not carry them, and the
usual scaling relations cannot supply them: those assume power-law dynamic
scaling, which is exactly what this fixed point does not have.

!!! warning "Two different `d`"
    `d = 2` is QAtlas's Euclidean convention (2D classical ≙ 1+1D quantum),
    matching `fetch(Universality(:IsingSDRG), CentralCharge(); d=2)`. The `d`
    in AbstractQAtlas's infinite-randomness relations (`ActivatedMomentGrowth`,
    `GriffithsSusceptibility`, `GriffithsSpecificHeat`) is the **spatial**
    dimension of the chain, `d = 1`. Passing `2` there gives a wrong answer
    that does not look wrong.
"""
function fetch(::Universality{:IsingSDRG}, ::CriticalExponents; d::Int=2, kwargs...)
    d == 2 || return error(
        "Universality{:IsingSDRG} CriticalExponents: only d=2 (1+1D) supported; got d=$d.",
    )
    return (
        β=(3 - sqrt(5)) / 2,
        x_m=(3 - sqrt(5)) / 4,
        β_s=1 // 1,
        x_m_s=1 // 2,
        ν=2 // 1,
        ν_typ=1 // 1,
        ψ=1 // 2,
        φ=(1 + sqrt(5)) / 2,
    )
end

function fetch(m::Universality{:IsingSDRG}, q::CriticalExponents, ::Infinite; kwargs...)
    return fetch(m, q; kwargs...)
end

"""
    fetch(::Universality{:IsingSDRG}, ::ActivatedExponent) -> Rational{Int}

`ψ = 1/2`, from the random walk the strong-disorder RG maps the chain onto.
The exponent of `ln(1/Δ) ∼ ξ^ψ`, which is what this fixed point has instead of
a dynamical exponent — see the [`DynamicalExponent`](@ref) method below.
"""
fetch(::Universality{:IsingSDRG}, ::ActivatedExponent; kwargs...) = 1 // 2

function fetch(m::Universality{:IsingSDRG}, q::ActivatedExponent, ::Infinite; kwargs...)
    return fetch(m, q; kwargs...)
end

"""
    fetch(::Universality{:IsingSDRG}, ::DynamicalExponent)

Always throws. No finite `z` exists at an infinite-randomness fixed point: the
gap closes as `ln(1/Δ) ∼ ξ^ψ`, so `−d(ln Δ)/d(ln ξ)` grows without bound rather
than settling on a value.

This is a refusal, not a gap in the table. Returning some large `z` would be
wrong in a way that reads as a measurement.
"""
function fetch(::Universality{:IsingSDRG}, ::DynamicalExponent; kwargs...)
    return error(
        "Universality{:IsingSDRG}: no finite dynamical exponent exists at an " *
        "infinite-randomness fixed point — the gap closes as ln(1/Δ) ~ ξ^ψ, not " *
        "Δ ~ ξ^{-z}, so -d(ln Δ)/d(ln ξ) grows without bound. Ask for " *
        "`fetch(Universality(:IsingSDRG), ActivatedExponent())` instead (= 1//2). " *
        "A finite z DOES exist off criticality, in the Griffiths phase, where it " *
        "varies continuously with the distance from the transition.",
    )
end
