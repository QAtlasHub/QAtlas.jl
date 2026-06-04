# ─────────────────────────────────────────────────────────────────────────────
# Mean-field (Landau) universality class
#
# Exact for d ≥ d_c (upper critical dimension; d_c = 4 for Ising-like).
# Serves as the baseline reference for all universality classes.
#
# References:
#   L. D. Landau, Zh. Eksp. Teor. Fiz. 7, 19 (1937).
# ─────────────────────────────────────────────────────────────────────────────

"""
    MeanField() <: AbstractQAtlasModel

Mean-field (Landau) universality class.  Exact for `d ≥ d_c` (upper
critical dimension).  Kept as a top-level alias so that existing
`fetch(MeanField(), ...)` callers keep working after the v0.13 API
redesign; the canonical form is `Universality(:MeanField)`.
"""
struct MeanField <: AbstractQAtlasModel end

"""
    fetch(::MeanField, ::CriticalExponents) -> NamedTuple

Mean-field critical exponents (exact, Rational{Int}):
α=0, β=1/2, γ=1, δ=3, ν=1/2, η=0.

Satisfies Rushbrooke, Widom, and Fisher scaling relations exactly.
"""
function fetch(::MeanField, ::CriticalExponents; kwargs...)
    return (α=0 // 1, β=1 // 2, γ=1 // 1, δ=3 // 1, ν=1 // 2, η=0 // 1)
end
