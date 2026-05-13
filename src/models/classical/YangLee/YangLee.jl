# ─────────────────────────────────────────────────────────────────────────────
# YangLee — non-unitary Virasoro minimal model M(5, 2) describing the
# Lee-Yang edge singularity (Yang-Lee 1952; Cardy 1985).
#
# Central charge
#
#     c(5, 2) = 1 - 6 (5 - 2)² / (5 · 2) = 1 - 54/10 = -22/5
#
# is negative — allowed because (5, 2) lies outside the Friedan-Qiu-Shenker
# unitarity series M(p+1, p) (p ≥ 2).  The Kac table has only two distinct
# primaries:
#
#     h_{1,1} = 0       (identity)
#     h_{1,2} = -1/5    (the famous "negative-dimension" Yang-Lee field)
#
# with the Kac symmetry (r, s) ↔ (p_prime - r, p - s) identifying
# (1, 1) ↔ (1, 4) and (1, 2) ↔ (1, 3).
#
# Yang-Lee (1952) showed that the partition function zeros of the Ising
# model in an imaginary magnetic field lie on a unit circle and pinch
# the positive real axis at an edge singularity.  Fisher (1978) and
# Cardy (1985) identified the universality class of this edge as a
# non-unitary CFT with c = -22/5, realised by M(5, 2).
#
# This Phase-1 entry registers `CentralCharge` and `ConformalWeights`,
# delegating to [`MinimalModel(5, 2)`](@ref).
#
# References:
#   - C. N. Yang and T. D. Lee, Phys. Rev. 87, 404 (1952);
#     T. D. Lee and C. N. Yang, Phys. Rev. 87, 410 (1952).
#   - J. L. Cardy, Phys. Rev. Lett. 54, 1354 (1985).
#   - M. E. Fisher, Phys. Rev. Lett. 40, 1610 (1978).
# ─────────────────────────────────────────────────────────────────────────────

"""
    YangLee() <: AbstractQAtlasModel

Yang-Lee CFT — the non-unitary Virasoro minimal model M(5, 2) describing
the universality class of the Lee-Yang edge singularity (Yang-Lee 1952;
Cardy 1985).

The model has no continuous parameters: M(5, 2) is fixed.

Quantities registered (Phase 1):

| Quantity                       | BC         | Method                         |
| ------------------------------ | ---------- | ------------------------------ |
| [`CentralCharge`](@ref)        | `Infinite` | delegated to MinimalModel(5,2) |
| [`ConformalWeights`](@ref)     | `Infinite` | delegated to MinimalModel(5,2) |

The famous negative-dimension primary `h_{1,2} = -1/5` is the edge
exponent and is what makes the theory non-unitary (c = -22/5 < 0).

# References

- C. N. Yang and T. D. Lee, *Phys. Rev.* **87**, 404 (1952).
- J. L. Cardy, *Phys. Rev. Lett.* **54**, 1354 (1985).
"""
struct YangLee <: AbstractQAtlasModel end

# ═══════════════════════════════════════════════════════════════════════════════
# Central charge via MinimalModel(5, 2) delegation
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::YangLee, ::CentralCharge, ::Infinite; kwargs...) -> Rational{Int}

Central charge of the Yang-Lee CFT (non-unitary M(5, 2)):

    c = -22/5

delegated to [`MinimalModel(5, 2)`](@ref).  Returned as an exact
`Rational{Int}` (`-22//5`).

# References

- J. L. Cardy, *Phys. Rev. Lett.* **54**, 1354 (1985).
- C. N. Yang and T. D. Lee, *Phys. Rev.* **87**, 404 (1952).
"""
function fetch(::YangLee, ::CentralCharge, ::Infinite; kwargs...)
    return QAtlas.fetch(QAtlas.MinimalModel(5, 2), CentralCharge())
end

# ═══════════════════════════════════════════════════════════════════════════════
# Conformal weights (Kac formula) via MinimalModel(5, 2) delegation
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::YangLee, ::ConformalWeights, ::Infinite; r::Integer=1, s::Integer=1, kwargs...)
        -> Rational{Int}

Kac-formula conformal weight `h_{r,s}` of the Yang-Lee CFT, delegated
to [`MinimalModel(5, 2)`](@ref).

For M(p, p_prime) = M(5, 2) the Kac-table index range is
`r ∈ [1, p_prime - 1] = [1, 1]` (only `r = 1`) and
`s ∈ [1, p - 1] = [1, 4]`, with the Kac symmetry
`(r, s) ↔ (p_prime - r, p - s)` identifying

    (1, 1) ↔ (1, 4)   with   h = 0           (identity)
    (1, 2) ↔ (1, 3)   with   h = -1/5         (Yang-Lee primary)

so the theory has only two distinct primaries.

Out-of-range `(r, s)` throws `DomainError` (inherited from
`MinimalModel`).

# References

- J. L. Cardy, *Phys. Rev. Lett.* **54**, 1354 (1985).
"""
function fetch(
    ::YangLee, ::ConformalWeights, ::Infinite; r::Integer=1, s::Integer=1, kwargs...
)
    return QAtlas.fetch(QAtlas.MinimalModel(5, 2), ConformalWeights(); r=r, s=s)
end
