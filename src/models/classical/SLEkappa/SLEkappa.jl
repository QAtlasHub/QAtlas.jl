# ─────────────────────────────────────────────────────────────────────────────
# SLEkappa — Schramm-Loewner Evolution SLE_κ (random conformal curves).
#
# SLE_κ (Schramm 2000) is the unique conformally-invariant one-parameter
# family of random planar curves; κ ∈ [0, ∞) sets all geometric and
# CFT-related universal data.  The two scalars exposed by this Phase-1
# entry are
#
#   * Hausdorff dimension (Beffara 2008):
#         d_H(SLE_κ) = min(2, 1 + κ/8),    κ ∈ [0, 8],  capped at 2 for κ ≥ 8.
#
#   * SLE-CFT correspondence (Bauer-Bernard 2006; Cardy 2005):
#         c(κ) = (3 κ − 8)(6 − κ) / (2 κ),
#     symmetric under the duality κ ↔ 16/κ that swaps SLE_κ with its
#     outer-boundary curve SLE_{16/κ}.  The map yields the canonical
#     central charges
#
#         κ = 2  →  c = -2     (loop-erased random walk)
#         κ = 8/3 → c = 0      (self-avoiding walk)
#         κ = 3  →  c = 1/2    (Ising spin cluster boundary)
#         κ = 4  →  c = 1      (GFF level lines)
#         κ = 6  →  c = 0      (percolation cluster boundary)
#         κ = 8  →  c = -2     (uniform-spanning-tree Peano curve).
#
# This Phase-1 entry registers `CentralCharge` only.  A
# `FractalDimension` quantity (or generic `HausdorffDimension`)
# requires a new core struct in `src/core/quantities.jl` and is
# tracked as a follow-up phase.
#
# References:
#   - O. Schramm, Israel J. Math. 118, 221 (2000).
#   - V. Beffara, Annals Probab. 36, 1421 (2008) — d_H = 1 + κ/8.
#   - M. Bauer, D. Bernard, Phys. Rep. 432, 115 (2006) — SLE-CFT.
#   - J. Cardy, Annals Phys. 318, 81 (2005) — physicist review.
# ─────────────────────────────────────────────────────────────────────────────

"""
    SLEkappa(; κ::Real = 6.0) <: AbstractQAtlasModel

Schramm-Loewner Evolution SLE_κ (random planar conformally-invariant
curve) at parameter `κ > 0`.  Default `κ = 6` is the percolation
cluster-boundary fixed point — the original Schramm prediction that
launched the field.

Quantities registered (Phase 1):

| Quantity                       | BC         | Method                              |
| ------------------------------ | ---------- | ----------------------------------- |
| [`CentralCharge`](@ref)        | `Infinite` | analytic (SLE-CFT correspondence)   |

A `FractalDimension` quantity (Beffara 2008 `d_H = min(2, 1 + κ/8)`)
requires a new core struct and is tracked as a follow-up phase.

# References

- O. Schramm, *Israel J. Math.* **118**, 221 (2000).
- V. Beffara, *Annals Probab.* **36**, 1421 (2008).
- M. Bauer, D. Bernard, *Phys. Rep.* **432**, 115 (2006).
"""
struct SLEkappa <: AbstractQAtlasModel
    κ::Float64
end
SLEkappa(; κ::Real=6.0) = SLEkappa(Float64(κ))

# ═══════════════════════════════════════════════════════════════════════════════
# Central charge from the SLE-CFT correspondence
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::SLEkappa, ::CentralCharge, ::Infinite; κ=m.κ) -> Float64

Central charge of the CFT dual to SLE_κ (Bauer-Bernard 2006):

    c(κ) = (3 κ − 8) (6 − κ) / (2 κ),    κ > 0.

The map is symmetric under κ ↔ 16/κ, encoding the SLE duality between
a chordal curve and its outer-boundary curve.  Canonical fixed points
return the expected rational values: `c(2) = -2`, `c(8/3) = 0`,
`c(3) = 1/2`, `c(4) = 1`, `c(6) = 0`, `c(8) = -2`.

`κ ≤ 0` raises `DomainError` (the SLE process is undefined at and
below zero diffusivity).

# References

- M. Bauer, D. Bernard, *Phys. Rep.* **432**, 115 (2006).
- J. Cardy, *Annals Phys.* **318**, 81 (2005).
"""
function fetch(m::SLEkappa, ::CentralCharge, ::Infinite; κ::Real=m.κ, kwargs...)
    κ > 0 || throw(DomainError(κ, "SLEkappa CentralCharge requires κ > 0; got κ = $κ."))
    return (3κ - 8) * (6 - κ) / (2κ)
end
