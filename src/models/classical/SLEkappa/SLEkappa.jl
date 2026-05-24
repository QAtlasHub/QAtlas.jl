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
# This entry registers two analytic quantities:
#   - `CentralCharge` (Phase 1) via the SLE-CFT correspondence,
#   - `FractalDimension` (Phase 2) via Beffara's theorem.
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

Quantities registered:

| Quantity                       | BC         | Method                                |
| ------------------------------ | ---------- | ------------------------------------- |
| [`CentralCharge`](@ref)        | `Infinite` | analytic (SLE-CFT correspondence)     |
| [`FractalDimension`](@ref)     | `Infinite` | analytic (Beffara 2008)               |

# References

- O. Schramm, *Israel J. Math.* **118**, 221 (2000).
- V. Beffara, *Annals Probab.* **36**, 1421 (2008).
- M. Bauer, D. Bernard, *Phys. Rep.* **432**, 115 (2006).
"""
# CONVENTION
#   Hamiltonian: see file-header description above
#   Observable:  per src/core/quantities.jl (matches the dispatch tag)
#   Reference:   docs/src/conventions.md (project-wide convention policy)
#   STATUS:      backfilled by PR (audit gate); per-field domain content
#                left to a follow-up - see issue tracker for the model-specific
#                Hamiltonian sign / observable normalisation.

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

# ═══════════════════════════════════════════════════════════════════════════════
# Hausdorff dimension d_H(κ) = min(2, 1 + κ/8)  (Beffara 2008, Phase 2)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::SLEkappa, ::FractalDimension, ::Infinite; κ=m.κ) -> Float64

Hausdorff dimension of the SLE_κ random curve (Beffara 2008):

    d_H(κ) = min(2, 1 + κ/8),    κ > 0.

The cap at `d_H = 2` for `κ ≥ 8` reflects the space-filling regime
(the SLE_κ curve becomes space-filling at κ = 8 and stays so beyond,
so its Hausdorff dimension cannot exceed the ambient plane).
Canonical fixed points evaluate to: `d_H(2)=5/4`, `d_H(8/3)=4/3`,
`d_H(3)=11/8`, `d_H(4)=3/2`, `d_H(6)=7/4`, `d_H(8)=2`.

`κ ≤ 0` raises `DomainError` (the SLE process is undefined at and
below zero diffusivity).

# References

- V. Beffara, *Annals Probab.* **36**, 1421 (2008).
"""
function fetch(m::SLEkappa, ::FractalDimension, ::Infinite; κ::Real=m.κ, kwargs...)
    κ > 0 || throw(DomainError(κ, "SLEkappa FractalDimension requires κ > 0; got κ = $κ."))
    return min(2.0, 1.0 + κ / 8.0)
end
