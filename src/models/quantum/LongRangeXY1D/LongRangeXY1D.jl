# ─────────────────────────────────────────────────────────────────────────────
# LongRangeXY1D — 1D XY chain with power-law-decaying interactions in a
# transverse field.
#
# Hamiltonian:
#
#     H = -J Σ_{i<j} (σᵢˣ σⱼˣ + σᵢʸ σⱼʸ) / |i-j|^α  -  h Σᵢ σᵢᶻ,
#
#     J > 0,   h ∈ ℝ,   α > 0.
#
# Like LongRangeIsing1D, this interpolates between the LMG-like
# all-to-all limit (α = 0) and the nearest-neighbour XX chain in a
# transverse field (α → ∞).  At α = ∞ the model is equivalent to the
# isotropic XX chain in a transverse field, exactly solvable via the
# Jordan-Wigner transformation (Lieb-Schultz-Mattis 1961; Pfeuty 1970).
# The mass gap then reads
#
#     Δ = 2 · max(0, |h| - 2J),
#
# with a gapless XX (Luttinger-liquid) phase for |h| < 2J and a
# polarised gapped phase for |h| > 2J.
#
# For finite α, the model is non-integrable and requires numerical
# treatment (DMRG / variational MPS).  Such regimes are deferred to
# Phase 2 of this file, following the LongRangeIsing1D phase-1 / phase-2
# convention established in PR #322 and the XY-chain power-law literature
# (Koffel-Lewenstein-Tagliacozzo 2012; Maghrebi-Gong-Gorshkov 2017).
#
# References:
#   - E. Lieb, T. Schultz, D. Mattis, "Two soluble models of an
#     antiferromagnetic chain", [LiebSchultzMattis1961](@cite).
#   - P. Pfeuty, "The one-dimensional Ising model with a transverse
#     field", [Pfeuty1970](@cite).
#   - T. Koffel, M. Lewenstein, L. Tagliacozzo, Phys. Rev. Lett. 109,
#     267203 (2012).
#   - M. F. Maghrebi, Z.-X. Gong, A. V. Gorshkov, Phys. Rev. Lett. 119,
#     023001 (2017).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Pauli σ (this file)
#   Observable:  Spin S = σ/2  (QAtlas-wide spin convention; see docs/src/conventions.md)

"""
    LongRangeXY1D(; J::Real = 1.0, h::Real = 0.0, α::Real = Inf)
        <: AbstractQAtlasModel

1D XY chain with power-law-decaying interactions in a transverse field:

    H = -J Σ_{i<j} (σᵢˣ σⱼˣ + σᵢʸ σⱼʸ) / |i-j|^α  -  h Σᵢ σᵢᶻ.

Constraints: `J > 0`, `α > 0`, `h ∈ ℝ`.

Quantities registered (Phase 1):

| Quantity              | BC         | Method                                  |
| --------------------- | ---------- | --------------------------------------- |
| [`MassGap`](@ref)     | `Infinite` | analytic, α = Inf NN XX limit (LSM)     |

Phase 1 exposes only the `α = Inf` nearest-neighbour XX-in-transverse-field
limit, whose mass gap is the closed form `Δ = 2·max(0, |h| - 2J)`
(Lieb-Schultz-Mattis 1961, isotropic XX limit of the XY chain; Pfeuty 1970 for the TFIM-context derivation).  Finite `α` raises
`DomainError` pending the Phase 2 DMRG follow-up (Maghrebi-Gong-Gorshkov
2017 for the XY chain power-law family).

# References

- E. Lieb, T. Schultz, D. Mattis, *Annals Phys.* **16**, 407 (1961).
- P. Pfeuty, *Annals Phys.* **57**, 79 (1970).
- M. F. Maghrebi, Z.-X. Gong, A. V. Gorshkov, *Phys. Rev. Lett.* **119**,
  023001 (2017).
"""
struct LongRangeXY1D <: AbstractQAtlasModel
    J::Float64
    h::Float64
    α::Float64
    function LongRangeXY1D(J::Real, h::Real, α::Real)
        J > 0 || throw(DomainError(J, "LongRangeXY1D requires J > 0; got J = $J."))
        α > 0 || throw(DomainError(α, "LongRangeXY1D requires α > 0; got α = $α."))
        return new(Float64(J), Float64(h), Float64(α))
    end
end
LongRangeXY1D(; J::Real=1.0, h::Real=0.0, α::Real=Inf) = LongRangeXY1D(J, h, α)

# ═══════════════════════════════════════════════════════════════════════════════
# Mass gap — α = ∞ NN XX-in-transverse-field closed form (Lieb-Schultz-Mattis 1961; Pfeuty 1970 TFIM context)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::LongRangeXY1D, ::MassGap, ::Infinite;
          J=m.J, h=m.h, α=m.α, kwargs...) -> Float64

Mass gap of the 1D long-range XY chain in a transverse field at the
α = ∞ nearest-neighbour limit:

    Δ = 2 · max(0, |h| - 2J).

Gapless (XX / Luttinger-liquid) for `|h| < 2J`, critical at `|h| = 2J`,
polarised gapped for `|h| > 2J`.  Finite `α` raises `DomainError` —
Phase 2 (DMRG; Maghrebi-Gong-Gorshkov 2017).

# References

- E. Lieb, T. Schultz, D. Mattis, *Annals Phys.* **16**, 407 (1961).
- P. Pfeuty, *Annals Phys.* **57**, 79 (1970).
"""
function fetch(
    m::LongRangeXY1D,
    ::MassGap,
    ::Infinite;
    J::Real=m.J,
    h::Real=m.h,
    α::Real=m.α,
    kwargs...,
)
    J > 0 || throw(DomainError(J, "LongRangeXY1D MassGap requires J > 0; got J = $J."))
    α > 0 || throw(DomainError(α, "LongRangeXY1D MassGap requires α > 0; got α = $α."))
    if !isinf(α)
        throw(
            DomainError(
                α,
                "LongRangeXY1D MassGap: closed form supported only at α = Inf (nearest-neighbour " *
                "isotropic XX in transverse field). Finite α (Maghrebi-Gong-Gorshkov 2017) deferred " *
                "to Phase 2. Got α = $α.",
            ),
        )
    end
    # α = ∞: nearest-neighbour XX in transverse field, JW closed form (Lieb-Schultz-Mattis 1961; Pfeuty 1970 TFIM context)
    return 2 * max(0.0, abs(h) - 2J)
end
