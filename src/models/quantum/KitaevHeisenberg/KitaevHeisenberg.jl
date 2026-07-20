# ─────────────────────────────────────────────────────────────────────────────
# KitaevHeisenberg — K-J-Γ honeycomb model (α-RuCl₃ family).
#
# Hamiltonian (Rau-Lee-Kee 2014 convention):
#
#     H = Σ_{⟨ij⟩_γ} [ K Sᵢ^γ Sⱼ^γ
#                    + J Sᵢ · Sⱼ
#                    + Γ (Sᵢ^α Sⱼ^β + Sᵢ^β Sⱼ^α) ],
#
# with γ ∈ {x, y, z} labelling the three honeycomb bond directions and
# (α, β) the two off-diagonal axes complementary to γ.
#
# The K-only limit (J = Γ = 0) reduces to Kitaev's celebrated
# exactly-solvable honeycomb model (Kitaev 2006), already exposed by
# QAtlas as [`KitaevHoneycomb`](@ref).  Switching on J (Heisenberg)
# or Γ (off-diagonal symmetric) immediately destroys the
# Z₂-flux integrability and leaves the model in the realm of DMRG /
# variational Monte Carlo / numerical ED.
#
# Phase 1 of this file therefore only exposes the K-only delegation:
# any (J, Γ) ≠ (0, 0) triple raises DomainError.  Future phases will
# add the Rau-Lee-Kee perturbative phase diagram (ferromagnetic /
# zigzag / Néel / Kitaev spin liquid) and the half-quantized thermal
# Hall conductance κ_xy / T = π/12 (in units of k_B²/ℏ) of the
# field-induced gapped Kitaev spin liquid (Kasahara 2018).
#
# References:
#   - A. Kitaev, "Anyons in an exactly solved model and beyond",
#     [Kitaev2006](@cite).
#   - G. Jackeli, G. Khaliullin, [JackeliKhaliullin2009](@cite).
#   - J. G. Rau, E. K.-H. Lee, H.-Y. Kee, Phys. Rev. Lett. 112,
#     077204 (2014).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Pauli σ (this file)
#   Observable:  Spin S = σ/2  (QAtlas-wide spin convention; see docs/src/conventions.md)

"""
    KitaevHeisenberg(; K::Real = 1.0, J::Real = 0.0, Γ::Real = 0.0)
        <: AbstractQAtlasModel

K-J-Γ honeycomb model (α-RuCl₃ family).  Three independent exchanges
on the three honeycomb bond directions:

    H = Σ_{⟨ij⟩_γ} [ K Sᵢ^γ Sⱼ^γ + J Sᵢ · Sⱼ + Γ (Sᵢ^α Sⱼ^β + Sᵢ^β Sⱼ^α) ].

Quantities registered (Phase 1):

| Quantity                       | BC         | Method                              |
| ------------------------------ | ---------- | ----------------------------------- |
| [`MassGap`](@ref)              | `Infinite` | delegated to KitaevHoneycomb (K=K)  |

Phase 1 exposes only the **K-only limit** `J = Γ = 0`, delegated to
the existing exactly-solvable [`KitaevHoneycomb`](@ref) entry.  Any
non-zero `J` or `Γ` raises `DomainError` pointing to the
DMRG / ED phase-2 follow-up.

# References

- A. Kitaev, *Annals Phys.* **321**, 2 (2006).
- G. Jackeli, G. Khaliullin, *Phys. Rev. Lett.* **102**, 017205 (2009).
- J. G. Rau, E. K.-H. Lee, H.-Y. Kee, *Phys. Rev. Lett.* **112**, 077204 (2014).
"""
struct KitaevHeisenberg <: AbstractQAtlasModel
    K::Float64
    J::Float64
    Γ::Float64
end
function KitaevHeisenberg(; K::Real=1.0, J::Real=0.0, Γ::Real=0.0)
    return KitaevHeisenberg(Float64(K), Float64(J), Float64(Γ))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Mass gap — K-only delegation to KitaevHoneycomb
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::KitaevHeisenberg, ::MassGap, ::Infinite;
          K=m.K, J=m.J, Γ=m.Γ, kwargs...) -> Float64

Single-particle gap of the K-J-Γ honeycomb model at the **K-only
point** `J = Γ = 0`.  Internally constructs the isotropic
[`KitaevHoneycomb`](@ref)`(Kx = K, Ky = K, Kz = K)` and forwards.
The Kitaev gapless A/B/C phase therefore returns `Δ = 0` at
isotropic |K|.

`J ≠ 0` or `Γ ≠ 0` raises `DomainError` — Phase 2.

# References

- A. Kitaev, *Annals Phys.* **321**, 2 (2006).
"""
function fetch(
    m::KitaevHeisenberg,
    ::MassGap,
    ::Infinite;
    K::Real=m.K,
    J::Real=m.J,
    Γ::Real=m.Γ,
    kwargs...,
)
    iszero(J) || throw(
        DomainError(
            J,
            "KitaevHeisenberg MassGap currently only handles the K-only limit J = 0 (delegated to KitaevHoneycomb); J ≠ 0 perturbations require DMRG/ED — Phase 2.",
        ),
    )
    iszero(Γ) || throw(
        DomainError(
            Γ,
            "KitaevHeisenberg MassGap currently only handles the Γ = 0 limit (delegated to KitaevHoneycomb); off-diagonal Γ ≠ 0 perturbations require DMRG/ED — Phase 2.",
        ),
    )
    return QAtlas.fetch(QAtlas.KitaevHoneycomb(; Kx=K, Ky=K, Kz=K), MassGap(), Infinite())
end
