# ─────────────────────────────────────────────────────────────────────────────
# Shastry–Sutherland — spin-½ AF on the Shastry–Sutherland lattice.
#
# Hamiltonian (Miyahara–Ueda 1999 convention):
#
#   H = J  Σ_{⟨i,j⟩} Sᵢ · Sⱼ      (square-lattice nearest-neighbour bonds)
#     + J' Σ_{⟨⟨i,j⟩⟩_d} Sᵢ · Sⱼ  (alternating-plaquette diagonal "dimer" bonds),
#                                  S = 1/2.
#
# At α = J / J' ≤ α_c ≈ 0.675 (Koga–Kawakami 2000) the ground state is
# *exactly* the product of nearest-neighbour singlets on the diagonal
# (J') bonds — the SrCu₂(BO₃)₂ realisation (Kageyama 1999).  On the
# dimer state each diagonal bond contributes ⟨S·S⟩_singlet = -3/4 and
# the square-lattice J bonds give zero by orthogonality between
# adjacent dimers (each spin is in exactly one singlet, the orthogonal
# component cancels in the matrix element), so the size-independent
# closed-form ground-state energy density is
#
#     E₀/N = -3 J' / 8.
#
# The triplet (triplon) gap is finite in the dimer phase (Δ ≈ J' for
# α → 0, decreasing with α; closed form known only as a series in α
# (Miyahara–Ueda 1999)).  Magnetisation plateaux at 1/8, 1/4, 1/3
# (Kageyama 1999 / Miyahara–Ueda 2003 / DMRG) require additional
# infrastructure and are tracked as a follow-up phase.
#
# This is the 2-D analog of the 1-D Majumdar–Ghosh chain (#166):
# both lock onto an exact nearest-neighbour-singlet dimer ground state
# at a single point in coupling space.
#
# References:
#   - B. S. Shastry, B. Sutherland, Physica B+[ShastrySutherland1981](@cite).
#   - S. Miyahara, K. Ueda, [MiyaharaUeda1999](@cite).
#   - A. Koga, N. Kawakami, [KogaKawakami2000](@cite).
#   - H. Kageyama et al., [Kageyama1999](@cite) — material.
# ─────────────────────────────────────────────────────────────────────────────

# Exact-dimer phase boundary alpha_c = J/J' below which the singlet-product
# state is the exact ground state.  Numerical value 0.675(5) from Koga-Kawakami
# (PRL 84, 4461, 2000); later series-expansion / iPEPS work (Corboz-Mila 2013;
# Boos-Toldin 2019) refines this within the same band.  Revisit if a tighter
# bound becomes the community consensus.

# CONVENTION
#   Hamiltonian: Spin S (this file)
#   Observable:  Spin S         (QAtlas-wide spin convention; see docs/src/conventions.md)

const _SS_DIMER_PHASE_ALPHA_CRIT = 0.675

"""
    ShastrySutherland(; J::Real = 0.0, Jp::Real = 1.0) <: AbstractQAtlasModel

Spin-½ Shastry–Sutherland model on the Shastry–Sutherland lattice
(corner-sharing-orthogonal-dimer geometry; realised in SrCu₂(BO₃)₂):

    H = J Σ_{⟨i,j⟩} Sᵢ · Sⱼ + J' Σ_{⟨⟨i,j⟩⟩_d} Sᵢ · Sⱼ.

The model has two named couplings: `J` for the square-lattice
nearest-neighbour bonds and `Jp` (J') for the alternating-plaquette
diagonal "dimer" bonds.  The default `J = 0, Jp = 1` is the trivial
isolated-dimer limit and lies safely inside the exact-dimer phase.

In the parameter window `α = J / Jp ≤ α_c ≈ 0.675` (Koga–Kawakami
2000) the ground state is the *exact* product of nearest-neighbour
singlets on the diagonal (`Jp`) bonds — the 2-D analog of the
Majumdar–Ghosh dimer state.

Quantities registered:

| Quantity                       | BC         | Method                |
| ------------------------------ | ---------- | --------------------- |
| [`Energy`](@ref) (`:per_site`) | `Infinite` | exact dimer GS        |

The triplon excitation spectrum, magnetisation plateaux (1/8, 1/4,
1/3) and large-α intermediate / plaquette phases require quantities
beyond a scalar ground-state energy and are tracked as a follow-up
phase.

# References

- B. S. Shastry, B. Sutherland, *Physica B+C* **108**, 1069 (1981).
- S. Miyahara, K. Ueda, *Phys. Rev. Lett.* **82**, 3701 (1999).
- A. Koga, N. Kawakami, *Phys. Rev. Lett.* **84**, 4461 (2000).
"""
struct ShastrySutherland <: AbstractQAtlasModel
    J::Float64
    Jp::Float64
end
ShastrySutherland(; J::Real=0.0, Jp::Real=1.0) = ShastrySutherland(Float64(J), Float64(Jp))

# ═══════════════════════════════════════════════════════════════════════════════
# Ground-state energy per site (exact dimer phase)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::ShastrySutherland, ::Energy{:per_site}, ::Infinite;
          J=m.J, Jp=m.Jp) -> Float64

Ground-state energy per site of the spin-½ Shastry–Sutherland model
in the exact-dimer phase `α = J / J' ≤ α_c ≈ 0.675`:

    E₀ / N = -3 J' / 8.

The result is *size-independent* and equal across the entire
exact-dimer parameter window — the J nearest-neighbour bonds
contribute zero on the singlet-product ground state.

`Jp > 0` is required (anti-ferromagnetic dimer bond) and `J / Jp`
must lie at or below the Koga–Kawakami α_c ≈ 0.675; otherwise a
`DomainError` is thrown.

# References

- B. S. Shastry, B. Sutherland, *Physica B+C* **108**, 1069 (1981).
- A. Koga, N. Kawakami, *Phys. Rev. Lett.* **84**, 4461 (2000).
"""
function fetch(
    m::ShastrySutherland,
    ::Energy{:per_site},
    ::Infinite;
    J::Real=m.J,
    Jp::Real=m.Jp,
    kwargs...,
)
    Jp > 0 || throw(
        DomainError(
            Jp,
            "ShastrySutherland Energy(:per_site) requires Jp > 0 (anti-ferromagnetic dimer bond); got Jp = $Jp.",
        ),
    )
    α = J / Jp
    α ≤ _SS_DIMER_PHASE_ALPHA_CRIT || throw(
        DomainError(
            α,
            "ShastrySutherland Energy(:per_site): exact dimer-phase closed form is valid only for J/J' ≤ α_c ≈ $(_SS_DIMER_PHASE_ALPHA_CRIT) (Koga–Kawakami 2000); got α = $α.",
        ),
    )
    return -3 * Jp / 8
end
