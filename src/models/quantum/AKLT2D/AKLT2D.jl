# ─────────────────────────────────────────────────────────────────────────────
# AKLT2D — 2D Affleck-Kennedy-Lieb-Tasaki VBS on the hexagonal (honeycomb) lattice
#
# Hamiltonian (J > 0 antiferromagnetic) with spin S = 3/2 on each vertex:
#
#   H = J Σ_{<i,j>} P₃(Sᵢ + Sⱼ),
#
# where the sum is over nearest-neighbour bonds of the honeycomb lattice
# and P₃ projects onto the total-spin S_{i,j} = 3 subspace of the bond.
# The unique ground state is the Valence Bond Solid (VBS) product of
# singlet covers on the honeycomb edges — a 2D PEPS of bond dimension 2
# and the prototypical 2D AKLT-class symmetry-protected topological (SPT)
# state (Affleck-Kennedy-Lieb-Tasaki 1988; Verstraete-Cirac 2004 PEPS).
#
# The Hamiltonian is frustration-free:  each P₃ ≥ 0, and the VBS state
# is annihilated by every bond projector (each bond contains only two
# spin-3/2 dimers, which cannot couple to total S = 3 in the VBS), so
#
#   e₀ / N = 0     (exact, J-independent)
#
# The bulk gap is conjectured to remain open in the thermodynamic limit,
# supported by rigorous finite-size lower bounds (Pomata-Wei 2020;
# Lemm-Sandvik-Wang 2020) — no strict TD-limit proof is known for the
# bare honeycomb AKLT model. The best numerical lower bound
# (Garcia-Saez-Murg-Wei 2013; Pomata-Wei 2020) is
#
#   Δ_2D-AKLT / J ≈ 0.10   (variational / DMRG; not in Phase 1).
#
# References:
#   I. Affleck, T. Kennedy, E. H. Lieb, H. Tasaki,
#     "Valence bond ground states in isotropic quantum antiferromagnets",
#     Commun. Math. Phys. 115, 477 (1988).
#   F. Verstraete and J. I. Cirac,
#     "Renormalization algorithms for quantum-many body systems in two
#      and higher dimensions", cond-mat/0407066 (2004) — PEPS realisation.
# ─────────────────────────────────────────────────────────────────────────────

struct AKLT2D <: AbstractQAtlasModel
    J::Float64
    function AKLT2D(J::Real)
        J > 0 || throw(DomainError(J, "AKLT2D requires J > 0; got J = $J."))
        return new(Float64(J))
    end
end
AKLT2D(; J::Real=1.0) = AKLT2D(J)

QAtlas.native_energy_granularity(::AKLT2D, ::Infinite) = :per_site

"""
    fetch(::AKLT2D, ::Energy{:per_site}, ::Infinite; J=m.J) -> Float64

Ground-state energy density of the 2D AKLT model on the honeycomb
lattice (Affleck-Kennedy-Lieb-Tasaki 1988):

    e_0 / N = 0

The Hamiltonian is a sum of non-negative spin-3 projectors, and the
valence-bond-solid (VBS) state is annihilated by every term →
frustration-free with exact zero ground-state energy density. The VBS
state is the prototypical 2D AKLT-class SPT and a paradigmatic
non-trivial PEPS (bond dimension 2).

# References

- I. Affleck, T. Kennedy, E. H. Lieb, H. Tasaki, *Commun. Math. Phys.* **115**, 477 (1988).
- F. Verstraete, J. I. Cirac, *cond-mat/0407066* (2004) — PEPS realisation.
"""
function fetch(m::AKLT2D, ::Energy{:per_site}, ::Infinite; J::Real=m.J, kwargs...)
    J > 0 || throw(DomainError(J, "AKLT2D Energy requires J > 0; got J = $J."))
    return 0.0
end
