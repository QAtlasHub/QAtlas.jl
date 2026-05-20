# ─────────────────────────────────────────────────────────────────────────────
# PpIp2DSC — 2-D p_x + i p_y chiral superconductor (Read-Green 2000).
#
# Mean-field BCS Hamiltonian for spinless fermions on a 2-D continuum
# (or lattice) with chiral p-wave pairing:
#
#     H = Σ_k [ ξ_k c†_k c_k + (Δ_k c†_k c†_{-k} + h.c.) ],
#     ξ_k = ε_k − μ,           Δ_k = Δ₀ (k_x + i k_y).
#
# Read–Green (2000) showed that the BdG ground state of this
# Hamiltonian is gapped in the bulk and falls into two distinct phases
# as a function of the chemical potential `μ`:
#
#   * Weak-pairing phase  (μ > 0): topologically non-trivial,
#     first Chern number `C = 1`, single chiral Majorana edge mode,
#     `c = 1/2` Ising-like boundary CFT, non-Abelian Ising anyons
#     bound to vortex cores.
#   * Strong-pairing phase (μ < 0): trivial (`C = 0`), no edge modes.
#
# This is the prototypical 2-D fermionic topological superconductor
# and is the candidate state for the ν = 5/2 fractional quantum Hall
# plateau (Moore–Read Pfaffian).
#
# Phase 1 of this entry registers the two parameter-independent
# bulk topological invariants of the weak-pairing phase:
#
#   * CentralCharge        c = 1/2   (chiral Majorana edge CFT)
#   * TopologicalInvariant C = 1     (first Chern number)
#
# Vortex-core Majorana zero-mode quantities, the BdG spectrum on a
# strip / disc, and the strong-pairing phase (`μ ≤ 0`) are tracked as
# Phase 2 and are intentionally excluded here.
#
# References:
#   - N. Read, D. Green, *Phys. Rev. B* **61**, 10267 (2000).
#   - A. Y. Kitaev, *Ann. Phys.* **321**, 2 (2006).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Fermion bilinears c†c
#   Observable:  Fermion (number n = c†c, bilinear ⟨c†_i c_j⟩); derived spin observables follow spin S = σ/2
#   Reference:   docs/src/conventions.md §Fermion convention

"""
    PpIp2DSC(; Δ₀::Real = 1.0, μ::Real = 1.0) <: AbstractQAtlasModel

2-D p_x + i p_y chiral superconductor (Read–Green 2000) for spinless
fermions, with BCS pairing amplitude `Δ₀ > 0` and chemical potential
`μ > 0` (weak-pairing topological phase).

Phase 1 exposes the parameter-independent bulk topological data of
the weak-pairing phase:

| Quantity                            | BC         | Value   |
| ----------------------------------- | ---------- | ------- |
| [`CentralCharge`](@ref)             | `Infinite` | `1/2`   |
| [`TopologicalInvariant`](@ref)      | `Infinite` | `1`     |

The strong-pairing trivial phase (`μ ≤ 0`) is excluded by the
constructor: only the topological branch is exposed in Phase 1.

# References

- N. Read, D. Green, *Phys. Rev. B* **61**, 10267 (2000).
- A. Y. Kitaev, *Ann. Phys.* **321**, 2 (2006).
"""
struct PpIp2DSC <: AbstractQAtlasModel
    Δ₀::Float64
    μ::Float64
    function PpIp2DSC(Δ₀::Real, μ::Real)
        Δ₀ > 0 || throw(DomainError(Δ₀, "PpIp2DSC requires Δ₀ > 0; got Δ₀ = $Δ₀."))
        μ > 0 || throw(
            DomainError(
                μ,
                "PpIp2DSC requires μ > 0 (weak-pairing topological phase); got μ = $μ.",
            ),
        )
        return new(Float64(Δ₀), Float64(μ))
    end
end
PpIp2DSC(; Δ₀::Real=1.0, μ::Real=1.0) = PpIp2DSC(Δ₀, μ)

# ═══════════════════════════════════════════════════════════════════════════════
# Chiral Majorana edge CFT — central charge c = 1/2
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::PpIp2DSC, ::CentralCharge, ::Infinite; kwargs...) -> Rational{Int}

Central charge of the chiral Majorana edge CFT of the 2-D p+ip
weak-pairing superconductor (Read–Green 2000):

    c = 1/2

(a single right-moving Majorana fermion, i.e. the chiral half of the
2-D Ising model). Parameter-independent within the weak-pairing
topological phase (`μ > 0`).

# References

- N. Read, D. Green, *Phys. Rev. B* **61**, 10267 (2000).
- A. Y. Kitaev, *Ann. Phys.* **321**, 2 (2006).
"""
function fetch(::PpIp2DSC, ::CentralCharge, ::Infinite; kwargs...)
    return 1 // 2
end

# ═══════════════════════════════════════════════════════════════════════════════
# First Chern number of the BdG ground state — C = 1
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::PpIp2DSC, ::TopologicalInvariant, ::Infinite; kwargs...) -> Int

First Chern number of the BdG ground-state band of the 2-D p+ip
weak-pairing superconductor (Read–Green 2000):

    C = 1

(in contrast to `C = 0` in the strong-pairing trivial phase at
`μ < 0`). By the bulk–boundary correspondence this drives a single
chiral Majorana edge mode.

# References

- N. Read, D. Green, *Phys. Rev. B* **61**, 10267 (2000).
"""
function fetch(::PpIp2DSC, ::TopologicalInvariant, ::Infinite; kwargs...)
    return 1
end
