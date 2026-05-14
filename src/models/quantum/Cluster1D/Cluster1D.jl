# ─────────────────────────────────────────────────────────────────────────────
# Cluster1D — 1D Z₂×Z₂ symmetry-protected topological (SPT) cluster Hamiltonian.
#
# Hamiltonian (Briegel-Raussendorf 2001 / cluster-state form):
#
#     H = -J Σ_i K_i,    K_i = σ^z_{i-1} σ^x_i σ^z_{i+1},   J > 0.
#
# The stabilisers `K_i` mutually commute (each pair shares two sites with a
# z/x mismatch, picking up two anticommutations → commute).  The Hilbert
# space decomposes into joint eigenspaces of {K_i}; the unique ground state
# is the +1 simultaneous eigenstate ("cluster state"), the canonical
# 1D Z₂×Z₂ bosonic SPT phase (Verresen-Moessner-Pollmann 2017; Else et al
# 2012).  Symmetry generators are P_even = ∏_{i even} σ^x_i and
# P_odd = ∏_{i odd} σ^x_i (the two Z₂ factors).
#
# Phase-1 closed-form observables:
#
#   * E_0 / N = -J            (each K_i contributes -J in the ground state)
#   * Δ      = 2J             (single-K_i flip: -J → +J)
#
# The cluster state is the seed resource of measurement-based quantum
# computation (one-way QC, Raussendorf-Briegel 2001 PRL).
#
# References:
#
#   - H. J. Briegel, R. Raussendorf, Phys. Rev. Lett. 86, 910 (2001) —
#     persistent entanglement, cluster state.
#   - R. Raussendorf, H. J. Briegel, Phys. Rev. Lett. 86, 5188 (2001) —
#     one-way quantum computer.
#   - D. V. Else, S. D. Bartlett, A. C. Doherty, New J. Phys. 14, 113016
#     (2012) — symmetry protection of the cluster Hamiltonian.
#   - R. Verresen, R. Moessner, F. Pollmann, Phys. Rev. B 96, 165124 (2017)
#     — relation to two coupled Kitaev chains.
# ─────────────────────────────────────────────────────────────────────────────

"""
    Cluster1D(; J::Real = 1.0) <: AbstractQAtlasModel

1D Z₂×Z₂ symmetry-protected topological cluster Hamiltonian
`H = -J Σ_i σ^z_{i-1} σ^x_i σ^z_{i+1}`.  Ground state is the cluster
state (Briegel-Raussendorf 2001).

Quantities registered (Phase 1):

| Quantity                       | BC         | Method                              |
| ------------------------------ | ---------- | ----------------------------------- |
| [`Energy`](@ref) (`:per_site`) | `Infinite` | analytic (stabiliser ground state)  |
| [`MassGap`](@ref)              | `Infinite` | analytic (single-stabiliser flip)   |

# References

- H. J. Briegel, R. Raussendorf, *Phys. Rev. Lett.* **86**, 910 (2001).
- D. V. Else, S. D. Bartlett, A. C. Doherty, *New J. Phys.* **14**, 113016 (2012).
- R. Verresen, R. Moessner, F. Pollmann, *Phys. Rev. B* **96**, 165124 (2017).
"""
struct Cluster1D <: AbstractQAtlasModel
    J::Float64
    function Cluster1D(J::Real)
        J > 0 || throw(DomainError(J, "Cluster1D requires J > 0; got J = $J."))
        return new(Float64(J))
    end
end
Cluster1D(; J::Real=1.0) = Cluster1D(J)

# Native energy granularity at the thermodynamic limit
QAtlas.native_energy_granularity(::Cluster1D, ::Infinite) = :per_site

# ═══════════════════════════════════════════════════════════════════════════════
# Ground-state energy density
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::Cluster1D, ::Energy{:per_site}, ::Infinite; J=m.J) -> Float64

Ground-state energy density `E_0 / N = -J`.  Every stabiliser
`K_i = σ^z_{i-1} σ^x_i σ^z_{i+1}` contributes -J in the cluster-state
ground state.
"""
function fetch(m::Cluster1D, ::Energy{:per_site}, ::Infinite; J::Real=m.J, kwargs...)
    J > 0 || throw(DomainError(J, "Cluster1D Energy requires J > 0; got J = $J."))
    return -float(J)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Mass gap
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::Cluster1D, ::MassGap, ::Infinite; J=m.J) -> Float64

Spectral gap `Δ = 2J`.  Flipping a single stabiliser eigenvalue
`K_i: +1 → -1` raises the energy by 2J; all other states differ by
integer multiples of 2J.
"""
function fetch(m::Cluster1D, ::MassGap, ::Infinite; J::Real=m.J, kwargs...)
    J > 0 || throw(DomainError(J, "Cluster1D MassGap requires J > 0; got J = $J."))
    return 2.0 * float(J)
end
