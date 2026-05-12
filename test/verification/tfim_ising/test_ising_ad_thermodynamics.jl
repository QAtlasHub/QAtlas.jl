# ─────────────────────────────────────────────────────────────────────────────
# Verification: thermodynamic quantities from the partition function via AD
#
# Check that macroscopic observables (internal energy, heat capacity, average
# bond energy) obtained by automatic differentiation of the transfer-matrix
# log Z agree with direct Boltzmann ensemble averages computed by brute-force
# enumeration of all 2^N spin configurations.
#
# Identities being verified:
#
#   ⟨E⟩   = -∂(log Z)/∂β
#   C_v   =  β² · ∂²(log Z)/∂β²  =  β² · (⟨E²⟩ − ⟨E⟩²)
#   ⟨ΣJσσ⟩ = (1/β) · ∂(log Z)/∂J
#
# Every equality must hold to machine precision because Z is the same
# mathematical object on both sides — the AD path exercises QAtlas's
# transfer-matrix Z, and the direct path enumerates configurations under
# the same PBC bond list (see `test/util/classical_partition.jl`).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, ForwardDiff, Test

# log Z as a generic function of (β, J) so ForwardDiff.Dual propagates.
function _log_Z_tm(Lx::Int, Ly::Int, β, J)
    return log(QAtlas.fetch(IsingSquare(), PartitionFunction(); Lx=Lx, Ly=Ly, β=β, J=J))
end

"""
    _direct_moments(Lx, Ly, J, β) -> (⟨E⟩, ⟨E²⟩, ⟨Σσσ⟩)

Directly compute the ensemble averages needed for the thermodynamic
identities by brute-force enumeration. Uses the same PBC bond list as
[`exact_partition`](@ref) so the comparison against the transfer-matrix
observables is exact by construction.
"""
function _direct_moments(Lx::Int, Ly::Int, J::Float64, β::Float64)
    N = Lx * Ly
    bond_pairs = square_pbc_bond_pairs(Lx, Ly)
    Z = 0.0
    E1 = 0.0
    E2 = 0.0
    S1 = 0.0  # ⟨Σ σ_i σ_j⟩ (no J factor)
    for σ_idx in 0:(2 ^ N - 1)
        σ = Int[((σ_idx >> k) & 1) == 1 ? 1 : -1 for k in 0:(N - 1)]
        bond_sum = sum(σ[a] * σ[b] for (a, b) in bond_pairs)
        E = -J * bond_sum
        w = exp(-β * E)
        Z += w
        E1 += E * w
        E2 += E^2 * w
        S1 += bond_sum * w
    end
    return (E1 / Z, E2 / Z, S1 / Z)
end

const J_ISING = 1.0

@testset "IsingSquare — thermodynamics from log Z (ForwardDiff)" begin
    for (Lx, Ly) in [(2, 2), (2, 3), (3, 3)]
        @testset "$(Lx)×$(Ly) square PBC" begin
            for β in [0.1, 0.3, 0.5, 1.0]
                # --- AD from transfer matrix ---
                #   E = -∂(log Z)/∂β
                E_ad = -ForwardDiff.derivative(βv -> _log_Z_tm(Lx, Ly, βv, J_ISING), β)

                #   C_v = β² · ∂²(log Z)/∂β²
                d2logZ_dβ2 = ForwardDiff.derivative(
                    βo -> ForwardDiff.derivative(βi -> _log_Z_tm(Lx, Ly, βi, J_ISING), βo),
                    β,
                )
                Cv_ad = β^2 * d2logZ_dβ2

                #   ⟨Σ σσ⟩ = (1/β) · ∂(log Z)/∂J
                bond_sum_ad =
                    (1 / β) *
                    ForwardDiff.derivative(Jv -> _log_Z_tm(Lx, Ly, β, Jv), J_ISING)

                # --- Direct ensemble averages ---
                E_direct, E2_direct, bond_sum_direct = _direct_moments(Lx, Ly, J_ISING, β)
                Cv_direct = β^2 * (E2_direct - E_direct^2)

                @test E_ad ≈ E_direct rtol = 1e-10
                @test Cv_ad ≈ Cv_direct rtol = 1e-10
                @test bond_sum_ad ≈ bond_sum_direct rtol = 1e-10
            end
        end
    end
end
