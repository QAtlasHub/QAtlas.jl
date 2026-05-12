# ─────────────────────────────────────────────────────────────────────────────
# Verification: classical 2D Ising partition function
#
# Cross-validate QAtlas's transfer-matrix result against a brute-force
# 2^N enumeration. The bond list is built locally (see
# `test/util/classical_partition.jl`) so the test is independent of the
# bond-listing convention of any particular `Lattice2D` release — both sides
# use the same PBC sum convention by construction.
#
# Test sizes: 2×2, 2×3, 3×3 (up to 9 sites = 512 configurations).
# β values: sweep from 0 through a value near the infinite-volume Tc
#            (βc = ln(1 + √2)/2 ≈ 0.4407) and into the ordered phase.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

const J_ISING = 1.0

@testset "IsingSquare — transfer-matrix vs brute-force" begin
    for (Lx, Ly) in [(2, 2), (2, 3), (3, 3)]
        @testset "$(Lx)×$(Ly) square PBC" begin
            for β in [0.0, 0.1, 0.2, 0.44, 1.0, 2.0]
                Z_bf = exact_partition(Lx, Ly, J_ISING, β)
                Z_tm = QAtlas.fetch(
                    IsingSquare(), PartitionFunction(); Lx=Lx, Ly=Ly, β=β, J=J_ISING
                )
                @test Z_tm ≈ Z_bf rtol = 1e-10
            end
        end
    end
end
