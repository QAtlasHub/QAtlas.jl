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

# ─────────────────────────────────────────────────────────────────────────────
# Migrated from raw @test to verify()-first (PR #449 phase 2 zero-legacy).
# Each (Lx, Ly, β) pair in the brute-force vs transfer-matrix sweep becomes
# an :ed_finite_size verify card; the sweep is the full pre-migration grid
# {(2,2), (2,3), (3,3)} × {0.0, 0.1, 0.2, 0.44, 1.0, 2.0} = 18 pins.
# ─────────────────────────────────────────────────────────────────────────────
@testset "IsingSquare — transfer-matrix vs brute-force (verify cards)" begin
    for (Lx, Ly) in [(2, 2), (2, 3), (3, 3)]
        for β in [0.0, 0.1, 0.2, 0.44, 1.0, 2.0]
            verify(
                IsingSquare(; Lx=Lx, Ly=Ly, J=J_ISING),
                PartitionFunction(),
                PBC(0);
                route=:ed_finite_size,
                fetch_kw=(; β=β, Lx=Lx, Ly=Ly, J=J_ISING),
                independent=exact_partition(Lx, Ly, J_ISING, β),
                agree_within=1e-10,
                at=["Lx=$(Lx)", "Ly=$(Ly)", "β=$(β)"],
                refs=[
                    "Brute-force Σ_σ exp(-βE) over all 2^N configurations (independent enumeration) vs transfer-matrix Z",
                ],
            )
        end
    end
end
