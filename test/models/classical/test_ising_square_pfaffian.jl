# ─────────────────────────────────────────────────────────────────────────────
# Lattice 非依存の standalone test for IsingSquare partition function.
#
# Tests special values that can be verified without building a lattice:
#   β = 0 → Z = 2^N  (all configurations equally weighted)
#   J = 0 → Z = 2^N  (no interactions, Z is independent of β)
#   positivity + finiteness checks
#
# No Lattice2D dependency.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "IsingSquare special values (Lattice 非依存)" begin
    @testset "β = 0 → Z = 2^N for various (Lx, Ly)" begin
        for (Lx, Ly) in [(2, 2), (2, 3), (3, 3), (3, 4), (4, 4)]
            N = Lx * Ly
            Z = QAtlas.fetch(IsingSquare(), PartitionFunction(); Lx=Lx, Ly=Ly, β=0.0)
            @test Z ≈ 2.0^N atol = 1e-8
        end
    end

    @testset "J = 0 → Z = 2^N (no interactions, any β)" begin
        for β in [0.1, 0.5, 1.0, 5.0]
            for (Lx, Ly) in [(2, 2), (2, 3), (3, 3)]
                N = Lx * Ly
                Z = QAtlas.fetch(
                    IsingSquare(), PartitionFunction(); Lx=Lx, Ly=Ly, β=β, J=0.0
                )
                @test Z ≈ 2.0^N atol = 1e-8
            end
        end
    end

    @testset "Z > 0 and finite for standard parameters" begin
        for (Lx, Ly) in [(2, 2), (2, 3), (3, 3)]
            for β in [0.0, 0.2, 0.44, 1.0, 2.0]
                Z = QAtlas.fetch(IsingSquare(), PartitionFunction(); Lx=Lx, Ly=Ly, β=β)
                @test Z > 0
                @test isfinite(Z)
            end
        end
    end

    @testset "Z monotone increasing in β (ferromagnetic J > 0)" begin
        # For J > 0 the ground state has E_gs < 0, so Z ≥ 2·exp(β|E_gs|)
        # grows with β. The Boltzmann sum is dominated by the two
        # all-up/all-down configurations as β → ∞.
        for (Lx, Ly) in [(2, 2), (3, 3)]
            betas = [0.0, 0.1, 0.5, 1.0, 2.0]
            Zs = [
                QAtlas.fetch(IsingSquare(), PartitionFunction(); Lx=Lx, Ly=Ly, β=β) for
                β in betas
            ]
            for k in 1:(length(Zs) - 1)
                @test Zs[k] <= Zs[k + 1]
            end
        end
    end

    @testset "Symmetry: swap Lx ↔ Ly gives same Z" begin
        # Z is symmetric under transposing the lattice (same model).
        for β in [0.0, 0.3, 1.0]
            Z_23 = QAtlas.fetch(IsingSquare(), PartitionFunction(); Lx=2, Ly=3, β=β)
            Z_32 = QAtlas.fetch(IsingSquare(), PartitionFunction(); Lx=3, Ly=2, β=β)
            @test Z_23 ≈ Z_32 rtol = 1e-10
        end
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "IsingSquare PartitionFunction — verification cards" begin
    # Brute-force exact_partition (independent of the transfer-matrix
    # / Pfaffian path) on small PBC lattices.
    for (L, β, J) in ((3, 0.3, 1.0), (3, 0.5, 1.0), (4, 0.2, 1.0))
        verify(
            IsingSquare(; Lx=L, Ly=L, J=J),
            PartitionFunction(),
            PBC(0);
            route=:ed_finite_size,
            fetch_kw=(; β=β, Lx=L, Ly=L, J=J),
            independent=exact_partition(L, L, J, β),
            agree_within=1e-6,
            refs=["Brute-force Σ_σ exp(-βE) over 2^(L²) configs (square_pbc_bond_pairs)"],
        )
    end

    # β = 0: every configuration has weight 1 => Z = 2^(Lx·Ly)
    verify(
        IsingSquare(; Lx=3, Ly=3, J=1.0),
        PartitionFunction(),
        PBC(0);
        route=:limiting_case,
        fetch_kw=(; β=0.0, Lx=3, Ly=3, J=1.0),
        independent=2.0^9,
        agree_within=1e-6,
        refs=["β = 0: all 2^N configs weight 1 => Z = 2^N"],
    )
end
