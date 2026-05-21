# =============================================================================
# TFIM Z-axis Infinite — N_proxy pass-through (SusceptibilityZZ + ZZStructureFactor, heavy N=80 Pfaffian)
#
# Split out of test/models/quantum/TFIM/test_TFIM_zaxis.jl (3.4 min on s05).
# Helpers _build_tfim_dense, _op_site, _SZ come from
# test/util/tfim_dense_ed.jl via runtests.jl ambient include.
# =============================================================================

@testset "TFIM Z-axis Infinite observables — N_proxy pass-through" begin

    # ───────────────────────────────────────────────────────────────────────
    # Layer 4: SusceptibilityZZ Infinite uses the OBC large-N proxy.
    # The default N_proxy = 80 should match an explicit OBC(80) call exactly.
    # ───────────────────────────────────────────────────────────────────────
    @testset "SusceptibilityZZ Infinite is the OBC N_proxy proxy" begin
        χ_inf = QAtlas.fetch(TFIM(; J=1.0, h=0.5), SusceptibilityZZ(), Infinite(); beta=1.0)
        χ_obc = QAtlas.fetch(
            TFIM(; J=1.0, h=0.5), SusceptibilityZZ(), OBC(; N=80); beta=1.0
        )
        @test χ_inf ≈ χ_obc atol = 1e-12
        # Custom N_proxy: disordered phase, smaller chain.
        χ_inf_40 = QAtlas.fetch(
            TFIM(; J=1.0, h=1.5), SusceptibilityZZ(), Infinite(); beta=1.0, N_proxy=40
        )
        χ_obc_40 = QAtlas.fetch(
            TFIM(; J=1.0, h=1.5), SusceptibilityZZ(), OBC(; N=40); beta=1.0
        )
        @test χ_inf_40 ≈ χ_obc_40 atol = 1e-12
    end

    # ───────────────────────────────────────────────────────────────────────
    # Layer 5: ZZStructureFactor Infinite reduces to OBC large-N proxy.
    # ───────────────────────────────────────────────────────────────────────
    @testset "ZZStructureFactor Infinite is the OBC N_proxy proxy" begin
        S_inf = QAtlas.fetch(
            TFIM(; J=1.0, h=0.7), ZZStructureFactor(), Infinite(); beta=1.0, q=π / 2
        )
        S_obc = QAtlas.fetch(
            TFIM(; J=1.0, h=0.7), ZZStructureFactor(), OBC(; N=80); beta=1.0, q=π / 2
        )
        @test S_inf ≈ S_obc atol=1e-12
    end
end
