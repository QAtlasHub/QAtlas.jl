# =============================================================================
# TFIM Z-axis Infinite — closed-form layers (m_z, CorrelationLength, ZZ connected vs static)
#
# Split out of test/models/quantum/TFIM/test_TFIM_zaxis.jl (3.4 min on s05).
# Helpers _build_tfim_dense, _op_site, _SZ come from
# test/util/tfim_dense_ed.jl via runtests.jl ambient include.
# =============================================================================

@testset "TFIM Z-axis Infinite observables — closed-form layers" begin

    # ───────────────────────────────────────────────────────────────────────
    # Layer 1: Pfeuty 1970 spontaneous magnetisation.
    # m_z = (1 - (h/J)²)^{1/8}  for h < J,  0 otherwise.
    # ───────────────────────────────────────────────────────────────────────
    @testset "SpontaneousMagnetization & MagnetizationZ closed form" begin
        J = 1.0
        for h in (0.0, 0.3, 0.5, 0.8, 0.999)
            mz_expected = (1 - (h / J)^2)^(1 / 8)
            mz_q = QAtlas.fetch(TFIM(; J=J, h=h), MagnetizationZ(), Infinite())
            mz_s = QAtlas.fetch(TFIM(; J=J, h=h), SpontaneousMagnetization(), Infinite())
            @test mz_q ≈ mz_expected atol=1e-12
            @test mz_s ≈ mz_expected atol=1e-12
        end
        # Critical and disordered phases — m_z = 0 exactly.
        for h in (1.0, 1.0001, 2.0, 5.0)
            @test QAtlas.fetch(TFIM(; J=1.0, h=h), MagnetizationZ(), Infinite()) == 0.0
            @test QAtlas.fetch(
                TFIM(; J=1.0, h=h), SpontaneousMagnetization(), Infinite()
            ) == 0.0
        end
        # Non-unit J: scale invariance under (J, h) → (αJ, αh).
        for α in (0.5, 2.0)
            mz1 = QAtlas.fetch(TFIM(; J=1.0, h=0.6), MagnetizationZ(), Infinite())
            mz2 = QAtlas.fetch(TFIM(; J=α, h=α * 0.6), MagnetizationZ(), Infinite())
            @test mz1 ≈ mz2 atol=1e-12
        end
    end

    # ───────────────────────────────────────────────────────────────────────
    # Layer 2: Correlation length closed form ξ = 1/(2|h-J|).
    # ───────────────────────────────────────────────────────────────────────
    @testset "CorrelationLength = 1/(2|h - J|)" begin
        J = 1.0
        for h in (0.3, 0.7, 1.5, 2.0)
            ξ_expected = 1 / (2 * abs(h - J))
            ξ_q = QAtlas.fetch(TFIM(; J=J, h=h), CorrelationLength(), Infinite())
            @test ξ_q ≈ ξ_expected atol=1e-12
        end
        # Critical point: ξ = ∞.
        @test isinf(QAtlas.fetch(TFIM(; J=1.0, h=1.0), CorrelationLength(), Infinite()))
    end

    # ───────────────────────────────────────────────────────────────────────
    # Layer 3: ZZ connected correlator at OBC equals the static one (Z₂).
    # ───────────────────────────────────────────────────────────────────────
    @testset "ZZCorrelation{:connected} = ZZCorrelation{:static} (Z₂)" begin
        for h in (0.5, 1.0, 1.5)
            for β in (Inf, 2.0, 0.5)
                for r in 1:4
                    v_static = QAtlas.fetch(
                        TFIM(; J=1.0, h=h),
                        ZZCorrelation{:static}(),
                        OBC(; N=10);
                        beta=β,
                        i=2,
                        j=2 + r,
                    )
                    v_conn = QAtlas.fetch(
                        TFIM(; J=1.0, h=h),
                        ZZCorrelation{:connected}(),
                        OBC(; N=10);
                        beta=β,
                        i=2,
                        j=2 + r,
                    )
                    @test v_conn ≈ v_static atol=1e-12
                end
            end
        end
    end
end
