# =============================================================================
# Tests for TFIM Z-axis (longitudinal) Infinite-system observables
# (TFIM_zaxis.jl): MagnetizationZ, SpontaneousMagnetization, CorrelationLength,
# SusceptibilityZZ, ZZStructureFactor (all at Infinite()).
#
# Layers:
#   1. Pfeuty (1970) closed-form spontaneous magnetisation
#         m_z = (1 - (h/J)²)^{1/8}    h < J
#               0                     h ≥ J
#   2. Correlation length ξ = 1/(2|h-J|)  (∞ at h = J)
#   3. ZZCorrelation{:connected} reduces to {:static} since ⟨σᶻ⟩ = 0 by Z₂
#   4. SusceptibilityZZ Infinite is the OBC large-N proxy and converges
#      to the OBC value at the chosen N_proxy.
#   5. ZZStructureFactor Infinite likewise matches the OBC large-N proxy.
# =============================================================================

@testset "TFIM Z-axis Infinite observables" begin

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
    # The TFIM ground state and any thermal state have ⟨σᶻ_i⟩ = 0 by the
    # Z₂ symmetry σᶻ → −σᶻ of the Hamiltonian, so connected = static.
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

    # ───────────────────────────────────────────────────────────────────────
    # Layer 4: SusceptibilityZZ Infinite uses the OBC large-N proxy.
    # The default N_proxy = 80 should match an explicit OBC(80) call exactly.
    # ───────────────────────────────────────────────────────────────────────
    @testset "SusceptibilityZZ Infinite is the OBC N_proxy proxy" begin
        # One representative (h, β) at the default N_proxy and one at a
        # smaller N_proxy is enough — the "Infinite ↔ OBC" identity is a
        # pass-through, not a physical claim, so over-sweeping (h, β) just
        # multiplies the cost.  Default proxy: ordered phase, finite β.
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
    # Single (q, β) point — same pass-through identity as Layer 4.
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
