# =============================================================================
# Standalone tests for Heisenberg1D infinite-volume spinon kinematics defined
# in `src/models/quantum/Heisenberg/Heisenberg_spinon.jl` (issue #154 phase 1).
#
# Coverage:
#   * Single-spinon dispersion ε(k) = (π J / 2) |sin k| at k = 0, π/2, π
#   * 2-spinon lower edge ε_L(q) coincides with the single-spinon dispersion
#   * 2-spinon upper edge ε_U(q) = π J |sin(q/2)|
#   * Boundary values at the gapless point q = π (lower = 0, upper = π J)
#   * Müller-ansatz S^{zz}(q, ω) is 0 outside [ε_L, ε_U] and finite inside
#   * S^{zz}(q, ω) → ∞ as ω → ε_L⁺ (integrable square-root singularity)
#   * Linear J-scaling of all three energy quantities
# =============================================================================

using QAtlas, Test

@testset "Heisenberg1D spinon kinematics (Phase 1, closed-form)" begin
    model = Heisenberg1D()

    # ─────────────────────────────────────────────────────────────────────
    # Single-spinon dispersion
    # ─────────────────────────────────────────────────────────────────────
    @testset "single-spinon dispersion ε(k) = (π J / 2) |sin k|" begin
        @test heisenberg_spinon_dispersion(model, 0.0) ≈ 0.0 atol = 1e-14
        @test heisenberg_spinon_dispersion(model, π / 2) ≈ π / 2 atol = 1e-14
        @test heisenberg_spinon_dispersion(model, π) ≈ 0.0 atol = 1e-13

        # |sin| symmetry: ε(π − k) = ε(k)
        for k in (0.1, 0.7, 1.3)
            @test heisenberg_spinon_dispersion(model, k) ≈
                heisenberg_spinon_dispersion(model, π - k) atol = 1e-14
        end
    end

    # ─────────────────────────────────────────────────────────────────────
    # 2-spinon continuum edges (des Cloizeaux–Pearson 1962)
    # ─────────────────────────────────────────────────────────────────────
    @testset "2-spinon lower edge ε_L(q) ≡ ε(q)" begin
        for q in range(0.0, π; length=17)
            @test heisenberg_two_spinon_lower_edge(model, q) ≈
                heisenberg_spinon_dispersion(model, q) atol = 1e-14
        end
    end

    @testset "2-spinon upper edge ε_U(q) = π J |sin(q/2)|" begin
        @test heisenberg_two_spinon_upper_edge(model, 0.0) ≈ 0.0 atol = 1e-14
        @test heisenberg_two_spinon_upper_edge(model, π / 2) ≈ π * sin(π / 4) atol = 1e-14
        @test heisenberg_two_spinon_upper_edge(model, π) ≈ π atol = 1e-13   # gapless top
    end

    @testset "boundary at q = π: lower = 0, upper = π J" begin
        @test heisenberg_two_spinon_lower_edge(model, π) ≈ 0.0 atol = 1e-13
        @test heisenberg_two_spinon_upper_edge(model, π) ≈ π atol = 1e-13
    end

    @testset "ε_U(q) ≥ ε_L(q) on (0, π)" begin
        # Strictly above except at gapless points 0, π where both vanish.
        for q in range(0.05, π - 0.05; length=20)
            @test heisenberg_two_spinon_upper_edge(model, q) >
                heisenberg_two_spinon_lower_edge(model, q) - 1e-14
        end
    end

    # ─────────────────────────────────────────────────────────────────────
    # Müller-ansatz S^{zz}(q, ω): support, finiteness, edge singularity
    # ─────────────────────────────────────────────────────────────────────
    @testset "Müller S^{zz}: zero outside [ε_L, ε_U]" begin
        q = 2π / 3
        εL = heisenberg_two_spinon_lower_edge(model, q)
        εU = heisenberg_two_spinon_upper_edge(model, q)
        # Below lower edge
        for ω in (0.0, 0.5 * εL, εL)
            @test QAtlas.fetch(Heisenberg1D(), ZZStructureFactor(), Infinite(); q=q, ω=ω) ==
                0.0
        end
        # Above upper edge
        for ω in (εU, εU + 0.1, 10 * εU)
            @test QAtlas.fetch(Heisenberg1D(), ZZStructureFactor(), Infinite(); q=q, ω=ω) ==
                0.0
        end
    end

    @testset "Müller S^{zz}: finite & positive inside" begin
        q = π / 2
        εL = heisenberg_two_spinon_lower_edge(model, q)
        εU = heisenberg_two_spinon_upper_edge(model, q)
        # mid-band sample
        ωmid = (εL + εU) / 2
        S = QAtlas.fetch(Heisenberg1D(), ZZStructureFactor(), Infinite(); q=q, ω=ωmid)
        @test isfinite(S)
        @test S > 0
        # closed-form check: 1 / (2 √(ω² − ε_L²))
        @test S ≈ 1 / (2 * sqrt(ωmid^2 - εL^2)) atol = 1e-14
    end

    @testset "Müller S^{zz}: integrable singularity ω → ε_L⁺" begin
        q = π / 2
        εL = heisenberg_two_spinon_lower_edge(model, q)
        # Sequence of points approaching ε_L from above; S should grow like
        # 1 / √(ω − ε_L) and exceed any finite constant for ω close enough.
        S_values = Float64[]
        for δ in (1e-2, 1e-4, 1e-6)
            push!(
                S_values,
                QAtlas.fetch(
                    Heisenberg1D(), ZZStructureFactor(), Infinite(); q=q, ω=εL + δ
                ),
            )
        end
        # Strictly monotone divergence
        @test S_values[1] < S_values[2] < S_values[3]
        # Exceeds 1e2 well before machine precision
        @test S_values[3] > 1e2
    end

    @testset "Müller S^{zz}: explicit Müller method symbol works" begin
        q = π / 2
        ω = 1.0
        S_default = QAtlas.fetch(Heisenberg1D(), ZZStructureFactor(), Infinite(); q=q, ω=ω)
        S_muller = QAtlas.fetch(
            Heisenberg1D(), ZZStructureFactor(), Infinite(); q=q, ω=ω, method=:muller
        )
        @test S_default == S_muller
    end

    @testset "Exact 2-spinon S^{zz}: support, finiteness, and limits" begin
        q = π / 2
        εL = heisenberg_two_spinon_lower_edge(model, q)
        εU = heisenberg_two_spinon_upper_edge(model, q)

        # Zero outside support
        for ω in (0.0, 0.5 * εL, εL)
            @test QAtlas.fetch(
                Heisenberg1D(),
                ZZStructureFactor(),
                Infinite();
                q=q,
                ω=ω,
                method=:exact_2spinon,
            ) == 0.0
        end
        for ω in (εU, εU + 0.1, 10.0 * εU)
            @test QAtlas.fetch(
                Heisenberg1D(),
                ZZStructureFactor(),
                Infinite();
                q=q,
                ω=ω,
                method=:exact_2spinon,
            ) == 0.0
        end

        # Finite inside
        ωmid = (εL + εU) / 2
        S_exact = QAtlas.fetch(
            Heisenberg1D(),
            ZZStructureFactor(),
            Infinite();
            q=q,
            ω=ωmid,
            method=:exact_2spinon,
        )
        @test isfinite(S_exact)
        @test S_exact > 0

        # Divergence at ε_L
        S_near = QAtlas.fetch(
            Heisenberg1D(),
            ZZStructureFactor(),
            Infinite();
            q=q,
            ω=εL + 1e-5,
            method=:exact_2spinon,
        )
        @test S_near > S_exact
        @test S_near > 1.0
    end

    @testset "Müller S^{zz}: Phase-2 :caux_hagemans raises informative error" begin
        @test_throws ErrorException QAtlas.fetch(
            Heisenberg1D(),
            ZZStructureFactor(),
            Infinite();
            q=π / 2,
            ω=1.0,
            method=:caux_hagemans,
        )
        @test_throws ErrorException QAtlas.fetch(
            Heisenberg1D(),
            ZZStructureFactor(),
            Infinite();
            q=π / 2,
            ω=1.0,
            method=:bogus_method,
        )
    end

    # ─────────────────────────────────────────────────────────────────────
    # J-scaling
    # ─────────────────────────────────────────────────────────────────────
    @testset "linear J-scaling of dispersion + edges" begin
        for J in (0.5, 2.0, 3.7)
            for q in (0.3, π / 2, 2.1)
                @test heisenberg_spinon_dispersion(model, q; J=J) ≈
                    J * heisenberg_spinon_dispersion(model, q) atol = 1e-14
                @test heisenberg_two_spinon_lower_edge(model, q; J=J) ≈
                    J * heisenberg_two_spinon_lower_edge(model, q) atol = 1e-14
                @test heisenberg_two_spinon_upper_edge(model, q; J=J) ≈
                    J * heisenberg_two_spinon_upper_edge(model, q) atol = 1e-14
            end
        end
    end

    @testset "S^{zz} dimensional scaling under (J, ω) → (λ J, λ ω)" begin
        # ε_L, ε_U scale linearly with J, and inside the continuum
        # S^{zz}(q, ω; J) = 1 / (2 √(ω² − ε_L²)) carries one inverse power
        # of energy.  Therefore S^{zz}(q, λ ω; λ J) = (1/λ) S^{zz}(q, ω; J).
        q = π / 2
        ω = 1.0
        S1 = QAtlas.fetch(Heisenberg1D(), ZZStructureFactor(), Infinite(); q=q, ω=ω, J=1.0)
        for λ in (0.5, 2.0, 3.0)
            S_scaled = QAtlas.fetch(
                Heisenberg1D(), ZZStructureFactor(), Infinite(); q=q, ω=λ * ω, J=λ * 1.0
            )
            @test S_scaled ≈ S1 / λ atol = 1e-14
        end
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "Heisenberg1D spinon ZZStructureFactor — verification cards" begin
    # Independent re-derivation of the Müller ansatz from the
    # des Cloizeaux-Pearson two-spinon dispersion bounds (not read from src):
    #   ε_L(q) = (π J / 2) |sin q|,  ε_U(q) = π J |sin(q/2)|
    #   S^zz(q,ω) = Θ[ω-ε_L] Θ[ε_U-ω] / (2 √(ω² - ε_L²))
    let J = 1.0, q = π / 2, ω = 2.0
        εL = (π * J / 2) * abs(sin(q))
        εU = π * J * abs(sin(q / 2))
        @assert εL < ω < εU "pick (q,ω) inside the two-spinon continuum"
        S_indep = 1 / (2 * sqrt(ω^2 - εL^2))
        verify(
            Heisenberg1D(),
            ZZStructureFactor(),
            Infinite();
            route=:second_closed_form,
            fetch_kw=(; q=q, ω=ω, J=J, method=:muller),
            independent=S_indep,
            agree_within=1e-12,
            refs=["Müller-Thomas-Beck-Bonner 1981; des Cloizeaux-Pearson 1962 dispersion"],
        )
    end

    # Support boundary: S = 0 just above the upper continuum edge ε_U
    let J = 1.0, q = π / 2
        εU = π * J * abs(sin(q / 2))
        verify(
            Heisenberg1D(),
            ZZStructureFactor(),
            Infinite();
            route=:second_closed_form,
            fetch_kw=(; q=q, ω=εU + 0.5, J=J, method=:muller),
            independent=0.0,
            agree_within=1e-14,
            refs=["Two-spinon continuum has compact support: S=0 for ω > ε_U(q)"],
        )
    end
end
