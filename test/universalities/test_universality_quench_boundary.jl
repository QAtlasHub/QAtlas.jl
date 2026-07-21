# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: Post-quench entanglement dynamics and Affleck-Ludwig
# boundary / residual terms at the Universality{C} level.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Universality :: Quench Dynamics & Boundary Entropies" begin
    # ── 1. Post-Quench Entanglement (QuenchEntanglementEntropy) ─────────────
    @testset "Infinite BC global quench" begin
        # c = 1/2 (Ising)
        model = Universality(:Ising)
        ℓ = 10.0
        v = 1.5
        beta_eff = 2.0
        c = 0.5

        # t = 1.0 (t < ℓ/(2v) = 10/3 = 3.333...) -> linear growth
        # S = 2 * v * t * pi * c / (6 * beta_eff) = 2 * 1.5 * 1.0 * pi * 0.5 / 12 = 0.125 * pi
        t1 = 1.0
        S1 = QAtlas.fetch(
            model,
            QuenchEntanglementEntropy(),
            Infinite();
            ℓ=ℓ,
            t=t1,
            v=v,
            beta_eff=beta_eff,
        )
        @test S1 ≈ (2 * v * t1) * π * c / (6 * beta_eff) atol=1e-12

        # t = 5.0 (t > ℓ/(2v) = 3.333...) -> saturation
        # S = ℓ * pi * c / (6 * beta_eff) = 10 * pi * 0.5 / 12 = 5/12 * pi
        t2 = 5.0
        S2 = QAtlas.fetch(
            model,
            QuenchEntanglementEntropy(),
            Infinite();
            ℓ=ℓ,
            t=t2,
            v=v,
            beta_eff=beta_eff,
        )
        @test S2 ≈ ℓ * π * c / (6 * beta_eff) atol=1e-12
    end

    @testset "OBC global quench with boundary state" begin
        # c = 1 (XY, K=1)
        model = Universality(:XY)
        ℓ = 5.0
        v = 2.0
        beta_eff = 1.0
        c = 1.0

        # t = 1.0 (t < ℓ/v = 2.5) -> linear growth.
        # OBC block at the open boundary has ONE entanglement cut, so a single
        # quasiparticle front (speed v) fills the region by t = ℓ/v. The thermal
        # entropy density πc/(6β_eff) is the same as the bulk (it is the entropy
        # of ℓ sites, independent of the boundary); only the filling time differs.
        t1 = 1.0
        S1 = QAtlas.fetch(
            model, QuenchEntanglementEntropy(), OBC(); ℓ=ℓ, t=t1, v=v, beta_eff=beta_eff
        )
        @test S1 ≈ (v * t1) * π * c / (6 * beta_eff) atol=1e-12

        # t = 4.0 (t > ℓ/v = 2.5) -> saturation at the thermal entropy πcℓ/(6β_eff)
        t2 = 4.0
        S2 = QAtlas.fetch(
            model, QuenchEntanglementEntropy(), OBC(); ℓ=ℓ, t=t2, v=v, beta_eff=beta_eff
        )
        @test S2 ≈ ℓ * π * c / (6 * beta_eff) atol=1e-12

        # With boundary state / log g correction
        # For XY with boundary_state = :neumann (free), K = 1.0:
        # log g = 0.25 * log(1/2) = -0.25 * log(2)
        log_g = 0.25 * log(0.5)
        S1_bdy = QAtlas.fetch(
            model,
            QuenchEntanglementEntropy(),
            OBC();
            ℓ=ℓ,
            t=t1,
            v=v,
            beta_eff=beta_eff,
            boundary_state=:neumann,
            K=1.0,
        )
        @test S1_bdy ≈ S1 + log_g atol=1e-12

        S1_direct = QAtlas.fetch(
            model,
            QuenchEntanglementEntropy(),
            OBC();
            ℓ=ℓ,
            t=t1,
            v=v,
            beta_eff=beta_eff,
            log_g=0.5,
        )
        @test S1_direct ≈ S1 + 0.5 atol=1e-12
    end

    @testset "PBC global quench (revivals)" begin
        # c = 1 (Heisenberg, SU(2)_1 WZW)
        model = Universality(:Heisenberg)
        L = 20.0
        ℓ = 6.0
        v = 1.0
        beta_eff = 3.0
        c = 1.0

        # ℓ_c = min(6.0, 14.0) = 6.0
        # t_revival = L/v = 20.0
        # Time points:
        # t = 1.0: 0 <= t < ℓ_c/2v = 3.0 -> S = 2 * v * t * coeff
        # t = 5.0: 3.0 <= t < (L-ℓ_c)/2v = 7.0 -> S = ℓ_c * coeff
        # t = 9.0: 7.0 <= t < L/2v = 10.0 -> S = (L - 2 * v * t) * coeff
        coeff = π * c / (6 * beta_eff)

        S_t1 = QAtlas.fetch(
            model,
            QuenchEntanglementEntropy(),
            PBC();
            ℓ=ℓ,
            L=L,
            t=1.0,
            v=v,
            beta_eff=beta_eff,
        )
        @test S_t1 ≈ 2.0 * coeff atol=1e-12

        S_t5 = QAtlas.fetch(
            model,
            QuenchEntanglementEntropy(),
            PBC();
            ℓ=ℓ,
            L=L,
            t=5.0,
            v=v,
            beta_eff=beta_eff,
        )
        @test S_t5 ≈ 6.0 * coeff atol=1e-12

        S_t9 = QAtlas.fetch(
            model,
            QuenchEntanglementEntropy(),
            PBC();
            ℓ=ℓ,
            L=L,
            t=9.0,
            v=v,
            beta_eff=beta_eff,
        )
        @test S_t9 ≈ (20.0 - 18.0) * coeff atol=1e-12

        # Check periodicity by adding L/v = 20.0 to times
        S_t21 = QAtlas.fetch(
            model,
            QuenchEntanglementEntropy(),
            PBC();
            ℓ=ℓ,
            L=L,
            t=21.0,
            v=v,
            beta_eff=beta_eff,
        )
        @test S_t21 ≈ S_t1 atol=1e-12

        S_t49 = QAtlas.fetch(
            model,
            QuenchEntanglementEntropy(),
            PBC();
            ℓ=ℓ,
            L=L,
            t=49.0,
            v=v,
            beta_eff=beta_eff,
        )
        @test S_t49 ≈ S_t9 atol=1e-12
    end

    # ── 2. Boundary Entropy (BoundaryEntropy) for XY and Heisenberg ──────────
    @testset "BoundaryEntropy values" begin
        # Ising c = 1/2
        @test QAtlas.fetch(
            Universality(:Ising), BoundaryEntropy(), Infinite(); boundary_state=:free
        ) ≈ 0.0 atol=1e-12
        @test QAtlas.fetch(
            Universality(:Ising), BoundaryEntropy(), Infinite(); boundary_state=:fixed_up
        ) ≈ -log(2) / 2 atol=1e-12

        # XY c = 1
        # log g_N = 0.25 * log(K/2), log g_D = -0.25 * log(2K)
        @test QAtlas.fetch(
            Universality(:XY), BoundaryEntropy(), Infinite(); boundary_state=:neumann, K=1.0
        ) ≈ 0.25 * log(0.5) atol=1e-12
        @test QAtlas.fetch(
            Universality(:XY),
            BoundaryEntropy(),
            Infinite();
            boundary_state=:dirichlet,
            K=1.0,
        ) ≈ -0.25 * log(2.0) atol=1e-12
        @test QAtlas.fetch(
            Universality(:XY), BoundaryEntropy(), Infinite(); boundary_state=:neumann, K=2.0
        ) ≈ 0.0 atol=1e-12
        @test QAtlas.fetch(
            Universality(:XY),
            BoundaryEntropy(),
            Infinite();
            boundary_state=:dirichlet,
            K=2.0,
        ) ≈ -0.25 * log(4.0) atol=1e-12

        # Heisenberg c = 1 (SU(2)_1 WZW). The two Cardy states j=0, j=1/2 are
        # g-degenerate: g_a = S_{0a}/sqrt(S_{00}) with S_{00}=S_{0,1/2}=1/sqrt(2)
        # gives g = 2^{-1/4} for both, log g = -(1/4) log 2. (Matches the :XY
        # compact boson at the self-dual radius K=1.)
        @test QAtlas.fetch(
            Universality(:Heisenberg), BoundaryEntropy(), Infinite(); boundary_state=:free
        ) ≈ -log(2) / 4 atol=1e-12
        @test QAtlas.fetch(
            Universality(:Heisenberg), BoundaryEntropy(), Infinite(); boundary_state=:fixed
        ) ≈ -log(2) / 4 atol=1e-12
        # cross-check: :XY at K=1 (self-dual) gives the same g = 2^{-1/4}
        @test QAtlas.fetch(
            Universality(:XY), BoundaryEntropy(), Infinite(); boundary_state=:neumann, K=1.0
        ) ≈ -log(2) / 4 atol=1e-12
        @test QAtlas.fetch(
            Universality(:XY),
            BoundaryEntropy(),
            Infinite();
            boundary_state=:dirichlet,
            K=1.0,
        ) ≈ -log(2) / 4 atol=1e-12
    end

    # ── 3. OBC Entanglement with Boundary Entropy ────────────────────────────
    @testset "OBC equilibrium entanglement boundary corrections" begin
        model = Universality(:Ising)
        ℓ = 4.0
        L = 8.0

        # Base S_vn and S_renyi
        S_vn_base = QAtlas.fetch(model, VonNeumannEntropy(), OBC(); ℓ=ℓ, L=L)
        S_ren_base = QAtlas.fetch(model, RenyiEntropy(2.0), OBC(); ℓ=ℓ, L=L)

        # Added boundary state :fixed_up (log g = -log(2)/2)
        log_g = -log(2) / 2
        S_vn_fixed = QAtlas.fetch(
            model, VonNeumannEntropy(), OBC(); ℓ=ℓ, L=L, boundary_state=:fixed_up
        )
        @test S_vn_fixed ≈ S_vn_base + log_g atol=1e-12

        S_ren_fixed = QAtlas.fetch(
            model, RenyiEntropy(2.0), OBC(); ℓ=ℓ, L=L, boundary_state=:fixed_up
        )
        @test S_ren_fixed ≈ S_ren_base + log_g atol=1e-12

        # Added explicit log_g = 0.5
        S_vn_direct = QAtlas.fetch(model, VonNeumannEntropy(), OBC(); ℓ=ℓ, L=L, log_g=0.5)
        @test S_vn_direct ≈ S_vn_base + 0.5 atol=1e-12
    end
end
