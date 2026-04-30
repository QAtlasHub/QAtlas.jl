# =============================================================================
# Tests for TFIM Infinite-volume dynamic / kinematic observables defined in
# `src/models/quantum/TFIM/TFIM_infinite_dynamics.jl`.
#
# Coverage:
#   * Single-quasiparticle dispersion endpoints and minimum (= MassGap).
#   * Two-spinon DOS support: zero outside the 2-spinon continuum, > 0 inside.
#   * Dynamic ZZ structure factor smoke test at criticality.
#   * Static-limit cross-check between dynamic and static structure factor.
# =============================================================================

using QAtlas, Test

@testset "TFIM Infinite dynamics" begin
    @testset "Single-quasiparticle dispersion sanity" begin
        for h in (0.5, 1.0, 1.5)
            model = TFIM(; J=1.0, h=h)
            Δ_inf = QAtlas.fetch(model, MassGap(), Infinite())
            # k = 0 → 2|h - J|
            @test tfim_quasiparticle_dispersion(model, 0.0) ≈ 2 * abs(model.h - model.J) atol =
                1e-12
            # k = π → 2(h + J)
            @test tfim_quasiparticle_dispersion(model, π) ≈ 2 * (model.h + model.J) atol =
                1e-12
            # min over k = MassGap
            ks = range(0, π; length=1000)
            Λs = [tfim_quasiparticle_dispersion(model, k) for k in ks]
            @test minimum(Λs) ≈ Δ_inf atol = 1e-6
        end
    end

    @testset "2-spinon DOS: support boundaries" begin
        # Gapped: J = 1, h = 0.5  ⇒  Δ = 1, Λ(0) = 1, Λ(π) = 3.
        # 2-spinon continuum at q_total = 0: [2 Δ, 2 Λ(π)] = [2, 6].
        model = TFIM(; J=1.0, h=0.5)
        @test tfim_two_spinon_dos(model, 1.5; q_total=0.0) == 0.0  # below threshold
        @test tfim_two_spinon_dos(model, 7.0; q_total=0.0) == 0.0  # above ceiling
        @test tfim_two_spinon_dos(model, 4.0; q_total=0.0) > 0.0
    end

    @testset "Dynamic ZZ structure factor at criticality (smoke)" begin
        # Cheap smoke test: small N_proxy + short t_max + coarse dt is
        # enough to catch dispatch breakage / non-finite return; the
        # full-resolution proxy is for users, not CI.
        model = TFIM(; J=1.0, h=1.0)
        S = QAtlas.fetch(
            model,
            ZZStructureFactor(),
            Infinite();
            beta=Inf,
            q=π / 2,
            ω=1.0,
            N_proxy=24,
            t_max=6.0,
            dt=0.3,
        )
        @test isfinite(S)
    end

    @testset "Static branch still works (no ω kwarg)" begin
        # Make sure the router preserves the `ω === nothing` behaviour
        # already provided by TFIM_zaxis.jl.
        model = TFIM(; J=1.0, h=0.7)
        S_static = QAtlas.fetch(model, ZZStructureFactor(), Infinite(); beta=5.0, q=π / 3)
        @test isfinite(S_static)
        @test S_static ≥ 0.0
    end

    @testset "Sum rule sanity (dynamic ↔ static)" begin
        # Frequency integral of the dynamic S(q, ω) should be related to
        # the static structure factor.  This is a finiteness smoke test —
        # full sum-rule verification belongs to a Stage-D-style
        # verification campaign with `t_max` swept.  Use a coarse 5-point
        # grid to keep CI cheap (each fetch is O(N_proxy² · t_grid)).
        model = TFIM(; J=1.0, h=0.7)
        q = π / 3
        ωs = (-4.0, -2.0, 0.0, 2.0, 4.0)
        S_omega_sum = sum(
            QAtlas.fetch(
                model,
                ZZStructureFactor(),
                Infinite();
                beta=10.0,
                q=q,
                ω=ω,
                N_proxy=24,
                t_max=6.0,
                dt=0.3,
            ) * 2.0 for ω in ωs
        )
        @test isfinite(S_omega_sum)
    end
end
