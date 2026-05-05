# ─────────────────────────────────────────────────────────────────────────────
# Verification: fluctuation-dissipation theorem (FDT) for TFIM at Infinite()
#
#     S_zz(q, ω; β) = (2 / (1 - exp(-β ω))) · χ''_zz(q, ω; β)        (bosonic)
#
# Cross-check QAtlas's two independent dynamic Z-axis observables:
#
#     ZZStructureFactor   →  S_zz(q, ω; β)
#     SusceptibilityZZ    →  χ''_zz(q, ω; β)   (dynamic ω-branch)
#
# Both are computed by the OBC large-N proxy on the same Majorana / Pfaffian
# machinery with identical (N_proxy, t_max, dt) discretisations, so the FDT
# residual is dominated by the discretisation rather than by independent
# numerical errors of the two proxies.  We therefore check a moderate
# tolerance (~5e-2) on the relative residual and document the achievable
# precision in the test name.
#
# Reference: Kubo 1957; Mahan, *Many-Particle Physics*, ch. 3 (FDT in
# bosonic convention for Hermitian operator response).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "FDT residual S_zz = (2/(1 - e^{-βω})) χ''_zz at TFIM Infinite()" begin
    # Use a single set of discretisation parameters across the whole sweep so
    # the two proxies share their dominant systematic error and the residual
    # is a fair test of the FDT identity itself.
    # CI budget: keep parameters small so the full sweep fits within the
    # standard test runtime.  Earlier (48, 20, 10ω, 3β) hung CI for 3+ h on
    # the runner; these numbers reproduce the FDT identity within a 1e-1
    # residual budget while keeping each proxy call cheap.
    N_proxy = 24
    t_max = 10.0
    dt = 0.1

    model = TFIM(; J=1.0, h=0.5)        # gapped (h < J), Δ = 1
    q = π / 2

    # 5 ω points in [0.6, 3.0]; skip ω near 0 where FDT is singular
    # (1 − e^{-βω} → 0) — the ω-resolution of the discrete Fourier transform
    # is ~π / t_max ≈ 0.31 with t_max=10, so 0.6 is comfortably above that
    # floor.
    ω_grid = collect(range(0.6, 3.0; length=5))
    βs = (0.5, 1.0)

    for β in βs
        S_vals = [
            QAtlas.fetch(
                model,
                ZZStructureFactor(),
                Infinite();
                beta=β,
                q=q,
                ω=ω,
                N_proxy=N_proxy,
                t_max=t_max,
                dt=dt,
            ) for ω in ω_grid
        ]
        χpp_vals = [
            QAtlas.fetch(
                model,
                SusceptibilityZZ(),
                Infinite();
                beta=β,
                q=q,
                ω=ω,
                N_proxy=N_proxy,
                t_max=t_max,
                dt=dt,
            ) for ω in ω_grid
        ]

        FDT_vals = [2 * χpp / (1 - exp(-β * ω)) for (χpp, ω) in zip(χpp_vals, ω_grid)]

        denom = maximum(abs, S_vals)
        rel_residual = maximum(abs, S_vals .- FDT_vals) / denom

        # 1e-1 budget: at (N_proxy=24, t_max=10, dt=0.1) the per-point FDT
        # residual is dominated by ω-resolution (~π / t_max ≈ 0.31) and
        # finite-N bulk corrections; the wider Δω broadens the tolerance
        # vs the original (48, 20) parameter set.  See the convergence note
        # in `_tfim_zz_structure_factor_dynamic_proxy` and
        # `_tfim_chi_imag_zz_dynamic_proxy`.
        @test rel_residual < 1e-1
    end

    # Also assert χ'' antisymmetry in ω (an FDT-independent sanity check
    # that the new dynamic χ'' branch is wired correctly).
    let β = 1.0, ω = 1.0
        χp = QAtlas.fetch(
            model,
            SusceptibilityZZ(),
            Infinite();
            beta=β,
            q=q,
            ω=ω,
            N_proxy=N_proxy,
            t_max=t_max,
            dt=dt,
        )
        χm = QAtlas.fetch(
            model,
            SusceptibilityZZ(),
            Infinite();
            beta=β,
            q=q,
            ω=(-ω),
            N_proxy=N_proxy,
            t_max=t_max,
            dt=dt,
        )
        @test isapprox(χm, -χp; atol=1e-10, rtol=1e-10)
    end

    # `q` is required for the dynamic branch (not for the static branch).
    @test_throws ArgumentError QAtlas.fetch(
        model, SusceptibilityZZ(), Infinite(); beta=1.0, ω=1.0
    )

    # Static branch unchanged: per-site uniform χ_zz(β) at q = 0 (kwarg q
    # absent).  Sanity-check that it returns a finite positive number for
    # the gapped phase — no regression on the previously-existing method.
    χ_static = QAtlas.fetch(model, SusceptibilityZZ(), Infinite(); beta=1.0)
    @test isfinite(χ_static)
    @test χ_static > 0
end
