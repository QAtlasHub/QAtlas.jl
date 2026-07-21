# ─────────────────────────────────────────────────────────────────────────────
# TFIM FDT residual at beta = 1.0
#
# Split out of test/verification/tfim_ising/test_tfim_fdt.jl (6.2 min on s03)
# so each beta slice runs on its own shard. Same FDT physics:
#     S_zz(q, omega; beta) = (2 / (1 - exp(-beta * omega))) * chi''_zz(q, omega; beta)
# Reference: Kubo 1957; Mahan, *Many-Particle Physics*, ch. 3.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "FDT residual at TFIM Infinite(), β = 1.0 (verify cards per ω)" begin
    # CI-budget proxy parameters identical to the pre-split file.
    N_proxy = 24
    t_max = 10.0
    dt = 0.1

    model = TFIM(; J=1.0, h=0.5)        # gapped (h < J), Δ = 1
    q = π / 2
    ω_grid = collect(range(0.6, 3.0; length=5))
    β = 1.0

    # Migrated from aggregate `@test rel_residual < 1e-1` to per-ω verify
    # :sum_rule cards (PR #449 phase 2 zero-legacy). Each ω becomes one
    # card: subject = S_zz(q,ω;β), independent = 2 χ''(q,ω;β) / (1-e^{-βω}).
    for ω in ω_grid
        χpp = QAtlas.fetch(
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
        verify(
            model,
            DynamicalSpinStructureFactor(:z, :z),
            Infinite();
            route=:sum_rule,
            fetch_kw=(; beta=β, q=q, ω=ω, N_proxy=N_proxy, t_max=t_max, dt=dt),
            independent=2 * χpp / (1 - exp(-β * ω)),
            agree_within=5e-2,
            at=["β=$(β)", "q=$(q)", "ω=$(ω)"],
            refs=[
                "Kubo FDT: S_zz(q,ω;β) = 2 χ''_zz(q,ω;β) / (1 − e^{-βω}); proxy params N_proxy=24, t_max=10, dt=0.1 (Mahan ch.3)",
            ],
        )
    end
end
