# ─────────────────────────────────────────────────────────────────────────────
# TFIM FDT residual at beta = 1.0
#
# Split out of test/verification/tfim_ising/test_tfim_fdt.jl (6.2 min on s03)
# so each beta slice runs on its own shard. Same FDT physics:
#     S_zz(q, omega; beta) = (2 / (1 - exp(-beta * omega))) * chi''_zz(q, omega; beta)
# Reference: Kubo 1957; Mahan, *Many-Particle Physics*, ch. 3.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "FDT residual at TFIM Infinite(), β = 1.0" begin
    # CI-budget proxy parameters identical to the pre-split file.
    N_proxy = 24
    t_max = 10.0
    dt = 0.1

    model = TFIM(; J=1.0, h=0.5)        # gapped (h < J), Δ = 1
    q = π / 2
    ω_grid = collect(range(0.6, 3.0; length=5))
    β = 1.0

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

    # 1e-1 budget at (N_proxy=24, t_max=10, dt=0.1) — see pre-split header.
    @test rel_residual < 1e-1
end
