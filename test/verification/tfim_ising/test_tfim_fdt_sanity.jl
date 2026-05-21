# ─────────────────────────────────────────────────────────────────────────────
# TFIM FDT sanity (chi-prime antisymmetry + static branch)
#
# Split out of test/verification/tfim_ising/test_tfim_fdt.jl (6.2 min on s03)
# so each beta slice runs on its own shard. Same FDT physics:
#     S_zz(q, omega; beta) = (2 / (1 - exp(-beta * omega))) * chi''_zz(q, omega; beta)
# Reference: Kubo 1957; Mahan, *Many-Particle Physics*, ch. 3.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM Infinite() SusceptibilityZZ dynamic-branch sanity" begin
    N_proxy = 24
    t_max = 10.0
    dt = 0.1
    model = TFIM(; J=1.0, h=0.5)
    q = π / 2

    # χ'' antisymmetry in ω (FDT-independent sanity).
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

    # Static branch: per-site χ_zz(β) at q = 0 (kwarg q absent).
    χ_static = QAtlas.fetch(model, SusceptibilityZZ(), Infinite(); beta=1.0)
    @test isfinite(χ_static)
    @test χ_static > 0
end
