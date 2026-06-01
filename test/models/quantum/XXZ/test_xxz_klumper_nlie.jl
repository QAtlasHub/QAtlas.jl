# test/models/quantum/XXZ/test_xxz_klumper_nlie.jl
#
# Tier-1 smoke for the Klümper NLIE dispatch on `XXZ1D` at `Infinite()`.
# Cross-checks:
#  (a) Δ = 0 (XX) still routes to the free-fermion closed form.
#  (b) Critical Δ ∈ (-1, 1) high-T limit: f - e_0 → -T ln 2 as β → 0.
#  (c) Δ = 0.99 / Δ = 1.5 endpoints return NaN with a warning.
# The NLIE grid build is ~90 s on a laptop (one-time cost, cached
# across β); we keep this test minimal so it does not gate CI.

using Test
using QAtlas
using QAtlas: XXZ1D, FreeEnergy, Energy, Infinite

@testset "XXZ Klümper NLIE — Infinite FreeEnergy" begin
    @testset "XX limit (Δ=0) routes to closed form" begin
        m = XXZ1D(; J=1.0, Δ=0.0)
        f = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=1.0)
        @test isfinite(f)
        @test f < 0
    end

    @testset "high-T limit at Δ=0.5: f - e_0 → -T ln 2" begin
        m = XXZ1D(; J=1.0, Δ=0.5)
        β = 0.001
        f = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
        e0 = QAtlas.fetch(m, Energy{:per_site}(), Infinite())
        excess = f - e0
        target = -log(2) / β
        @test isfinite(f)
        @test isapprox(excess, target; rtol=1e-3)
    end

    @testset "gapped Δ=1.5 returns NaN with warning" begin
        m = XXZ1D(; J=1.0, Δ=1.5)
        f = (@test_logs (:warn,) match_mode=:any QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=1.0))
        @test isnan(f)
    end

    @testset "near-endpoint Δ=0.999 also returns NaN" begin
        m = XXZ1D(; J=1.0, Δ=0.999)
        f = (@test_logs (:warn,) match_mode=:any QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=1.0))
        @test isnan(f)
    end
end
