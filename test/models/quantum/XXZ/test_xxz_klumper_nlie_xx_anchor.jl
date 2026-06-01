# test/models/quantum/XXZ/test_xxz_klumper_nlie_xx_anchor.jl
#
# Part of the Klümper-NLIE-backed FreeEnergy@Infinite test suite for
# `XXZ1D`, split from the original `test_xxz_klumper_nlie.jl` so each
# subset runs in its own CI shard.

using Test
using QAtlas
using QAtlas: XXZ1D, FreeEnergy, ThermalEntropy, SpecificHeat, Energy, Infinite

@testset "XXZ Klümper NLIE — XX regression anchor (Δ = 1e-6)" begin
    # Regression guard against dispatch-ordering bugs: at Δ → 0⁺ the
    # NLIE branch should agree with the closed-form XX result to ~1e-3.
    # Forces a fresh γ ≈ π/2 grid build (~2 min on CI workers).
    @testset "NLIE near XX (Δ = 1e-6) matches closed form" begin
        m_xx = XXZ1D(; J=1.0, Δ=0.0)
        m_eps = XXZ1D(; J=1.0, Δ=1e-6)
        β = 1.0
        f_xx = QAtlas.fetch(m_xx, FreeEnergy(), Infinite(); beta=β)
        f_eps = QAtlas.fetch(m_eps, FreeEnergy(), Infinite(); beta=β)
        @test isapprox(f_eps, f_xx; rtol=1e-3)
    end
end
