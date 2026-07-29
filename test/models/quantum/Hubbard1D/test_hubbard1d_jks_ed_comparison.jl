# test/models/quantum/Hubbard1D/test_hubbard1d_jks_ed_comparison.jl
# ED-based mid-T verification for the JKS NLIE.

using Test
using LinearAlgebra
using SparseArrays
using QAtlas
using QAtlas: Hubbard1D, FreeEnergy, Infinite
using QAtlas.Hubbard1DJKSNLIE: atomic_free_energy

include(joinpath(@__DIR__, "..", "..", "..", "util", "hubbard_ed.jl"))

@testset "Hubbard1D JKS NLIE — ED N=4 mid-T comparison (#523)" begin
    @testset "ED essentially equals atomic at β ≤ 0.2 (large U regime)" begin
        # For U=4, β ≤ 0.2: kinetic correction is t²/U ~ 0.25 with β prefactor,
        # so f_ED should equal f_atom to within ~1%%.
        for β in (0.001, 0.01, 0.05, 0.1, 0.2)
            f_ed = _ed_hubbard_free_energy(4, 1.0, 4.0, 2.0, β; pbc=false)
            f_atom = atomic_free_energy(β, 4.0, 2.0)
            @test isapprox(f_ed, f_atom; rtol=0.05)
        end
    end

    @testset "JKS at high T matches ED to <0.5%%" begin
        m = Hubbard1D(; t=1.0, U=4.0, μ=2.0)
        for β in (1e-4, 1e-3)
            f_ed = _ed_hubbard_free_energy(4, 1.0, 4.0, 2.0, β; pbc=false)
            f_jks = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
            @test isapprox(f_jks, f_ed; rtol=0.005)
        end
    end

    @testset "JKS mid-T deviates from ED >5%% — current implementation bug regression guard" begin
        # Stage G.2: ED proved mid-T ratio drop is a JKS formula bug, not physics.
        # This test PINS the buggy behaviour so any future fix is detected.
        m = Hubbard1D(; t=1.0, U=4.0, μ=2.0)
        for β in (0.1, 0.2)
            f_ed = _ed_hubbard_free_energy(4, 1.0, 4.0, 2.0, β; pbc=false)
            f_jks = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
            relative_error = abs(f_jks - f_ed) / abs(f_ed)
            @test_broken relative_error < 0.05  # currently fails ~25-50%%
        end
    end
end
