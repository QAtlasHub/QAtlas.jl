# test/models/quantum/Hubbard1D/test_hubbard1d_jks_thermo_sweep.jl
#
# Detailed multi-temperature benchmark for the JKS NLIE FreeEnergy@Infinite
# dispatch (#523). Mirrors the structure of
# test/models/quantum/XXZ/test_xxz_klumper_nlie_thermo.jl so the two
# finite-T NLIE ports can be cross-compared at a glance.
#
# Part 1/2 (the temperature sweeps): split off the cross-validation half
# (test_hubbard1d_jks_thermo_sweep2.jl) so CI can run the two halves on
# separate shards — the single file was ~9 min from the repeated NLIE solves.

using Test
using QAtlas
using QAtlas: Hubbard1D, FreeEnergy, Infinite
using QAtlas.Hubbard1DJKSNLIE: atomic_free_energy

@testset "Hubbard1D JKS NLIE — temperature sweep (#523, part 1/2)" begin
    @testset "Atomic limit at 7 temperature points (rtol scaled with β)" begin
        # Kinetic corrections O(β t / U) grow with β; allow looser rtol
        # away from β → 0 to track the physical kinetic shift.
        m = Hubbard1D(; t=1.0, U=4.0, μ=2.0)
        cases = (
            (1e-5, 0.02),
            (1e-4, 0.02),
            (1e-3, 0.02),
            (1e-2, 0.05),
            (5e-2, 0.15),
            (1e-1, 0.30),
            (2e-1, 0.50),
        )
        for (β, rtol) in cases
            f = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
            f_atom = atomic_free_energy(β, 4.0, 2.0)
            @test isfinite(f)
            @test isapprox(f, f_atom; rtol=rtol)
        end
    end

    @testset "f(β) monotone non-decreasing in β" begin
        m = Hubbard1D(; t=1.0, U=4.0, μ=2.0)
        betas = (1e-3, 5e-3, 1e-2, 5e-2, 1e-1)
        fs = [QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β) for β in betas]
        @test all(isfinite, fs)
        @test all(diff(fs) .>= -1e-3 * maximum(abs.(fs)))
    end
end
