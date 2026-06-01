# test/models/quantum/XXZ/test_xxz_klumper_nlie_solver_unit.jl
#
# Part of the Klümper-NLIE-backed FreeEnergy@Infinite test suite for
# `XXZ1D`, split from the original `test_xxz_klumper_nlie.jl` so each
# subset runs in its own CI shard.

using Test
using QAtlas
using QAtlas: XXZ1D, FreeEnergy, ThermalEntropy, SpecificHeat, Energy, Infinite

@testset "XXZ Klümper NLIE — solver unit (small grid)" begin
    @testset "solve_klumper_nlie at maxiter=1 reports non-convergence" begin
        # Direct exercise of the !sol.converged branch in the wrapper.
        # The wrapper hard-codes maxiter=400 so we cannot trigger this
        # path via the public fetch. We call the internal solver directly
        # to confirm the contract: residual > tol and converged == false.
        γ = acos(0.5)
        grid = QAtlas.XXZKlumperNLIE.build_grid(γ; N=64, L_factor=10.0, ε_shift=0.5)
        β̃ = 1.0
        sol = QAtlas.XXZKlumperNLIE.solve_klumper_nlie(grid, β̃; α=0.4, maxiter=1, tol=1e-9)
        @test !sol.converged
        @test sol.residual > 1e-9
        @test sol.iterations == 1
    end
end
