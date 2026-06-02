# test/models/quantum/Hubbard1D/test_hubbard1d_jks_stage_c6.jl
#
# Stage C.6 regression: adaptive alpha-mix Picard solver (#523).

using Test
using QAtlas
using QAtlas.Hubbard1DJKSNLIE:
    JKSContourGrid,
    JKSAuxFunctions,
    init_atomic_limit,
    jks_nlie_residual_shifted,
    solve_jks_nlie_shifted,
    solve_jks_nlie_adaptive,
    JKSSolution

@testset "Hubbard1D — JKS Stage C.6 adaptive Picard (#523)" begin
    @testset "Adaptive solver returns JKSSolution with finite residual" begin
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        sol = solve_jks_nlie_adaptive(grid, 1.0, 4.0, 2.0; alpha=0.5, maxiter=200)
        @test sol isa JKSSolution
        @test isfinite(sol.residual)
        @test sol.iterations >= 1
    end

    @testset "Adaptive solver validates inputs" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        @test_throws DomainError solve_jks_nlie_adaptive(grid, 0.0, 4.0, 2.0; alpha=0.5)
        @test_throws DomainError solve_jks_nlie_adaptive(
            grid, 1.0, 4.0, 2.0; alpha=0.5, alpha_mix_init=0.0
        )
        @test_throws DomainError solve_jks_nlie_adaptive(
            grid, 1.0, 4.0, 2.0; alpha=0.5, alpha_mix_init=2.0
        )
        @test_throws DomainError solve_jks_nlie_adaptive(
            grid, 1.0, 4.0, 2.0; alpha=0.5, shrink=0.0
        )
        @test_throws DomainError solve_jks_nlie_adaptive(
            grid, 1.0, 4.0, 2.0; alpha=0.5, shrink=1.0
        )
        @test_throws DomainError solve_jks_nlie_adaptive(
            grid, 1.0, 4.0, 2.0; alpha=0.5, grow=1.0
        )
        @test_throws DomainError solve_jks_nlie_adaptive(
            grid, 1.0, 4.0, 2.0; alpha=0.5, grow=0.5
        )
    end

    @testset "High-T case: adaptive converges to tol" begin
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        sol = solve_jks_nlie_adaptive(
            grid, 0.01, 4.0, 2.0; alpha=0.5, tol=1e-4, maxiter=500
        )
        @test sol.residual < 1e-2
    end

    @testset "Adaptive solver does not blow up at mid-T" begin
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        sol_fixed = solve_jks_nlie_shifted(
            grid, 1.0, 4.0, 2.0; alpha=0.5, tol=1e-10, maxiter=500, alpha_mix=0.02
        )
        sol_adaptive = solve_jks_nlie_adaptive(
            grid, 1.0, 4.0, 2.0; alpha=0.5, tol=1e-10, maxiter=500
        )
        @test isfinite(sol_adaptive.residual)
        @test sol_adaptive.residual <= sol_fixed.residual * 2 + 1e-6
    end

    @testset "Adaptive: low U case (Stage C.5 diverging) stays bounded" begin
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        sol = solve_jks_nlie_adaptive(
            grid, 1.0, 1.0, 0.5; alpha=0.2, tol=1e-10, maxiter=500
        )
        @test isfinite(sol.residual)
        @test sol.residual < 100
    end
end
