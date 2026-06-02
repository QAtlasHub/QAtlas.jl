# test/models/quantum/Hubbard1D/test_hubbard1d_jks_stage_c8.jl
#
# Stage C.8 regression: beta-continuation Newton solver (#523).

using Test
using QAtlas
using QAtlas.Hubbard1DJKSNLIE:
    JKSContourGrid,
    JKSAuxFunctions,
    init_atomic_limit,
    jks_nlie_residual_shifted,
    solve_jks_nlie_newton,
    solve_jks_nlie_newton_from,
    solve_jks_nlie_continuation,
    JKSSolution

@testset "Hubbard1D — JKS Stage C.8 continuation (#523)" begin
    @testset "solve_jks_nlie_newton_from accepts init aux" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        aux_init = init_atomic_limit(grid, 0.01, 4.0, 2.0)
        sol = solve_jks_nlie_newton_from(
            aux_init, grid, 0.01, 4.0, 2.0; alpha=0.5, tol=1e-6, maxiter=10
        )
        @test sol isa JKSSolution
        @test sol.converged
        @test sol.residual < 1e-3
    end

    @testset "Continuation: high-T target trivially converges" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        sol = solve_jks_nlie_continuation(grid, 0.01, 4.0, 2.0; alpha=0.5, tol=1e-6)
        @test sol.converged
        @test sol.residual < 1e-3
    end

    @testset "Continuation reaches mid-T beta=0.1" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        sol = solve_jks_nlie_continuation(
            grid, 0.1, 4.0, 2.0; alpha=0.5, tol=1e-6, beta_start=0.01
        )
        @test isfinite(sol.residual)
    end

    @testset "Continuation: beta_target = beta_start works" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        sol = solve_jks_nlie_continuation(
            grid, 0.01, 4.0, 2.0; alpha=0.5, tol=1e-6, beta_start=0.01
        )
        @test sol.converged
    end

    @testset "Continuation validates inputs" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        @test_throws DomainError solve_jks_nlie_continuation(grid, 0.0, 4.0, 2.0; alpha=0.5)
        @test_throws DomainError solve_jks_nlie_continuation(
            grid, 1.0, 4.0, 2.0; alpha=0.5, beta_start=0.0
        )
        @test_throws ArgumentError solve_jks_nlie_continuation(
            grid, 0.001, 4.0, 2.0; alpha=0.5, beta_start=0.1
        )
        @test_throws DomainError solve_jks_nlie_continuation(
            grid, 1.0, 4.0, 2.0; alpha=0.5, shrink_factor=0.0
        )
        @test_throws DomainError solve_jks_nlie_continuation(
            grid, 1.0, 4.0, 2.0; alpha=0.5, grow_factor=1.0
        )
    end

    @testset "Continuation does not blow up at intermediate beta" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        sol = solve_jks_nlie_continuation(grid, 0.5, 4.0, 2.0; alpha=0.5, tol=1e-6)
        @test isfinite(sol.residual)
        @test sol.residual < 100
    end
end
