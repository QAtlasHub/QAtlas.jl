# test/models/quantum/Hubbard1D/test_hubbard1d_jks_stage_c7.jl
#
# Stage C.7 regression: damped Newton solver (#523).

using Test
using QAtlas
using QAtlas.Hubbard1DJKSNLIE:
    JKSContourGrid,
    JKSAuxFunctions,
    init_atomic_limit,
    jks_nlie_residual_shifted,
    jks_jacobian_b_finite_diff,
    solve_jks_nlie_newton,
    JKSSolution

@testset "Hubbard1D — JKS Stage C.7 Newton solver (#523)" begin
    @testset "Jacobian has correct shape + finite entries" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        aux = init_atomic_limit(grid, 1.0, 4.0, 2.0)
        J = jks_jacobian_b_finite_diff(aux, grid, 1.0, 4.0, 2.0, 0.5)
        @test size(J) == (8, 8)
        @test all(isfinite, J)
    end

    @testset "Newton solver returns JKSSolution" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        sol = solve_jks_nlie_newton(grid, 0.5, 4.0, 2.0; alpha=0.5, maxiter=10)
        @test sol isa JKSSolution
        @test isfinite(sol.residual)
        @test sol.iterations >= 1
    end

    @testset "Newton solver validates inputs" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        @test_throws DomainError solve_jks_nlie_newton(grid, 0.0, 4.0, 2.0; alpha=0.5)
        @test_throws DomainError solve_jks_nlie_newton(grid, -1.0, 4.0, 2.0; alpha=0.5)
    end

    @testset "Newton: high-T regime converges quickly" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        sol = solve_jks_nlie_newton(grid, 0.01, 4.0, 2.0; alpha=0.5, tol=1e-6, maxiter=30)
        @test sol.residual < 1e-3
    end

    @testset "Newton at mid-T: residual stays bounded" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        sol = solve_jks_nlie_newton(grid, 1.0, 4.0, 2.0; alpha=0.5, tol=1e-10, maxiter=30)
        @test isfinite(sol.residual)
        @test sol.residual < 100
    end
end
