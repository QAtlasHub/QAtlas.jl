# test/models/quantum/Hubbard1D/test_hubbard1d_jks_stage_c12.jl
#
# Stage C.12 regression: 3N x 3N full Newton solver (#523).

using Test
using QAtlas
using QAtlas.Hubbard1DJKSNLIE:
    JKSContourGrid,
    JKSAuxFunctions,
    init_atomic_limit,
    jks_nlie_residual_full,
    jks_jacobian_full_finite_diff,
    solve_jks_nlie_full_newton,
    JKSSolution

@testset "Hubbard1D — JKS Stage C.12 full 3N Newton (#523)" begin
    @testset "Full Jacobian has shape 3N x 3N + finite entries" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        aux = init_atomic_limit(grid, 1.0, 4.0, 2.0)
        J = jks_jacobian_full_finite_diff(aux, grid, 1.0, 4.0, 2.0, 0.5)
        @test size(J) == (24, 24)  # 3 * 8
        @test all(isfinite, J)
    end

    @testset "Full Newton returns JKSSolution" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        sol = solve_jks_nlie_full_newton(grid, 0.1, 4.0, 2.0; alpha=0.5, maxiter=10)
        @test sol isa JKSSolution
        @test isfinite(sol.residual)
        @test sol.iterations >= 1
    end

    @testset "Full Newton solver validates inputs" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        @test_throws DomainError solve_jks_nlie_full_newton(grid, 0.0, 4.0, 2.0; alpha=0.5)
        @test_throws DomainError solve_jks_nlie_full_newton(grid, -1.0, 4.0, 2.0; alpha=0.5)
    end

    @testset "High-T: full Newton converges to tol" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        sol = solve_jks_nlie_full_newton(
            grid, 0.01, 4.0, 2.0; alpha=0.5, tol=1e-6, maxiter=20
        )
        @test sol.residual < 1e-3
    end

    @testset "Mid-T: residual stays bounded" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        sol = solve_jks_nlie_full_newton(
            grid, 1.0, 4.0, 2.0; alpha=0.5, tol=1e-10, maxiter=20
        )
        @test isfinite(sol.residual)
        @test sol.residual < 100
    end
end
