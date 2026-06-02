# test/models/quantum/Hubbard1D/test_hubbard1d_jks_stage_c5.jl
#
# Stage C.5 regression: physical kernel regularization with alpha shift (#523).

using Test
using QAtlas
using QAtlas.Hubbard1DJKSNLIE:
    JKSContourGrid,
    JKSAuxFunctions,
    init_atomic_limit,
    build_kernel_matrix,
    build_kernel_matrix_shifted,
    jks_nlie_residual_shifted,
    solve_jks_nlie_shifted,
    JKSSolution

@testset "Hubbard1D — JKS Stage C.5 kernel alpha-shift + improved solver (#523)" begin
    @testset "Shifted kernel diagonal is O(1)" begin
        # eps=1e-10 build_kernel_matrix gives K[j,j] ~ 1e9; alpha-shifted
        # version with alpha ~ 0.5 should give K[j,j] of order O(1).
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        K_eps = build_kernel_matrix(grid, 1)
        K_shifted = build_kernel_matrix_shifted(grid, 1, 0.5)
        @test abs(K_eps[8, 8]) > 1e6
        @test abs(K_shifted[8, 8]) < 1.0
    end

    @testset "Shifted kernel off-diagonal is O(1)" begin
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        K_shifted = build_kernel_matrix_shifted(grid, 1, 0.5)
        @test abs(K_shifted[5, 10]) < 1.0
    end

    @testset "build_kernel_matrix_shifted validates inputs" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        @test_throws DomainError build_kernel_matrix_shifted(grid, 0, 0.5)
        @test_throws DomainError build_kernel_matrix_shifted(grid, 1, 0.0)
        @test_throws DomainError build_kernel_matrix_shifted(grid, 1, -0.1)
    end

    @testset "Shifted residual scale is much smaller than C.4" begin
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        beta, U, mu, alpha = 1.0, 4.0, 2.0, 0.5
        aux = init_atomic_limit(grid, beta, U, mu)
        res = jks_nlie_residual_shifted(aux, grid, beta, U, mu, alpha)
        @test all(isfinite, res)
        @test maximum(abs.(res)) < 100  # massive improvement from 4e8
    end

    @testset "Shifted residual validation" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        aux = init_atomic_limit(grid, 1.0, 1.0, 0.5)
        @test_throws DomainError jks_nlie_residual_shifted(aux, grid, 0.0, 1.0, 0.5, 0.1)
        @test_throws DomainError jks_nlie_residual_shifted(aux, grid, 1.0, 1.0, 0.5, 0.0)
        aux_wrong = JKSAuxFunctions(4)
        @test_throws DimensionMismatch jks_nlie_residual_shifted(
            aux_wrong, grid, 1.0, 1.0, 0.5, 0.1
        )
    end

    @testset "Shifted solver returns finite residual" begin
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        sol = solve_jks_nlie_shifted(
            grid, 1.0, 4.0, 2.0; alpha=0.5, tol=1e-10, maxiter=100, alpha_mix=0.05
        )
        @test isfinite(sol.residual)
    end

    @testset "Shifted solver validates inputs" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        @test_throws DomainError solve_jks_nlie_shifted(grid, 0.0, 4.0, 2.0; alpha=0.5)
        @test_throws DomainError solve_jks_nlie_shifted(
            grid, 1.0, 4.0, 2.0; alpha=0.5, alpha_mix=0.0
        )
        @test_throws DomainError solve_jks_nlie_shifted(
            grid, 1.0, 4.0, 2.0; alpha=0.5, alpha_mix=2.0
        )
    end
end
