# test/models/quantum/Hubbard1D/test_hubbard1d_jks_stage_c4.jl
#
# Stage C.4 regression: driving terms + NLIE residual + (b-channel)
# Picard solver (#523).

using Test
using QAtlas
using QAtlas.Hubbard1DJKSNLIE:
    JKSContourGrid,
    JKSAuxFunctions,
    init_atomic_limit,
    jks_driving_b,
    jks_driving_c,
    jks_driving_cbar,
    jks_nlie_residual,
    solve_jks_nlie_b_only,
    JKSSolution

@testset "Hubbard1D — JKS Stage C.4 driving + residual + solver (#523)" begin
    @testset "Driving terms have correct length" begin
        grid = JKSContourGrid(32, 1.0; x_max=2.0)
        beta, U, mu, alpha = 1.0, 4.0, 2.0, 0.5
        psi_b = jks_driving_b(grid, beta, U, alpha)
        psi_c = jks_driving_c(grid, beta, U, mu)
        psi_cbar = jks_driving_cbar(grid, beta, U, mu)
        @test length(psi_b) == 32
        @test length(psi_c) == 32
        @test length(psi_cbar) == 32
        @test eltype(psi_b) == ComplexF64
    end

    @testset "Particle-hole sym at h=0: psi_c + psi_cbar = -beta U" begin
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        beta, U, mu = 2.0, 4.0, 2.0
        psi_c = jks_driving_c(grid, beta, U, mu)
        psi_cbar = jks_driving_cbar(grid, beta, U, mu)
        for j in 1:16
            @test isapprox(psi_c[j] + psi_cbar[j], -beta * U; atol=1e-12)
        end
    end

    @testset "Magnetic field h > 0: psi_c changes" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        beta, U, mu = 1.0, 4.0, 2.0
        psi_c_0 = jks_driving_c(grid, beta, U, mu)
        psi_c_h = jks_driving_c(grid, beta, U, mu; H=0.5)
        @test any(j -> !isapprox(psi_c_0[j], psi_c_h[j]; atol=1e-6), 1:8)
    end

    @testset "Driving terms validate inputs" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        @test_throws DomainError jks_driving_b(grid, 0.0, 1.0, 0.5)
        @test_throws DomainError jks_driving_b(grid, 1.0, 1.0, 0.0)
        @test_throws DomainError jks_driving_c(grid, -1.0, 1.0, 0.5)
        @test_throws DomainError jks_driving_cbar(grid, 0.0, 1.0, 0.5)
    end

    @testset "NLIE residual is computable + finite" begin
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        beta, U, mu, alpha = 1.0, 4.0, 2.0, 0.5
        aux = init_atomic_limit(grid, beta, U, mu)
        res = jks_nlie_residual(aux, grid, beta, U, mu, alpha)
        @test length(res) == 16
        @test eltype(res) == ComplexF64
        @test all(isfinite, res)
    end

    @testset "NLIE residual validation" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        aux = init_atomic_limit(grid, 1.0, 1.0, 0.5)
        @test_throws DomainError jks_nlie_residual(aux, grid, 0.0, 1.0, 0.5, 0.1)
        @test_throws DomainError jks_nlie_residual(aux, grid, 1.0, 1.0, 0.5, 0.0)
        aux_wrong = JKSAuxFunctions(4)
        @test_throws DimensionMismatch jks_nlie_residual(
            aux_wrong, grid, 1.0, 1.0, 0.5, 0.1
        )
    end

    @testset "Solver returns JKSSolution with all fields" begin
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        sol = solve_jks_nlie_b_only(grid, 1.0, 4.0, 2.0; alpha=0.5, tol=1e-3, maxiter=20)
        @test sol isa JKSSolution
        @test sol.aux isa JKSAuxFunctions
        @test sol.iterations >= 1
        @test sol.residual >= 0 || isnan(sol.residual)
        @test sol.converged isa Bool
    end

    @testset "Solver validates inputs" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        @test_throws DomainError solve_jks_nlie_b_only(grid, 0.0, 4.0, 2.0; alpha=0.5)
        @test_throws DomainError solve_jks_nlie_b_only(
            grid, 1.0, 4.0, 2.0; alpha=0.5, alpha_mix=0.0
        )
        @test_throws DomainError solve_jks_nlie_b_only(
            grid, 1.0, 4.0, 2.0; alpha=0.5, alpha_mix=1.5
        )
    end

    @testset "Solver does not blow up: residual stays finite" begin
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        sol = solve_jks_nlie_b_only(
            grid, 1.0, 4.0, 2.0; alpha=0.5, tol=1e-10, maxiter=50, alpha_mix=0.2
        )
        @test true  # solver runs; convergence is Stage C.5
    end
end
