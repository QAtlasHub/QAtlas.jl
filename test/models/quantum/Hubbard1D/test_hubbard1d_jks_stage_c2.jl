# test/models/quantum/Hubbard1D/test_hubbard1d_jks_stage_c2.jl
#
# Stage C.2 regression for the JKS elementary function phi + free-energy
# evaluator (#523).

using Test
using QAtlas
using QAtlas.Hubbard1DJKSNLIE:
    jks_log_phi,
    jks_phi,
    jks_log_z_deriv,
    free_energy_jks,
    JKSAuxFunctions,
    JKSContourGrid,
    init_atomic_limit,
    atomic_free_energy

@testset "Hubbard1D — JKS Stage C.2 phi + FE evaluator (#523)" begin
    @testset "jks_log_phi formula on |s| > 1" begin
        # ln phi(s) = -2 beta |s| sqrt(1 - 1/s^2)
        for s in (1.5, 2.0, 5.0, -3.0), beta in (0.1, 1.0, 5.0)
            v = jks_log_phi(s, beta)
            expected = -2 * beta * abs(s) * sqrt(1 - 1/s^2)
            @test isapprox(v, expected; atol=1e-14)
        end
    end

    @testset "jks_log_phi returns -Inf on the cut |s| <= 1" begin
        @test jks_log_phi(0.5, 1.0) == -Inf
        @test jks_log_phi(-0.99, 1.0) == -Inf
        @test jks_log_phi(0.0, 1.0) == -Inf
        @test jks_log_phi(1.0, 1.0) == -Inf
    end

    @testset "jks_phi = exp(jks_log_phi)" begin
        for s in (1.5, 3.0), beta in (0.1, 1.0)
            @test jks_phi(s, beta) ≈ exp(jks_log_phi(s, beta))
        end
        @test jks_phi(0.5, 1.0) == 0.0  # exp(-Inf)
    end

    @testset "jks_log_phi validates beta > 0" begin
        @test_throws DomainError jks_log_phi(2.0, 0.0)
        @test_throws DomainError jks_log_phi(2.0, -1.0)
    end

    @testset "jks_log_z_deriv = 1 / (s * sqrt(1 - 1/s^2))  (Stage C.24 paper-precise eq 23)" begin
        # Real axis, |s| > 1: real value with sign(s) factor.
        # Old broken form 1/sqrt(s^2-1) was missing this factor on s < 0;
        # Stage C.24 fixed jks_log_z_deriv per paper eq (23).
        for s in (1.5, 2.0, 5.0)
            v = jks_log_z_deriv(s + 0im)
            expected = 1 / (s * sqrt(1 - 1/s^2 + 0im))
            @test isapprox(v, expected; atol=1e-14)
        end
        # s = -3: with sign correction, value is -1/sqrt(8) (was +1/sqrt(8)
        # in the pre-C.24 broken implementation).
        v_neg = jks_log_z_deriv(-3.0 + 0im)
        @test isapprox(v_neg, -1 / sqrt(8); atol=1e-14)
        # Inside the cut |s| < 1: imaginary, |value| = 1 / sqrt(1 - s^2)
        v_inside = jks_log_z_deriv(0.5 + 0im)
        @test abs(real(v_inside)) < 1e-14
        @test isapprox(abs(imag(v_inside)), 1 / sqrt(1 - 0.25); atol=1e-14)
    end

    @testset "Free energy evaluator returns finite real number" begin
        grid = JKSContourGrid(64, 1.0; x_max=2.0)
        beta, U, mu = 1.0, 4.0, 2.0
        aux = init_atomic_limit(grid, beta, U, mu)
        f_jks = free_energy_jks(aux, grid, beta, U; mu=mu)
        @test isfinite(f_jks)
        @test f_jks isa Float64
    end

    @testset "Free energy at U = 0 is negative (atomic-limit baseline)" begin
        grid = JKSContourGrid(64, 1e-3; x_max=1.5)  # eta -> 0 limit
        beta, U, mu = 1.0, 0.0, 0.0
        aux = init_atomic_limit(grid, beta, U, mu)
        f_jks = free_energy_jks(aux, grid, beta, U; mu=mu)
        f_atomic = atomic_free_energy(beta, U, mu)
        @test isfinite(f_jks)
        @test f_atomic ≈ -log(4) atol=1e-12
    end

    @testset "Free energy validation: beta > 0 required" begin
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        aux = init_atomic_limit(grid, 1.0, 1.0, 0.5)
        @test_throws DomainError free_energy_jks(aux, grid, 0.0, 1.0; mu=0.5)
        @test_throws DomainError free_energy_jks(aux, grid, -1.0, 1.0; mu=0.5)
    end

    @testset "Free energy: aux/grid length mismatch raises" begin
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        aux_wrong = JKSAuxFunctions(8)  # not 16
        @test_throws DimensionMismatch free_energy_jks(aux_wrong, grid, 1.0, 1.0; mu=0.5)
    end
end
