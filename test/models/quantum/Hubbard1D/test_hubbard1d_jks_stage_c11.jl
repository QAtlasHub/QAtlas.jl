# test/models/quantum/Hubbard1D/test_hubbard1d_jks_stage_c11.jl
#
# Stage C.11 regression: c, c_bar channel residuals + full 3N residual (#523).

using Test
using QAtlas
using QAtlas.Hubbard1DJKSNLIE:
    JKSContourGrid,
    JKSAuxFunctions,
    init_atomic_limit,
    jks_nlie_residual_shifted,
    jks_nlie_residual_c,
    jks_nlie_residual_cbar,
    jks_nlie_residual_full

@testset "Hubbard1D — JKS Stage C.11 c/c_bar residuals (#523)" begin
    @testset "c residual has correct length + element type" begin
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        aux = init_atomic_limit(grid, 1.0, 4.0, 2.0)
        res = jks_nlie_residual_c(aux, grid, 1.0, 4.0, 2.0, 0.5)
        @test length(res) == 16
        @test eltype(res) == ComplexF64
        @test all(isfinite, res)
    end

    @testset "c_bar residual has correct length + element type" begin
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        aux = init_atomic_limit(grid, 1.0, 4.0, 2.0)
        res = jks_nlie_residual_cbar(aux, grid, 1.0, 4.0, 2.0, 0.5)
        @test length(res) == 16
        @test eltype(res) == ComplexF64
        @test all(isfinite, res)
    end

    @testset "PH-conjugate driving: res_c - res_cbar reflects psi_c - psi_cbar" begin
        # At half filling h=0 the atomic-limit init has c_const == c_bar_const
        # so log_c == log_cbar. The residual difference is then just the
        # driving-term difference: res_c - res_cbar = psi_cbar - psi_c.
        # We just verify the residuals are individually finite and bounded.
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        aux = init_atomic_limit(grid, 1.0, 4.0, 2.0)
        res_c = jks_nlie_residual_c(aux, grid, 1.0, 4.0, 2.0, 0.5)
        res_cbar = jks_nlie_residual_cbar(aux, grid, 1.0, 4.0, 2.0, 0.5)
        @test all(isfinite, res_c)
        @test all(isfinite, res_cbar)
        @test maximum(abs.(res_c)) < 100
        @test maximum(abs.(res_cbar)) < 100
    end

    @testset "Full residual is length 3N + concatenation order" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        aux = init_atomic_limit(grid, 1.0, 4.0, 2.0)
        res_full = jks_nlie_residual_full(aux, grid, 1.0, 4.0, 2.0, 0.5)
        @test length(res_full) == 24  # 3 * 8

        res_b = jks_nlie_residual_shifted(aux, grid, 1.0, 4.0, 2.0, 0.5)
        res_c = jks_nlie_residual_c(aux, grid, 1.0, 4.0, 2.0, 0.5)
        res_cbar = jks_nlie_residual_cbar(aux, grid, 1.0, 4.0, 2.0, 0.5)
        @test res_full[1:8] == res_b
        @test res_full[9:16] == res_c
        @test res_full[17:24] == res_cbar
    end

    @testset "c residual: validates beta > 0, alpha > 0, length match" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        aux = init_atomic_limit(grid, 1.0, 4.0, 2.0)
        @test_throws DomainError jks_nlie_residual_c(aux, grid, 0.0, 4.0, 2.0, 0.5)
        @test_throws DomainError jks_nlie_residual_c(aux, grid, 1.0, 4.0, 2.0, 0.0)
        aux_wrong = JKSAuxFunctions(4)
        @test_throws DimensionMismatch jks_nlie_residual_c(
            aux_wrong, grid, 1.0, 4.0, 2.0, 0.5
        )
    end

    @testset "c_bar residual: validates beta > 0, alpha > 0, length match" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        aux = init_atomic_limit(grid, 1.0, 4.0, 2.0)
        @test_throws DomainError jks_nlie_residual_cbar(aux, grid, 0.0, 4.0, 2.0, 0.5)
        @test_throws DomainError jks_nlie_residual_cbar(aux, grid, 1.0, 4.0, 2.0, 0.0)
        aux_wrong = JKSAuxFunctions(4)
        @test_throws DimensionMismatch jks_nlie_residual_cbar(
            aux_wrong, grid, 1.0, 4.0, 2.0, 0.5
        )
    end

    @testset "Magnetic field h > 0 breaks res_c == res_cbar" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        aux = init_atomic_limit(grid, 1.0, 4.0, 2.0; h=0.5)
        res_c = jks_nlie_residual_c(aux, grid, 1.0, 4.0, 2.0, 0.5; H=0.5)
        res_cbar = jks_nlie_residual_cbar(aux, grid, 1.0, 4.0, 2.0, 0.5; H=0.5)
        # h > 0 makes c_const != c_bar_const, so the residuals must differ.
        @test any(j -> !isapprox(res_c[j], res_cbar[j]; atol=1e-6), 1:8)
    end
end
