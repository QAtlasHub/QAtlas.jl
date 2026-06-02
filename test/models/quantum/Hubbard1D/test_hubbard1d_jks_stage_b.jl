# test/models/quantum/Hubbard1D/test_hubbard1d_jks_stage_b.jl
#
# Stage B regression for the JKS kernel evaluator + contour grid +
# discrete convolution scaffold (#523).

using Test
using QAtlas
using QAtlas.Hubbard1DJKSNLIE:
    jks_kernel_K_n_concrete, JKSContourGrid, build_kernel_matrix, apply_kernel

@testset "Hubbard1D — JKS Stage B kernels (#523)" begin
    @testset "Kernel value at simple points (gamma = 2pi/3, n = 1)" begin
        gamma = 2pi/3
        # K_1(1.0) = gamma / (pi * 1 * (1 + 2i gamma)) — direct evaluation
        s = 1.0 + 0im
        K = jks_kernel_K_n_concrete(s, 1, gamma)
        expected = gamma / (pi * s * (s + 2im * gamma))
        @test isapprox(K, expected; atol=1e-14)
    end

    @testset "Kernel decays at infinity: |K_n(s)| ~ gamma / (pi |s|^2)" begin
        gamma = 2pi/3
        for n in (1, 2, 3)
            K_large = jks_kernel_K_n_concrete(100.0 + 0im, n, gamma)
            # |K_n(s)| ~ gamma / (pi s^2)
            @test abs(K_large) < gamma / (pi * 100.0^2) * 2
            @test abs(K_large) > gamma / (pi * 100.0^2) / 2
        end
    end

    @testset "Conjugation: K_n(-s) = conj(K_n(s)) for real s" begin
        # K_n(s) = gamma / (pi s (s + 2niγ))
        # K_n(-s) = gamma / (pi s (s - 2niγ)) = conj(K_n(s)) for real s.
        gamma = pi/2
        for s in (0.5, 1.0, 3.7)
            K_pos = jks_kernel_K_n_concrete(s + 0im, 1, gamma)
            K_neg = jks_kernel_K_n_concrete(-s + 0im, 1, gamma)
            @test isapprox(K_neg, conj(K_pos); atol=1e-14)
        end
    end

    @testset "Kernel pole at s = 0 returns Inf" begin
        @test isinf(real(jks_kernel_K_n_concrete(0.0 + 0im, 1, pi/2)))
    end

    @testset "Kernel input validation" begin
        @test_throws DomainError jks_kernel_K_n_concrete(1.0 + 0im, 0, pi/2)
        @test_throws DomainError jks_kernel_K_n_concrete(1.0 + 0im, -1, pi/2)
        @test_throws DomainError jks_kernel_K_n_concrete(1.0 + 0im, 1, 0.0)
        @test_throws DomainError jks_kernel_K_n_concrete(1.0 + 0im, 1, pi)
    end

    @testset "JKSContourGrid construction + invariants" begin
        g = JKSContourGrid(64, 2pi/3; x_max=5.0)
        @test g.N == 64
        @test isapprox(g.gamma, 2pi/3)
        @test isapprox(g.x_max, 5.0)
        @test length(g.x) == 64
        @test isapprox(first(g.x), -5.0; atol=1e-14)
        @test isapprox(last(g.x), 5.0; atol=1e-14)
        @test isapprox(g.dx, g.x[2] - g.x[1]; atol=1e-14)

        # Default x_max = 10
        g_default = JKSContourGrid(8, pi/2)
        @test isapprox(g_default.x_max, 10.0)
    end

    @testset "JKSContourGrid input validation" begin
        @test_throws DomainError JKSContourGrid(1, pi/2)
        @test_throws DomainError JKSContourGrid(64, 0.0)
        @test_throws DomainError JKSContourGrid(64, pi)
        @test_throws DomainError JKSContourGrid(64, pi/2; x_max=0.0)
        @test_throws DomainError JKSContourGrid(64, pi/2; x_max=-1.0)
    end

    @testset "build_kernel_matrix dimensions + finiteness" begin
        g = JKSContourGrid(16, 2pi/3; x_max=3.0)
        K1 = build_kernel_matrix(g, 1)
        @test size(K1) == (16, 16)
        @test eltype(K1) == ComplexF64
        @test all(isfinite, K1)
    end

    @testset "build_kernel_matrix is Toeplitz (depends on j - k only)" begin
        # K[j, k] depends only on x_j - x_k = (j - k) * dx, so each
        # anti-diagonal should be constant.
        g = JKSContourGrid(8, pi/2; x_max=2.0)
        K = build_kernel_matrix(g, 1)
        for diag_offset in -7:7
            vals = [K[j, j - diag_offset] for j in 1:8 if 1 <= j - diag_offset <= 8]
            length(vals) > 1 || continue
            @test all(isapprox(v, vals[1]; atol=1e-12) for v in vals)
        end
    end

    @testset "apply_kernel dimension check + linearity" begin
        g = JKSContourGrid(8, pi/2; x_max=2.0)
        K = build_kernel_matrix(g, 1)
        f = ComplexF64.(g.x)  # x as test vector
        out = apply_kernel(K, f)
        @test length(out) == 8
        @test eltype(out) == ComplexF64

        # Linearity check: K(2f) = 2 K(f)
        out_2f = apply_kernel(K, 2 .* f)
        @test all(isapprox(out_2f[j], 2 * out[j]; atol=1e-12) for j in 1:8)

        # Dimension mismatch raises
        @test_throws DimensionMismatch apply_kernel(K, ComplexF64[1.0])
    end

    @testset "build_kernel_matrix eps validation" begin
        g = JKSContourGrid(8, pi/2; x_max=2.0)
        @test_throws DomainError build_kernel_matrix(g, 1; eps=0.0)
        @test_throws DomainError build_kernel_matrix(g, 1; eps=-1e-10)
    end
end
