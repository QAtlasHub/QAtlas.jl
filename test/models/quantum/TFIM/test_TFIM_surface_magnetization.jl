# Peschel's exact surface magnetization, summed in closed form on the uniform
# chain.  Every check here is against something derived elsewhere: the sum the
# closed form replaces, the two exponents its limits must reproduce, and the
# 2^N ground state itself.

using Test
using QAtlas
using QAtlas: fetch
using AbstractQAtlas: SurfaceMagnetization, OBC
using LinearAlgebra: eigen, Symmetric

include(joinpath(@__DIR__, "..", "..", "..", "util", "verify.jl"))

# The sum the closed form replaces: Eq. (4.4) evaluated term by term.
_ms_sum(r, N) = (1 + sum(r^(2i) for i in 1:(N - 1); init=0.0))^-0.5

# The 2^N ground state, with no free-fermion step anywhere in it: fix the far
# end by removing its transverse field, then read |<0|σᶻ₁|1>| off the two-fold
# ground space.  H = -J Σ σᶻσᶻ - h Σ σˣ, QAtlas's convention.
function _ms_exact_diagonalization(J, h, N)
    D = 1 << N
    idx = 0:(D - 1)
    sz = [1 - 2 * ((i >> (j - 1)) & 1) for i in idx, j in 1:N]
    H = zeros(Float64, D, D)
    for i in 1:D
        H[i, i] = -J * sum(sz[i, j] * sz[i, j + 1] for j in 1:(N - 1))
    end
    for j in 1:(N - 1)                       # h_N = 0 fixes the far end
        for i in idx
            H[i + 1, (i ⊻ (1 << (j - 1))) + 1] -= h
        end
    end
    F = eigen(Symmetric(H))
    G = F.vectors[:, 1:2]
    sz1 = [1 - 2 * (i & 1) for i in idx]
    return maximum(abs, eigen(Symmetric(G' * (sz1 .* G))).values)
end

@testset "TFIM :: the closed form is the sum it replaces" begin
    for r in (0.3, 0.5, 0.9, 1.0, 1.1, 2.0, 3.0), N in (2, 5, 20, 200)
        m = fetch(TFIM(; J=1.0, h=r), SurfaceMagnetization(), OBC(N))
        @test m ≈ _ms_sum(r, N) rtol = 1e-12
    end
    # N = 1 has no bond, so the sum is empty and the fixed end IS the surface.
    @test fetch(TFIM(), SurfaceMagnetization(), OBC(1)) == 1.0
    @test_throws ArgumentError fetch(TFIM(), SurfaceMagnetization(), OBC(0))
    @test_throws ArgumentError fetch(TFIM(; J=0.0, h=1.0), SurfaceMagnetization(), OBC(4))
end

@testset "TFIM :: against the 2^N ground state, no free fermions involved" begin
    # The closed form comes from a Jordan-Wigner solution; exact diagonalization
    # shares none of that, so agreement is evidence and not a shared convention.
    for (J, h, N) in ((1.0, 0.5, 10), (1.0, 1.0, 10), (1.0, 1.5, 10), (1.0, 3.0, 8))
        @test fetch(TFIM(; J=J, h=h), SurfaceMagnetization(), OBC(N)) ≈
            _ms_exact_diagonalization(J, h, N) rtol = 1e-10
    end
end

@testset "TFIM :: the limits are the surface exponents, not a fit" begin
    # Ordered side, N -> infinity: m_s -> sqrt(1 - r^2), i.e. beta_s = 1/2.
    for r in (0.3, 0.5, 0.9)
        @test fetch(TFIM(; J=1.0, h=r), SurfaceMagnetization(), OBC(4000)) ≈ sqrt(1 - r^2) rtol =
            1e-12
    end
    # At criticality m_s = N^{-1/2} exactly, i.e. x_m^s = 1/2. Both are the 2D
    # Ising surface entries, and both are what the RANDOM chain has for x_m^s by
    # an unrelated argument, so this average cannot tell clean from infinitely
    # random; only the typical value can.
    for N in (4, 16, 64, 1024)
        @test fetch(TFIM(), SurfaceMagnetization(), OBC(N)) == 1 / sqrt(N)
    end
    slope(N1, N2) =
        log(
            fetch(TFIM(), SurfaceMagnetization(), OBC(N2)) /
            fetch(TFIM(), SurfaceMagnetization(), OBC(N1)),
        ) / log(N2 / N1)
    @test slope(64, 1024) ≈ -1 / 2
end

@testset "TFIM :: the disordered side stays representable, and shrinks" begin
    # r > 1 makes m_s exponentially small. Overflowing the geometric sum would
    # report 0; the log form keeps the number until the exponent runs out.
    m = [fetch(TFIM(; J=1.0, h=2.0), SurfaceMagnetization(), OBC(N)) for N in (10, 50, 200)]
    @test all(>(0), m)
    @test issorted(m; rev=true)
    # ln m_s ~ -N ln r, so the ratio of successive logs is the ratio of the Ns.
    @test log(m[3]) / log(m[2]) ≈ 200 / 50 rtol = 0.02
    @test fetch(TFIM(; J=1.0, h=2.0), SurfaceMagnetization(), OBC(5000)) == 0.0   # underflow, said plainly
end

@testset "TFIM :: verify card" begin
    verify(
        TFIM(),
        SurfaceMagnetization(),
        OBC(10);
        route=:ed_finite_size,
        independent=_ms_exact_diagonalization(1.0, 1.0, 10),
        agree_within=1e-10,
        refs=["Peschel1984"],
    )
end
