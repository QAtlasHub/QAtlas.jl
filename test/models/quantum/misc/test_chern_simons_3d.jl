# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: ChernSimons3D — Sugawara central charge for SU(N)_k.
#
# Verifies:
#   * Canonical c(SU(2)_k) = 3k/(k+2) for k = 1, 2, 3, 4
#   * c(SU(3)_1) = 2, c(SU(N)_1) = N - 1 (free-fermion realisation)
#   * Large-k limit c → N² - 1 (dim of SU(N))
#   * DomainError on N < 2 or k < 1
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "ChernSimons3D — SU(2)_k central charge equals 3k/(k+2)" begin
    for k in 1:6
        c = QAtlas.fetch(ChernSimons3D(; N=2, k=k), CentralCharge(), Infinite())
        @test c == Rational(3 * k, k + 2)
    end
end

@testset "ChernSimons3D — SU(N)_1 central charge equals N - 1" begin
    for N in 2:6
        c = QAtlas.fetch(ChernSimons3D(; N=N, k=1), CentralCharge(), Infinite())
        @test c == N - 1
    end
end

@testset "ChernSimons3D — Large-k limit c → N² - 1" begin
    N = 3
    for k in (10, 100, 1000)
        c = QAtlas.fetch(ChernSimons3D(; N=N, k=k), CentralCharge(), Infinite())
        @test c < N^2 - 1
        @test float(c) ≈ N^2 - 1 atol = (N^2 - 1) * N / (k + N)   # 1st-order correction
    end
end

@testset "ChernSimons3D — DomainError on invalid (N, k)" begin
    @test_throws DomainError ChernSimons3D(; N=1, k=1)
    @test_throws DomainError ChernSimons3D(; N=2, k=0)
    @test_throws DomainError ChernSimons3D(; N=2, k=-1)
end
