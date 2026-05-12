# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: GrossNeveu — UV c = N at g = 0.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "GrossNeveu — UV free-fermion c = N at g = 0" begin
    for N in 1:5
        @test QAtlas.fetch(GrossNeveu(; N=N, g=0.0), CentralCharge(), Infinite()) == N
    end
end

@testset "GrossNeveu — DomainError on N < 1" begin
    @test_throws DomainError GrossNeveu(; N=0, g=0.0)
end

@testset "GrossNeveu — DomainError on g ≠ 0 (Phase 2)" begin
    @test_throws DomainError QAtlas.fetch(
        GrossNeveu(; N=2, g=0.5), CentralCharge(), Infinite()
    )
end
