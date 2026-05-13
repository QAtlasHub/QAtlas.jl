# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: SYK — large-N IR Majorana ConformalWeight Δ_ψ = 1/q.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "SYK — IR Majorana ConformalWeight Δ_ψ = 1/q (Phase 1)" begin
    # Default q=4
    @test QAtlas.fetch(SYK(), ConformalWeights(), Infinite(); field=:ψ) == 1 // 4
    # Other even q
    @test QAtlas.fetch(SYK(; q=2), ConformalWeights(), Infinite(); field=:ψ) == 1 // 2
    @test QAtlas.fetch(SYK(; q=6), ConformalWeights(), Infinite(); field=:ψ) == 1 // 6
    @test QAtlas.fetch(SYK(; q=8), ConformalWeights(), Infinite(); field=:ψ) == 1 // 8
    # Default field = :ψ
    @test QAtlas.fetch(SYK(), ConformalWeights(), Infinite()) == 1 // 4
end

@testset "SYK — unitarity bound and monotonicity (Phase 1 cross-check)" begin
    # Unitarity bound: free Majorana in d=1 has Δ=1/2; q=2 saturates, q>2 strict.
    # Equivalently 1//q ≤ 1//2 for all valid q ≥ 2.
    for q in (2, 4, 6, 8, 10, 12)
        Δ = QAtlas.fetch(SYK(; q=q), ConformalWeights(), Infinite(); field=:ψ)
        @test Δ ≤ 1 // 2
    end
    # Monotonicity: Δ_ψ(q+2) < Δ_ψ(q) — strictly decreasing in q.
    for q in (2, 4, 6, 8, 10)
        Δq = QAtlas.fetch(SYK(; q=q), ConformalWeights(), Infinite(); field=:ψ)
        Δqp = QAtlas.fetch(SYK(; q=q + 2), ConformalWeights(), Infinite(); field=:ψ)
        @test Δqp < Δq
    end
end

@testset "SYK — rejects q ≤ 1 / q odd (Phase 1)" begin
    @test_throws DomainError SYK(; q=1)
    @test_throws DomainError SYK(; q=0)
    @test_throws DomainError SYK(; q=3)  # odd q invalid for Majorana
    @test_throws DomainError SYK(; q=5)
end

@testset "SYK — rejects unknown field (Phase 2 deferral)" begin
    m = SYK()
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); field=:O)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); field=:ε)
end
