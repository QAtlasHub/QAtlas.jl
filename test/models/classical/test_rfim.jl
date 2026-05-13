# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: RFIM — Imry-Ma T_c = 0 at d ≤ 2.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "RFIM — Imry-Ma T_c = 0 at d = 1, 2" begin
    for d in (1, 2)
        @test QAtlas.fetch(RFIM(; Δ=0.5), CriticalTemperature(), Infinite(); d=d) == 0.0
        @test QAtlas.fetch(RFIM(; J=2.0, Δ=1.5), CriticalTemperature(), Infinite(); d=d) ==
            0.0
    end
end

@testset "RFIM — DomainError on d ≥ 3 (Phase 2 numerical reference)" begin
    for d in (3, 4, 5)
        @test_throws DomainError QAtlas.fetch(
            RFIM(; Δ=1.0), CriticalTemperature(), Infinite(); d=d
        )
    end
end

@testset "RFIM — DomainError on Δ = 0 (use pure Ising)" begin
    @test_throws DomainError QAtlas.fetch(
        RFIM(; Δ=0.0), CriticalTemperature(), Infinite(); d=2
    )
end

@testset "RFIM — DomainError on d < 1" begin
    @test_throws DomainError QAtlas.fetch(
        RFIM(; Δ=1.0), CriticalTemperature(), Infinite(); d=0
    )
end
