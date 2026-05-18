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

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "RFIM — verification cards" begin
    # Imry-Ma 1975: the lower critical dimension of the RFIM is 2, so
    # for d <= 2 there is no finite-temperature order => Tc = 0.
    for d in (1, 2)
        verify(
            RFIM(; Δ=0.5),
            CriticalTemperature(),
            Infinite();
            route=:second_closed_form,
            fetch_kw=(; d=d),
            independent=0.0,
            agree_within=1e-12,
            refs=["Imry-Ma 1975: RFIM lower critical dimension is 2 => Tc = 0 for d ≤ 2"],
        )
    end
end
