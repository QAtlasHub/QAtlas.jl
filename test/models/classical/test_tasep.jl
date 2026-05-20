using QAtlas, Test

@testset "TASEP — SteadyStateCurrent j(ρ) = p ρ (1−ρ) (Phase 1)" begin
    # Default: p=1, ρ=0.5 → j_max = 1/4
    @test QAtlas.fetch(TASEP(), SteadyStateCurrent(), Infinite()) == 0.25
    # Symmetric in ρ ↔ 1−ρ (particle-hole)
    for ρ in (0.1, 0.3, 0.7, 0.9)
        @test QAtlas.fetch(TASEP(; ρ=ρ), SteadyStateCurrent(), Infinite()) ≈
            QAtlas.fetch(TASEP(; ρ=1-ρ), SteadyStateCurrent(), Infinite())
    end
    # Boundary: ρ = 0 or 1 → j = 0
    @test QAtlas.fetch(TASEP(; ρ=0.0), SteadyStateCurrent(), Infinite()) == 0.0
    @test QAtlas.fetch(TASEP(; ρ=1.0), SteadyStateCurrent(), Infinite()) == 0.0
    # Linear in p
    @test QAtlas.fetch(TASEP(; p=3.0, ρ=0.5), SteadyStateCurrent(), Infinite()) == 0.75
    # Closed-form at ρ = 0.3 (j = p·0.3·0.7 = 0.21)
    @test QAtlas.fetch(TASEP(; ρ=0.3), SteadyStateCurrent(), Infinite()) ≈ 0.21
end

@testset "TASEP — rejects p ≤ 0 or ρ ∉ [0,1] (Phase 1)" begin
    @test_throws DomainError TASEP(; p=0.0)
    @test_throws DomainError TASEP(; p=-1.0)
    @test_throws DomainError TASEP(; ρ=-0.1)
    @test_throws DomainError TASEP(; ρ=1.5)
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "TASEP — verification cards" begin
    # Fundamental diagram j(ρ) = p ρ (1 - ρ) (independent closed form)
    for (p, ρ) in ((1.0, 0.5), (1.0, 0.3), (0.7, 0.5), (1.0, 0.9))
        verify(
            TASEP(; p=p, ρ=ρ),
            SteadyStateCurrent(),
            Infinite();
            route=:second_closed_form,
            independent=p * ρ * (1 - ρ),
            agree_within=1e-12,
            refs=["TASEP mean-field steady state current j = p ρ (1-ρ)"],
        )
    end

    # Particle-hole boundary: j(0) = j(1) = 0
    verify(
        TASEP(; p=1.0, ρ=0.0),
        SteadyStateCurrent(),
        Infinite();
        route=:limiting_case,
        independent=0.0,
        agree_within=1e-12,
        refs=["Empty/full lattice carries no current: j(0)=j(1)=0"],
    )
end
# ── additional verification cards (#381 batch 2) ─────────────────────────
@testset "TASEP — additional cards (#381 batch 2)" begin
    # Mean-field steady-state current of TASEP at density ρ on the
    # half-filled ring is J = ρ(1-ρ), max J = 1/4 at ρ = 1/2 (Derrida 1998).
    # NOTE: TASEP() defaults to (p=1.0, ρ=0.5) per the struct constructor
    # (src/models/classical/TASEP/TASEP.jl), i.e. the half-filled maximal-
    # current point — exactly what this card is checking.
    verify(
        TASEP(),
        SteadyStateCurrent(),
        Infinite();
        route=:second_closed_form,
        independent=0.25,
        agree_within=1e-12,
        refs=["Derrida 1998: TASEP steady-state J = ρ(1-ρ), max at ρ=1/2 ⇒ J = 1/4"],
    )
end

