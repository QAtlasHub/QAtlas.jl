# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: SherringtonKirkpatrick — T_c = J.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "SherringtonKirkpatrick — T_c = J" begin
    for J in (0.5, 1.0, 2.5)
        @test QAtlas.fetch(
            SherringtonKirkpatrick(; J=J), CriticalTemperature(), Infinite()
        ) ≈ J
    end
end

@testset "SherringtonKirkpatrick — J ≤ 0 returns 0" begin
    @test QAtlas.fetch(
        SherringtonKirkpatrick(; J=0.0), CriticalTemperature(), Infinite()
    ) == 0.0
    @test QAtlas.fetch(
        SherringtonKirkpatrick(; J=-1.5), CriticalTemperature(), Infinite()
    ) == 0.0
end

@testset "SherringtonKirkpatrick — Parisi T=0 ground-state energy density (Phase 2)" begin
    # Default J=1 → e_0 ≈ -0.7631667
    e0 = QAtlas.fetch(SherringtonKirkpatrick(), Energy{:per_site}(), Infinite())
    @test e0 ≈ -0.7631667
    @test e0 < 0  # ground-state is negative for spin-glass mean-field
    # Scales linearly with J
    e0_3 = QAtlas.fetch(SherringtonKirkpatrick(; J=3.0), Energy{:per_site}(), Infinite())
    @test e0_3 ≈ 3 * e0
    # Identifies the Parisi/full-RSB value within Crisanti-Rizzo error bar
    @test isapprox(e0, -0.7631667; atol=1e-5)
end

@testset "SherringtonKirkpatrick — Energy rejects J ≤ 0 (Phase 2)" begin
    m = SherringtonKirkpatrick(; J=1.0)
    @test_throws DomainError QAtlas.fetch(m, Energy{:per_site}(), Infinite(); J=0.0)
    @test_throws DomainError QAtlas.fetch(m, Energy{:per_site}(), Infinite(); J=-1.0)
end
