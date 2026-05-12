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
    @test QAtlas.fetch(SherringtonKirkpatrick(; J=0.0), CriticalTemperature(), Infinite()) ==
        0.0
    @test QAtlas.fetch(SherringtonKirkpatrick(; J=-1.5), CriticalTemperature(), Infinite()) ==
        0.0
end
