using QAtlas, Test

@testset "LogarithmicCFT — CentralCharge c = 0 (Phase 1)" begin
    c = QAtlas.fetch(LogarithmicCFT(), CentralCharge(), Infinite())
    @test c == 0
    @test c == 0 // 1
    @test c isa Rational
end
