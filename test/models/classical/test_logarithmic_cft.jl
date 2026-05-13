using QAtlas, Test

@testset "LogarithmicCFT — CentralCharge c = 0 (Phase 1)" begin
    c = QAtlas.fetch(LogarithmicCFT(), CentralCharge(), Infinite())
    @test c == 0
    @test c == 0 // 1
    @test c isa Rational
    @test c isa Rational{Int}
    @test iszero(c)
    # Idempotent — no hidden state
    @test QAtlas.fetch(LogarithmicCFT(), CentralCharge(), Infinite()) === c
end
