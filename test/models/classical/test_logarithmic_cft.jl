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

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "LogarithmicCFT — verification cards" begin
    # Logarithmic CFT (c=0 / percolation-type): central charge is exactly 0.
    verify(
        LogarithmicCFT(),
        CentralCharge(),
        Infinite();
        route=:second_closed_form,
        independent=0.0,
        agree_within=1e-12,
        refs=["Logarithmic CFT: c = 0 (identity conformal dimension vanishes)"],
    )
end
# ── additional verification cards (#381 batch 3) ─────────────────────────
@testset "LogarithmicCFT — CentralCharge (#381 batch 3)" begin
    # Default LogarithmicCFT registers c = 0 (the canonical c=0 LogCFT
    # examples: percolation, dilute polymers / SAW; cf. Cardy 1999).
    verify(
        LogarithmicCFT(),
        CentralCharge(),
        Infinite();
        route=:second_closed_form,
        independent=0//1,
        agree_within=0,
        refs=["Cardy 1999: canonical c=0 LogCFTs (percolation, SAW)"],
    )
end
