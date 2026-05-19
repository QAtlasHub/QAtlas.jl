using QAtlas, Test

@testset "LongRangeXY1D — α=Inf NN XX limit MassGap (Phase 1)" begin
    @test QAtlas.fetch(LongRangeXY1D(), MassGap(), Infinite()) == 0.0        # h=0 gapless XX
    @test QAtlas.fetch(LongRangeXY1D(; h=1.0), MassGap(), Infinite()) == 0.0  # inside band
    @test QAtlas.fetch(LongRangeXY1D(; h=2.0), MassGap(), Infinite()) == 0.0  # critical
    @test QAtlas.fetch(LongRangeXY1D(; h=3.0), MassGap(), Infinite()) == 2.0
    @test QAtlas.fetch(LongRangeXY1D(; h=-3.0), MassGap(), Infinite()) == 2.0
    @test QAtlas.fetch(LongRangeXY1D(; J=0.5, h=2.0), MassGap(), Infinite()) == 2.0
end

@testset "LongRangeXY1D — finite α throws DomainError" begin
    @test_throws DomainError QAtlas.fetch(LongRangeXY1D(; α=2.0), MassGap(), Infinite())
    @test_throws DomainError QAtlas.fetch(LongRangeXY1D(; α=0.5), MassGap(), Infinite())
end

@testset "LongRangeXY1D — rejects J ≤ 0 or α ≤ 0" begin
    @test_throws DomainError LongRangeXY1D(; J=0.0)
    @test_throws DomainError LongRangeXY1D(; J=-1.0)
    @test_throws DomainError LongRangeXY1D(; α=0.0)
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "LongRangeXY1D — verification cards" begin
    # α = ∞ NN XX chain is gapless for |h| <= 2J (free-fermion band).
    for h in (0.0, 1.0, 2.0)
        verify(
            LongRangeXY1D(; h=h),
            MassGap(),
            Infinite();
            route=:second_closed_form,
            independent=0.0,
            agree_within=1e-10,
            refs=["NN XX chain: gapless for |h| <= 2J (free fermion)"],
        )
    end
end
