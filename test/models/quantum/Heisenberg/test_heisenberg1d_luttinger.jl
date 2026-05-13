using QAtlas, Test

@testset "Heisenberg1D — LuttingerParameter K = 1/2 (Phase 2, Luther-Peschel)" begin
    K = QAtlas.fetch(Heisenberg1D(), LuttingerParameter(), Infinite())
    @test K ≈ 0.5
    # J-independence
    for J in (0.5, 1.0, 3.0)
        @test QAtlas.fetch(Heisenberg1D(), LuttingerParameter(), Infinite(); J=J) ≈ 0.5
    end
    # Delegation invariant: matches XXZ1D at Δ=1
    @test K ≈ QAtlas.fetch(QAtlas.XXZ1D(; Δ=1.0), LuttingerParameter(), Infinite())
end
