using QAtlas, Test

@testset "Heisenberg1D — LuttingerParameter K = 1/2 (Phase 2, Luther-Peschel)" begin
    K = QAtlas.fetch(Heisenberg1D(), LuttingerParameter(), Infinite())
    # Strict ==: acos(1.0) == 0.0 exactly in IEEE, so π/(2π) == 0.5 exactly.
    @test K == 0.5
    # J-independence (strict ==)
    for J in (0.5, 1.0, 3.0)
        @test QAtlas.fetch(Heisenberg1D(), LuttingerParameter(), Infinite(); J=J) == 0.5
    end
    # Delegation invariant: bit-identical to XXZ1D at Δ=1
    @test K === QAtlas.fetch(QAtlas.XXZ1D(; Δ=1.0), LuttingerParameter(), Infinite())
end
