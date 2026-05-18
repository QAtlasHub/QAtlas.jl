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

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "Heisenberg1D LuttingerParameter — verification cards" begin
    # K = 1/2 at the SU(2) isotropic point (Luther-Peschel 1975)
    verify(
        Heisenberg1D(),
        LuttingerParameter(),
        Infinite();
        route=:limiting_case,
        independent=0.5,
        agree_within=1e-12,
        refs=["Luther-Peschel 1975: K=1/2 at SU(2) isotropic point, J-independent"],
    )

    # Delegation invariant: Heisenberg1D === XXZ1D(Delta=1)
    verify(
        Heisenberg1D(),
        LuttingerParameter(),
        Infinite();
        route=:delegation_invariant,
        independent=QAtlas.fetch(XXZ1D(; J=1.0, Δ=1.0), LuttingerParameter(), Infinite()),
        agree_within=1e-14,
        refs=["Heisenberg1D delegates to XXZ1D(Delta=1): two code paths must agree"],
    )
end
