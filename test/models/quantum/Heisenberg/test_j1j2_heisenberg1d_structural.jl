# test/models/quantum/Heisenberg/test_j1j2_heisenberg1d_structural.jl
#
# Phase-1 closed-form tests + delegation verify + MG-point closed-form
# verify cards for J1J2Heisenberg1D (issue #297). Trivially fast — no ED.
#
# Split out of test_j1j2_heisenberg1d.jl (9.9 min on s02). The ED sweep
# lives in the sibling test_j1j2_heisenberg1d_verify_ed.jl.

using QAtlas, Test

@testset "J1J2Heisenberg1D — j = 0 (pure Heisenberg, Bethe-Hulthén) (Phase 1)" begin
    e0 = QAtlas.fetch(J1J2Heisenberg1D(; J1=1.0, J2=0.0), Energy{:per_site}(), Infinite())
    @test e0 ≈ 0.25 - log(2)
    e0_3 = QAtlas.fetch(J1J2Heisenberg1D(; J1=3.0, J2=0.0), Energy{:per_site}(), Infinite())
    @test e0_3 ≈ 3 * (0.25 - log(2))
    @test e0 ≈ QAtlas.fetch(Heisenberg1D(), GroundStateEnergyDensity(), Infinite(); J=1.0)
end

@testset "J1J2Heisenberg1D — j = 1/2 (Majumdar-Ghosh) (Phase 1)" begin
    e0 = QAtlas.fetch(J1J2Heisenberg1D(; J1=1.0, J2=0.5), Energy{:per_site}(), Infinite())
    @test e0 ≈ -3 / 8
    @test QAtlas.fetch(J1J2Heisenberg1D(), Energy{:per_site}(), Infinite()) ≈ -3 / 8
    e0_2 = QAtlas.fetch(J1J2Heisenberg1D(; J1=2.0, J2=1.0), Energy{:per_site}(), Infinite())
    @test e0_2 ≈ -3 / 4
    @test e0 ≈ QAtlas.fetch(MajumdarGhosh(; J=1.0), GroundStateEnergyDensity(), Infinite())
end

@testset "J1J2Heisenberg1D — generic j throws DomainError (Phase 1)" begin
    for (J1, J2) in [(1.0, 0.25), (1.0, 0.75), (2.0, 0.3), (1.5, 1.0)]
        @test_throws DomainError QAtlas.fetch(
            J1J2Heisenberg1D(; J1=J1, J2=J2), Energy{:per_site}(), Infinite()
        )
    end
end

@testset "J1J2Heisenberg1D — rejects J1 ≤ 0 or J2 < 0 (Phase 1)" begin
    @test_throws DomainError J1J2Heisenberg1D(; J1=0.0, J2=0.5)
    @test_throws DomainError J1J2Heisenberg1D(; J1=-1.0, J2=0.5)
    @test_throws DomainError J1J2Heisenberg1D(; J1=1.0, J2=-0.1)
end

# ── Verification cards (delegation + closed-form, no ED) ──────────────────
@testset "J1J2Heisenberg1D — delegation verify cards" begin
    # j = 0 (pure Heisenberg): delegates to Heisenberg1D, e0 = J1(1/4 - log 2)
    verify(
        J1J2Heisenberg1D(; J1=1.0, J2=0.0),
        Energy(:per_site),
        Infinite();
        route=:delegation_invariant,
        independent=QAtlas.fetch(Heisenberg1D(), GroundStateEnergyDensity(), Infinite()),
        agree_within=1e-12,
        refs=["J1J2 at J2=0 delegates to Heisenberg1D (Hulthen 1938)"],
    )

    # j = 1/2 (Majumdar-Ghosh): delegates to MajumdarGhosh, e0 = -3 J1 / 8
    verify(
        J1J2Heisenberg1D(; J1=1.0, J2=0.5),
        Energy(:per_site),
        Infinite();
        route=:delegation_invariant,
        independent=QAtlas.fetch(
            MajumdarGhosh(; J=1.0), GroundStateEnergyDensity(), Infinite()
        ),
        agree_within=1e-12,
        refs=["J1J2 at J2=J1/2 delegates to MajumdarGhosh (exact dimer, -3J/8)"],
    )
end

# ── additional verification cards (#381 batch 3) ─────────────────────────
@testset "J1J2Heisenberg1D — Majumdar-Ghosh point Energy (#381 batch 3)" begin
    for J1 in (0.5, 1.0, 2.0)
        verify(
            J1J2Heisenberg1D(; J1=J1, J2=J1 / 2),
            Energy(:per_site),
            Infinite();
            route=:second_closed_form,
            independent=-3 * J1 / 8,
            agree_within=1e-14,
            refs=[
                "Majumdar-Ghosh 1969: at J2 = J1/2 GS is exact dimer-product state ⇒ e₀ = -3J1/8",
            ],
        )
    end
end
