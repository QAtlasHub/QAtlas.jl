using QAtlas, Test

@testset "TFIM CentralCharge — Infinite" begin
    # At the Ising critical point h = J, c = 1/2 (Ising CFT).
    @test QAtlas.fetch(TFIM(; J=1.0, h=1.0), CentralCharge(), Infinite()) == 0.5
    @test QAtlas.fetch(TFIM(; J=2.0, h=2.0), CentralCharge(), Infinite()) == 0.5
    @test QAtlas.fetch(TFIM(; J=0.5, h=0.5), CentralCharge(), Infinite()) == 0.5

    # Off-critical (gapped) phases return 0.0 — not NaN. The value 0 is the
    # natural CFT-side answer for a gapped chain (no low-energy CFT).
    @test QAtlas.fetch(TFIM(; J=1.0, h=0.0), CentralCharge(), Infinite()) == 0.0
    @test QAtlas.fetch(TFIM(; J=1.0, h=0.5), CentralCharge(), Infinite()) == 0.0
    @test QAtlas.fetch(TFIM(; J=1.0, h=2.0), CentralCharge(), Infinite()) == 0.0
    @test QAtlas.fetch(TFIM(; J=1.0, h=10.0), CentralCharge(), Infinite()) == 0.0

    # Result is plain Float64 — no NaN propagation.
    c = QAtlas.fetch(TFIM(; J=1.0, h=0.3), CentralCharge(), Infinite())
    @test c isa Float64
    @test !isnan(c)

    # Tolerance |h/J - 1| ≤ 1e-6 still detects criticality.
    @test QAtlas.fetch(TFIM(; J=1.0, h=1.0 + 1e-9), CentralCharge(), Infinite()) == 0.5
    @test QAtlas.fetch(TFIM(; J=1.0, h=1.0 + 1e-3), CentralCharge(), Infinite()) == 0.0
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "TFIM CentralCharge — verification cards" begin
    # Ising universality: c = 1/2 at the critical point h = J
    verify(
        TFIM(; J=1.0, h=1.0),
        CentralCharge(),
        Infinite();
        route=:literature_value,
        independent=0.5,
        agree_within=1e-9,
        refs=["2D Ising CFT (Belavin-Polyakov-Zamolodchikov 1984): c = 1/2"],
    )

    # Gapped phases (h != J) have no critical theory: c = 0
    for (J, h) in ((1.0, 0.5), (1.0, 2.0), (1.0, 10.0))
        verify(
            TFIM(; J=J, h=h),
            CentralCharge(),
            Infinite();
            route=:second_closed_form,
            independent=0.0,
            agree_within=1e-9,
            refs=["Gapped phase (h != J): no conformal sector, c = 0"],
        )
    end
end
