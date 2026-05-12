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
