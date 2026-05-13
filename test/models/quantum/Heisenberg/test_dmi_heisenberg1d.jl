# test/models/quantum/Heisenberg/test_dmi_heisenberg1d.jl
#
# Phase-1 closed-form tests for DMIHeisenberg1D (issue #298).  Only the
# D = 0 analytically-known point is exercised here:
#
#   - D = 0   → Bethe-Hulthén  (delegate Heisenberg1D)
#
# Generic D ≠ 0 (twisted-XXZ / spiral, Affleck-Oshikawa 1999) is
# deferred to Phase 2 and must raise DomainError.

using QAtlas, Test

@testset "DMIHeisenberg1D — D = 0 delegate to Heisenberg1D (Phase 1)" begin
    e0 = QAtlas.fetch(DMIHeisenberg1D(; J=1.0, D=0.0), Energy{:per_site}(), Infinite())
    @test e0 ≈ 0.25 - log(2)
    # Default IS D = 0
    @test QAtlas.fetch(DMIHeisenberg1D(), Energy{:per_site}(), Infinite()) ≈ 0.25 - log(2)
    # Linear in J
    @test QAtlas.fetch(DMIHeisenberg1D(; J=3.0, D=0.0), Energy{:per_site}(), Infinite()) ≈
        3 * (0.25 - log(2))
    # Delegation matches Heisenberg1D directly (legacy GroundStateEnergyDensity API)
    @test e0 ≈ QAtlas.fetch(Heisenberg1D(), GroundStateEnergyDensity(), Infinite(); J=1.0)
end

@testset "DMIHeisenberg1D — D ≠ 0 throws DomainError (Phase 2 deferral)" begin
    @test_throws DomainError QAtlas.fetch(
        DMIHeisenberg1D(; J=1.0, D=0.1), Energy{:per_site}(), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        DMIHeisenberg1D(; J=1.0, D=-0.5), Energy{:per_site}(), Infinite()
    )
end

@testset "DMIHeisenberg1D — rejects J ≤ 0 (Phase 1)" begin
    @test_throws DomainError DMIHeisenberg1D(; J=0.0)
    @test_throws DomainError DMIHeisenberg1D(; J=-1.5)
end
