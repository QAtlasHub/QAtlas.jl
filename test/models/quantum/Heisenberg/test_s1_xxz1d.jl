using QAtlas, Test

@testset "S1XXZ1D — Δ=1 delegate to S1Heisenberg1D (Phase 1)" begin
    Δ_gap = QAtlas.fetch(S1XXZ1D(; J=1.0, Δ=1.0), MassGap(), Infinite())
    @test Δ_gap > 0
    @test isapprox(Δ_gap, 0.41048; atol=1e-4)  # White 1992 DMRG Haldane gap
    # Delegation invariant
    @test Δ_gap ≈ QAtlas.fetch(S1Heisenberg1D(; J=1.0), MassGap(), Infinite())
    # Linear in J
    Δ3 = QAtlas.fetch(S1XXZ1D(; J=3.0, Δ=1.0), MassGap(), Infinite())
    @test Δ3 ≈ 3 * Δ_gap
end

@testset "S1XXZ1D — Δ ≠ 1 throws DomainError (Phase 2)" begin
    @test_throws DomainError QAtlas.fetch(S1XXZ1D(; Δ=0.5), MassGap(), Infinite())
    @test_throws DomainError QAtlas.fetch(S1XXZ1D(; Δ=2.0), MassGap(), Infinite())
    @test_throws DomainError QAtlas.fetch(S1XXZ1D(; Δ=-1.0), MassGap(), Infinite())
end

@testset "S1XXZ1D — rejects J ≤ 0" begin
    @test_throws DomainError S1XXZ1D(; J=0.0)
    @test_throws DomainError S1XXZ1D(; J=-1.0)
end
