using QAtlas, Test

# Phase 1 of S1AnisotropicD1D (#302): a single-ion-anisotropy spin-1
# chain that, at D = 0, delegates to the existing S1Heisenberg1D entry
# for the Haldane gap (White-Huse 1993 DMRG, Δ ≈ 0.41048 J).  Non-zero
# D throws DomainError — Phase 2 (Chen-Roncaglia 2008; Tzeng-Yang-Hsu 2017).

@testset "S1AnisotropicD1D — D=0 delegate (Phase 1)" begin
    m = S1AnisotropicD1D(; J=1.0, D=0.0)
    Δ = QAtlas.fetch(m, MassGap(), Infinite())
    @test Δ > 0
    # White-Huse 1993 DMRG value as encoded by S1Heisenberg1D.
    @test isapprox(Δ, 0.41048; atol=1e-6)
    # Delegation invariant: identical to direct S1Heisenberg1D fetch.
    Δ_delegate = QAtlas.fetch(S1Heisenberg1D(; J=1.0), MassGap(), Infinite())
    @test Δ ≈ Δ_delegate
    # J scaling propagates through the delegate.
    Δ2 = QAtlas.fetch(S1AnisotropicD1D(; J=2.0, D=0.0), MassGap(), Infinite())
    @test isapprox(Δ2, 2.0 * Δ; atol=1e-10)
end

@testset "S1AnisotropicD1D — D≠0 throws DomainError" begin
    @test_throws DomainError QAtlas.fetch(
        S1AnisotropicD1D(; J=1.0, D=0.1), MassGap(), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        S1AnisotropicD1D(; J=1.0, D=-0.5), MassGap(), Infinite()
    )
end

@testset "S1AnisotropicD1D — constructor guards" begin
    # J ≤ 0 rejected at construction time.
    @test_throws DomainError S1AnisotropicD1D(; J=0.0, D=0.0)
    @test_throws DomainError S1AnisotropicD1D(; J=-1.0, D=0.0)
    # D may be any real, including negative.
    @test S1AnisotropicD1D(; J=1.0, D=-1.5) isa S1AnisotropicD1D
end
