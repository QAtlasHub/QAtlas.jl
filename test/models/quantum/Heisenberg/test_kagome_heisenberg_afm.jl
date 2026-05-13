# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: KagomeHeisenbergAFM — DMRG reference values.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "KagomeHeisenbergAFM — DMRG reference energy density" begin
    for J in (0.5, 1.0, 2.0)
        e = QAtlas.fetch(KagomeHeisenbergAFM(; J=J), Energy(:per_site), Infinite())
        @test e ≈ -0.4386 * J atol = 1e-14
    end
end

@testset "KagomeHeisenbergAFM — DMRG reference spin gap" begin
    for J in (0.5, 1.0, 2.0)
        Δ = QAtlas.fetch(KagomeHeisenbergAFM(; J=J), MassGap(), Infinite())
        @test Δ ≈ 0.13 * J atol = 1e-14
    end
end

@testset "KagomeHeisenbergAFM — DomainError on J < 0" begin
    @test_throws DomainError QAtlas.fetch(
        KagomeHeisenbergAFM(; J=-1.0), Energy(:per_site), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        KagomeHeisenbergAFM(; J=-0.5), MassGap(), Infinite()
    )
end
