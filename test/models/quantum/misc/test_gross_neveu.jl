# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: GrossNeveu — UV c = N at g = 0.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "GrossNeveu — UV free-fermion c = N at g = 0" begin
    for N in 1:5
        @test QAtlas.fetch(GrossNeveu(; N=N, g=0.0), CentralCharge(), Infinite()) == N
    end
end

@testset "GrossNeveu — DomainError on N < 1" begin
    @test_throws DomainError GrossNeveu(; N=0, g=0.0)
end

@testset "GrossNeveu — DomainError on g ≠ 0 (Phase 2)" begin
    @test_throws DomainError QAtlas.fetch(
        GrossNeveu(; N=2, g=0.5), CentralCharge(), Infinite()
    )
end

@testset "GrossNeveu — dynamic mass m_F = Λ exp(-π/(N g²)) (Phase 2)" begin
    # Λ=1, N=1, g=1 → exp(-π) ≈ 0.04321391826377226
    @test QAtlas.fetch(GrossNeveu(; N=1, g=1.0), MassGap(), Infinite(); Λ=1.0) ≈ exp(-π)
    # Λ=1, N=2, g=1 → exp(-π/2)
    @test QAtlas.fetch(GrossNeveu(; N=2, g=1.0), MassGap(), Infinite(); Λ=1.0) ≈ exp(-π / 2)
    # Λ=2, N=4, g=0.5 → 2·exp(-π/(4·0.25)) = 2·exp(-π)
    @test QAtlas.fetch(GrossNeveu(; N=4, g=0.5), MassGap(), Infinite(); Λ=2.0) ≈ 2 * exp(-π)
    # Asymptotic-free UV limit: g large → m_F → Λ
    @test QAtlas.fetch(GrossNeveu(; N=1, g=10.0), MassGap(), Infinite(); Λ=1.0) ≈
        exp(-π / 100)
    # Linear in Λ
    @test QAtlas.fetch(GrossNeveu(; N=2, g=1.5), MassGap(), Infinite(); Λ=3.0) ≈
        3 * QAtlas.fetch(GrossNeveu(; N=2, g=1.5), MassGap(), Infinite(); Λ=1.0)
    # Weak-coupling limit: g → 0⁺ ⇒ m_F → 0⁺ (essential singularity at g=0; dimensional transmutation).
    @test QAtlas.fetch(GrossNeveu(; N=1, g=0.1), MassGap(), Infinite(); Λ=1.0) < 1e-100
    # g=0.05 ⇒ exp(-π/0.0025) ≈ exp(-1256) which underflows to 0.0 exactly.
    @test QAtlas.fetch(GrossNeveu(; N=1, g=0.05), MassGap(), Infinite(); Λ=1.0) == 0.0
end

@testset "GrossNeveu — MassGap rejects Λ, g, N out of domain (Phase 2)" begin
    m = GrossNeveu(; N=1, g=1.0)
    @test_throws DomainError QAtlas.fetch(m, MassGap(), Infinite(); Λ=0.0)
    @test_throws DomainError QAtlas.fetch(m, MassGap(), Infinite(); Λ=-1.5)
    @test_throws DomainError QAtlas.fetch(m, MassGap(), Infinite(); Λ=1.0, g=0.0)
    @test_throws DomainError QAtlas.fetch(m, MassGap(), Infinite(); Λ=1.0, g=-0.5)
    @test_throws DomainError QAtlas.fetch(m, MassGap(), Infinite(); Λ=1.0, N=0)
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "GrossNeveu — verification cards" begin
    for N in (1, 2, 5)
        verify(
            GrossNeveu(; N=N, g=0.0),
            CentralCharge(),
            Infinite();
            route=:second_closed_form,
            independent=N,
            agree_within=1e-10,
            refs=["UV free-fermion fixed point: c = N at g = 0"],
        )
    end
end
# ── additional verification cards (#381 batch) ─────────────────────────────
@testset "GrossNeveu — additional verification cards (#381 batch)" begin
    # Large-N dynamic mass m_F = Λ exp(-π / (N g²)); at N=1, g=1, Λ=1 → exp(-π).
    verify(
        GrossNeveu(; N=1, g=1.0),
        MassGap(),
        Infinite();
        route=:second_closed_form,
        independent=exp(-π),
        agree_within=1e-12,
        refs=["Gross-Neveu 1974; Andrei-Lowenstein 1979: m_F = Λ exp(-π/(N g²))"],
        fetch_kw=(; Λ=1.0),
    )
end

