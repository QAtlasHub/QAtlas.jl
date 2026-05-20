# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: SchwingerModel — massless Schwinger m_γ = e/√π.
#
# Verifies:
#   * m_γ = e/√π exactly at m = 0 for several e values
#   * Scaling: m_γ(2e) = 2 m_γ(e)
#   * DomainError on e ≤ 0
#   * DomainError on m ≠ 0 (massive case is Phase 2)
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "SchwingerModel — massless m_γ = e/√π" begin
    for e in (0.5, 1.0, 2.0, 3.7)
        mγ = QAtlas.fetch(SchwingerModel(; e=e, m=0.0), MassGap(), Infinite())
        @test mγ ≈ e / sqrt(π) atol = 1e-14
    end
end

@testset "SchwingerModel — linear e-scaling" begin
    m1 = QAtlas.fetch(SchwingerModel(; e=1.0), MassGap(), Infinite())
    m2 = QAtlas.fetch(SchwingerModel(; e=2.0), MassGap(), Infinite())
    @test m2 ≈ 2 * m1 atol = 1e-14
end

@testset "SchwingerModel — DomainError on e ≤ 0" begin
    @test_throws DomainError QAtlas.fetch(SchwingerModel(; e=0.0), MassGap(), Infinite())
    @test_throws DomainError QAtlas.fetch(SchwingerModel(; e=-1.0), MassGap(), Infinite())
end

@testset "SchwingerModel — DomainError on m ≠ 0 (Phase 2)" begin
    @test_throws DomainError QAtlas.fetch(
        SchwingerModel(; e=1.0, m=0.1), MassGap(), Infinite()
    )
end

@testset "SchwingerModel — ChiralCondensate (massless Schwinger anomaly, Phase 2)" begin
    cond = QAtlas.fetch(SchwingerModel(; e=1.0), ChiralCondensate(), Infinite())
    @test cond ≈ -exp(MathConstants.eulergamma) / (2 * π^(3/2)) atol = 1e-14
    @test cond < 0  # conventional sign
    # Linear in e
    cond2 = QAtlas.fetch(SchwingerModel(; e=2.0), ChiralCondensate(), Infinite())
    @test cond2 ≈ 2 * cond atol = 1e-14
    # Verified numerical value at e=1 (Julia MathConstants.eulergamma)
    @test cond ≈ -0.1599288349216857 atol = 1e-12
end

@testset "SchwingerModel — ChiralCondensate rejects e ≤ 0 (Phase 2)" begin
    @test_throws DomainError QAtlas.fetch(
        SchwingerModel(; e=1.0), ChiralCondensate(), Infinite(); e=0.0
    )
    @test_throws DomainError QAtlas.fetch(
        SchwingerModel(; e=1.0), ChiralCondensate(), Infinite(); e=-1.5
    )
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "SchwingerModel — verification cards" begin
    for e in (1.0, 2.0, 3.0)
        verify(
            SchwingerModel(; e=e, m=0.0),
            MassGap(),
            Infinite();
            route=:second_closed_form,
            independent=e / sqrt(pi),
            agree_within=1e-9,
            refs=["Massless Schwinger model: m_gamma = e / sqrt(pi) (exact)"],
        )
    end
end

# ── additional verification card (#381 batch) ─────────────────────────────
@testset "SchwingerModel — ChiralCondensate closed form (#381 batch)" begin
    # Coleman-Jackiw-Susskind 1975 (Ann. Phys. 93, 267) eq. 6.13 and
    # Lowenstein-Swieca 1971 (Ann. Phys. 68, 172): in the massless 1+1-D
    # Schwinger model the axial U(1) anomaly drives a finite vacuum
    # condensate ⟨ψ̄ψ⟩ = -exp(γ_E) · e / (2π^{3/2}). Linear in e.
    # (Schwinger 1962 derived the photon mass gap, NOT this condensate;
    # do not re-attribute the condensate result to Schwinger 1962.)
    γ_E = MathConstants.eulergamma
    # e-sweep matches the MassGap sibling card (#381 batch): {0.5, 1.0, 2.0, 3.7}.
    for e in (0.5, 1.0, 2.0, 3.7)
        verify(
            SchwingerModel(; e=e, m=0.0),
            ChiralCondensate(),
            Infinite();
            route=:second_closed_form,
            independent=-exp(γ_E) * e / (2 * π^(3/2)),
            agree_within=1e-12,
            refs=["Coleman-Jackiw-Susskind 1975 (Ann. Phys. 93, 267) eq. 6.13: ⟨ψ̄ψ⟩ = -exp(γ_E)·e/(2π^{3/2}); Lowenstein-Swieca 1971 (Ann. Phys. 68, 172)"],
        )
    end
end

