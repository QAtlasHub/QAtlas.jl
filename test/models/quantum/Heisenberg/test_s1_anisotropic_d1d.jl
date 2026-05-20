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

@testset "S1AnisotropicD1D — tiny D (1e-13) is non-zero (iszero strictness)" begin
    # Regression: boundary check must be exact iszero(D), not isapprox(...; atol=1e-12).
    # Tiny but non-zero D crosses Gaussian/Ising transitions and has no closed form.
    @test_throws DomainError QAtlas.fetch(
        S1AnisotropicD1D(; J=1.0, D=1e-13), MassGap(), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        S1AnisotropicD1D(; J=1.0, D=-1e-13), MassGap(), Infinite()
    )
end

@testset "S1AnisotropicD1D — constructor guards" begin
    # J ≤ 0 rejected at construction time.
    @test_throws DomainError S1AnisotropicD1D(; J=0.0, D=0.0)
    @test_throws DomainError S1AnisotropicD1D(; J=-1.0, D=0.0)
    # D may be any real, including negative.
    @test S1AnisotropicD1D(; J=1.0, D=-1.5) isa S1AnisotropicD1D
end

# ─────────────────────────────────────────────────────────────────────────────
# Energy{:per_site} — D = 0 delegate to S1Heisenberg1D (White-Huse 1993).
# ─────────────────────────────────────────────────────────────────────────────

@testset "S1AnisotropicD1D — D=0 Energy{:per_site} delegate" begin
    for J in (0.5, 1.0, 2.0, 3.7)
        m = S1AnisotropicD1D(; J=J, D=0.0)
        e = QAtlas.fetch(m, Energy{:per_site}(), Infinite())
        e_delegate = QAtlas.fetch(S1Heisenberg1D(; J=J), Energy{:per_site}(), Infinite())
        @test e ≈ e_delegate
        @test isapprox(e, -1.40148403897 * J; atol=1e-10)
    end
end

@testset "S1AnisotropicD1D — Energy{:per_site} D≠0 throws DomainError" begin
    @test_throws DomainError QAtlas.fetch(
        S1AnisotropicD1D(; J=1.0, D=0.1), Energy{:per_site}(), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        S1AnisotropicD1D(; J=1.0, D=-0.5), Energy{:per_site}(), Infinite()
    )
end

@testset "S1AnisotropicD1D — Energy{:per_site} tiny D (1e-13) strict iszero" begin
    # Regression: boundary check must be exact iszero(D), not isapprox.
    @test_throws DomainError QAtlas.fetch(
        S1AnisotropicD1D(; J=1.0, D=1e-13), Energy{:per_site}(), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        S1AnisotropicD1D(; J=1.0, D=-1e-13), Energy{:per_site}(), Infinite()
    )
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "S1AnisotropicD1D — verification cards" begin
    # D=0 delegates to S1Heisenberg1D (Haldane chain, White-Huse 1993).
    verify(
        S1AnisotropicD1D(; J=1.0, D=0.0),
        Energy(:per_site),
        Infinite();
        route=:delegation_invariant,
        independent=QAtlas.fetch(S1Heisenberg1D(; J=1.0), Energy(:per_site), Infinite()),
        agree_within=1e-12,
        refs=["S1AnisotropicD1D(D=0) delegates to S1Heisenberg1D: code paths must agree"],
    )

    verify(
        S1AnisotropicD1D(; J=1.0, D=0.0),
        MassGap(),
        Infinite();
        route=:delegation_invariant,
        independent=QAtlas.fetch(S1Heisenberg1D(; J=1.0), MassGap(), Infinite()),
        agree_within=1e-12,
        refs=["S1AnisotropicD1D(D=0) Haldane gap delegates to S1Heisenberg1D"],
    )

    # J-scaling linear at D=0
    verify(
        S1AnisotropicD1D(; J=2.0, D=0.0),
        Energy(:per_site),
        Infinite();
        route=:delegation_invariant,
        independent=2.0 *
                    QAtlas.fetch(S1Heisenberg1D(; J=1.0), Energy(:per_site), Infinite()),
        agree_within=1e-10,
        refs=["Linear J scaling: e(2J) = 2 e(J) for spin-1 Heisenberg"],
    )
end
# ── additional verification cards (#381 batch 5) ─────────────────────────
@testset "S1AnisotropicD1D — D=0 isotropic Haldane reduction (#381 batch 5)" begin
    # At D=0 the S=1 single-ion anisotropy chain delegates to the
    # SU(2)-symmetric S=1 Heisenberg AF (the Haldane chain).
    # DMRG (White 1993): e₀ ≈ -1.40148 J, Haldane gap ≈ 0.41048 J.
    # route=:delegation_invariant captures the actual test semantics
    # (code-path delegation at D=0) more accurately than :limiting_case.
    # J-scan covers default J=1 and one off-default point J=2 to verify
    # the delegation respects linear J scaling.
    for J in (1.0, 2.0)
        verify(
            S1AnisotropicD1D(; J=J, D=0.0),
            Energy(:per_site),
            Infinite();
            route=:delegation_invariant,
            independent=-1.40148403897 * J,
            agree_within=1e-6,
            refs=["White 1993 PRL 69 2863: S=1 Heisenberg DMRG e₀ ≈ -1.4014840 · J (D=0 delegation; J-linear scaling)"],
        )
        verify(
            S1AnisotropicD1D(; J=J, D=0.0),
            MassGap(),
            Infinite();
            route=:delegation_invariant,
            independent=0.41048 * J,
            agree_within=1e-3,
            refs=["White 1993 / Wang-Qin-Hu 2012: S=1 Heisenberg Haldane gap ≈ 0.41048 · J (D=0 delegation; J-linear scaling)"],
        )
    end
end

