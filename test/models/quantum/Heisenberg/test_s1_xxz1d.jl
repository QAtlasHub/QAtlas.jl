using QAtlas, Test

@testset "S1XXZ1D — Δ=1 delegate to S1Heisenberg1D (Phase 1)" begin
    Δ_gap = QAtlas.fetch(S1XXZ1D(; J=1.0, Δ=1.0), MassGap(), Infinite())
    @test Δ_gap > 0
    @test isapprox(Δ_gap, 0.41048; atol=1e-4)  # White-Huse 1993 DMRG Haldane gap (PRB 48, 3844)
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

@testset "S1XXZ1D — Δ=1 Energy{:per_site} delegate to S1Heisenberg1D (Phase 1)" begin
    for J in (0.5, 1.0, 2.0, 3.7)
        e_xxz = QAtlas.fetch(S1XXZ1D(; J=J, Δ=1.0), Energy(:per_site), Infinite())
        e_hei = QAtlas.fetch(S1Heisenberg1D(; J=J), Energy(:per_site), Infinite())
        @test e_xxz ≈ e_hei
        @test e_xxz < 0
    end
    # Literature value at J=1
    @test isapprox(
        QAtlas.fetch(S1XXZ1D(; J=1.0, Δ=1.0), Energy(:per_site), Infinite()),
        -1.40148403897;
        atol=1e-10,
    )
end

@testset "S1XXZ1D — Δ ≠ 1 Energy{:per_site} throws DomainError (Phase 2)" begin
    @test_throws DomainError QAtlas.fetch(S1XXZ1D(; Δ=0.5), Energy(:per_site), Infinite())
    @test_throws DomainError QAtlas.fetch(S1XXZ1D(; Δ=2.0), Energy(:per_site), Infinite())
    @test_throws DomainError QAtlas.fetch(S1XXZ1D(; Δ=-1.0), Energy(:per_site), Infinite())
    # isone(Δ) strictness regression
    @test_throws DomainError QAtlas.fetch(
        S1XXZ1D(; Δ=1.0 + 1e-13), Energy(:per_site), Infinite()
    )
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "S1XXZ1D — verification cards" begin
    # Delta=1 delegates to S1Heisenberg1D (Haldane chain, White-Huse 1993).
    verify(
        S1XXZ1D(; J=1.0, Δ=1.0),
        Energy(:per_site),
        Infinite();
        route=:delegation_invariant,
        independent=QAtlas.fetch(S1Heisenberg1D(; J=1.0), Energy(:per_site), Infinite()),
        agree_within=1e-12,
        refs=["S1XXZ1D(Delta=1) delegates to S1Heisenberg1D: code paths must agree"],
    )

    verify(
        S1XXZ1D(; J=1.0, Δ=1.0),
        MassGap(),
        Infinite();
        route=:delegation_invariant,
        independent=QAtlas.fetch(S1Heisenberg1D(; J=1.0), MassGap(), Infinite()),
        agree_within=1e-12,
        refs=["S1XXZ1D(Delta=1) Haldane gap delegates to S1Heisenberg1D"],
    )

    # J-scaling linear at the isotropic point
    verify(
        S1XXZ1D(; J=3.0, Δ=1.0),
        Energy(:per_site),
        Infinite();
        route=:delegation_invariant,
        independent=3.0 *
                    QAtlas.fetch(S1Heisenberg1D(; J=1.0), Energy(:per_site), Infinite()),
        agree_within=1e-10,
        refs=["Linear J scaling: e(3J) = 3 e(J) for spin-1 Heisenberg point"],
    )
end
# ── additional verification cards (#381 batch 5) ─────────────────────────
@testset "S1XXZ1D — Δ=1 isotropic Haldane reduction (#381 batch 5)" begin
    # At Δ=1 the spin-1 XXZ chain delegates to the SU(2)-symmetric S=1
    # Heisenberg AF (the Haldane chain). DMRG (White 1993, PRL 69 2863):
    #   e₀ ≈ -1.401484038971   per site
    #   Δ_Haldane ≈ 0.41048      gap (Wang-Qin-Hu 2012 refinement)
    # route=:delegation_invariant captures the actual test semantics
    # (code-path delegation at Δ=1) more accurately than :limiting_case.
    # J-scan covers default J=1 and one off-default point J=2 to verify
    # the delegation respects linear J scaling.
    for J in (1.0, 2.0)
        verify(
            S1XXZ1D(; J=J, Δ=1.0),
            Energy(:per_site),
            Infinite();
            route=:delegation_invariant,
            independent=-1.40148403897 * J,
            agree_within=1e-6,
            refs=[
                "White 1993 PRL 69 2863: S=1 Heisenberg DMRG e₀ ≈ -1.4014840 · J (Δ=1 delegation; J-linear scaling)",
            ],
        )
        verify(
            S1XXZ1D(; J=J, Δ=1.0),
            MassGap(),
            Infinite();
            route=:delegation_invariant,
            independent=0.41048 * J,
            agree_within=1e-3,
            refs=[
                "White 1993 / Wang-Qin-Hu 2012: S=1 Heisenberg Haldane gap ≈ 0.41048 · J (Δ=1 delegation; J-linear scaling)",
            ],
        )
    end
end
