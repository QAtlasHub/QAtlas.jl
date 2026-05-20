using QAtlas, Test

@testset "MixedFieldIsing1D — h_z = 0 delegate to TFIM (Phase 1)" begin
    # Paramagnetic side h_x > J: Δ = 2(h_x − J)
    @test QAtlas.fetch(
        MixedFieldIsing1D(; J=1.0, h_x=2.0, h_z=0.0), MassGap(), Infinite()
    ) ≈ 2.0
    @test QAtlas.fetch(
        MixedFieldIsing1D(; J=1.0, h_x=3.0, h_z=0.0), MassGap(), Infinite()
    ) ≈ 4.0
    # Ferromagnetic side h_x < J: Δ = 2(J − h_x)
    @test QAtlas.fetch(
        MixedFieldIsing1D(; J=1.0, h_x=0.5, h_z=0.0), MassGap(), Infinite()
    ) ≈ 1.0
    @test QAtlas.fetch(
        MixedFieldIsing1D(; J=2.0, h_x=1.0, h_z=0.0), MassGap(), Infinite()
    ) ≈ 2.0
    # Quantum critical point h_x = J: gap closes
    @test QAtlas.fetch(
        MixedFieldIsing1D(; J=1.0, h_x=1.0, h_z=0.0), MassGap(), Infinite()
    ) ≈ 0.0 atol = 1e-12
    # Delegation invariant: MixedFieldIsing1D at h_z = 0 matches TFIM directly.
    # TFIM stores (J, h) in struct fields, so the call signature is
    # `fetch(TFIM(; J, h=h_x), MassGap(), Infinite())`.
    Δ_mf = QAtlas.fetch(MixedFieldIsing1D(; J=1.5, h_x=0.7, h_z=0.0), MassGap(), Infinite())
    Δ_tfim = QAtlas.fetch(TFIM(; J=1.5, h=0.7), MassGap(), Infinite())
    @test Δ_mf ≈ Δ_tfim
    # Default constructor sits on the Phase-1 point (J = h_x = 1, h_z = 0; QCP)
    @test QAtlas.fetch(MixedFieldIsing1D(), MassGap(), Infinite()) ≈ 0.0 atol = 1e-12
end

@testset "MixedFieldIsing1D — h_z ≠ 0 throws DomainError (non-integrable)" begin
    @test_throws DomainError QAtlas.fetch(
        MixedFieldIsing1D(; J=1.0, h_x=1.0, h_z=0.1), MassGap(), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        MixedFieldIsing1D(; J=1.0, h_x=1.0, h_z=-0.5), MassGap(), Infinite()
    )
end

@testset "MixedFieldIsing1D — rejects J ≤ 0 (Phase 1)" begin
    @test_throws DomainError MixedFieldIsing1D(; J=0.0)
    @test_throws DomainError MixedFieldIsing1D(; J=-1.0)
end

@testset "MixedFieldIsing1D — strict h_z=0 boundary (no isapprox)" begin
    # A subnormal h_z must NOT silently delegate to TFIM. The Phase-1
    # delegate boundary is strict (iszero), matching the KitaevHeisenberg
    # convention. Even h_z = 1e-13 selects the non-integrable regime.
    @test_throws DomainError QAtlas.fetch(
        MixedFieldIsing1D(; J=1.0, h_x=1.0, h_z=1e-13), MassGap(), Infinite()
    )
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "MixedFieldIsing1D — verification cards" begin
    # h_z = 0 reduces to TFIM: gap Δ = 2|h_x - J| (Pfeuty).
    for (J, hx) in ((1.0, 2.0), (1.0, 3.0), (2.0, 1.0))
        verify(
            MixedFieldIsing1D(; J=J, h_x=hx, h_z=0.0),
            MassGap(),
            Infinite();
            route=:delegation_invariant,
            independent=2 * abs(hx - J),
            agree_within=1e-9,
            refs=["h_z=0 delegates to TFIM: Δ = 2|h_x - J| (Pfeuty 1970)"],
        )
    end
end
# ── additional verification cards (#381 batch 6) ─────────────────────────
@testset "MixedFieldIsing1D — TFIM critical line MassGap = 0 (#381 batch 6)" begin
    # At h_z=0 the mixed-field Ising H = -J Σ σ^z σ^z - h_x Σ σ^x reduces
    # to TFIM, which is gapless on the critical line h_x = J. The struct
    # uses fields (J, h_x, h_z); MFI at h_x=J, h_z=0 gives Δ = 0.
    # Single J point: 2|h_x - J| = 0.0 identically at h_x = J for any J;
    # a multi-J sweep would collapse to identical zero residuals on the
    # same closed-form code path and adds no discriminating power.
    let J = 1.0
        verify(
            MixedFieldIsing1D(; J=J, h_x=J, h_z=0.0),
            MassGap(),
            Infinite();
            route=:limiting_case,
            independent=0.0,
            agree_within=1e-9,
            refs=["TFIM critical line h_x=J (Pfeuty 1970) ⇒ MFI at h_z=0, h_x=J gives MassGap = 0"],
        )
    end
end

