# ─────────────────────────────────────────────────────────────────────────────
# Test: ExtendedHubbard1D — Phase 1, V = 0 delegate to Hubbard1D.
#
# Cases (issue #294):
#   1. V = 0 delegate: ChargeGap matches Hubbard1D(μ = U/2) for several
#      (t, U) points along the Lieb–Wu curve.
#   2. V ≠ 0 raises DomainError (Phase 2 deferred).
#   3. t ≤ 0 in the constructor raises DomainError.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "ExtendedHubbard1D — V=0 delegate to Hubbard1D (Phase 1)" begin
    for (t, U) in [(1.0, 4.0), (1.0, 8.0), (0.5, 2.0)]
        Δ_eh = QAtlas.fetch(ExtendedHubbard1D(; t=t, U=U, V=0.0), ChargeGap(), Infinite())
        Δ_h = QAtlas.fetch(Hubbard1D(; t=t, U=U, μ=U / 2), ChargeGap(), Infinite())
        @test Δ_eh ≈ Δ_h
    end
end

@testset "ExtendedHubbard1D — V ≠ 0 throws DomainError (Phase 2)" begin
    @test_throws DomainError QAtlas.fetch(
        ExtendedHubbard1D(; V=0.5), ChargeGap(), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        ExtendedHubbard1D(; V=-1.0), ChargeGap(), Infinite()
    )
end

@testset "ExtendedHubbard1D — rejects t ≤ 0" begin
    @test_throws DomainError ExtendedHubbard1D(; t=0.0)
    @test_throws DomainError ExtendedHubbard1D(; t=-1.0)
end
# ── additional verification cards (#381 batch 6) ─────────────────────────
@testset "ExtendedHubbard1D — ChargeGap V=0 delegation to Hubbard1D (#381 batch 6)" begin
    # At V=0 the t-U-V Extended Hubbard chain reduces to the Lieb-Wu
    # Hubbard1D at half-filling. The src delegates ChargeGap directly to
    # Hubbard1D(t, U, U/2). Cross-check via independent delegate fetch.
    for (t, U) in ((1.0, 1.0), (1.0, 4.0), (0.5, 4.0))
        ref = QAtlas.fetch(Hubbard1D(; t=t, U=U, μ=U/2), ChargeGap(), Infinite())
        verify(
            ExtendedHubbard1D(; t=t, U=U, V=0.0),
            ChargeGap(),
            Infinite();
            route=:delegation_invariant,
            independent=ref,
            agree_within=1e-12,
            refs=["ExtendedHubbard1D at V=0 ≡ Lieb-Wu Hubbard1D at half-filling ⇒ ChargeGap delegation equality"],
        )
    end
end

