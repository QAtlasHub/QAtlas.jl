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

@testset "KagomeHeisenbergAFM — TopologicalEntanglementEntropy = log 2 (Phase 2, Z₂)" begin
    γ = QAtlas.fetch(KagomeHeisenbergAFM(), TopologicalEntanglementEntropy(), Infinite())
    @test γ ≈ log(2.0)
    # J-independence: γ is topological, not energy-scale-dependent
    @test γ == QAtlas.fetch(
        KagomeHeisenbergAFM(; J=3.7), TopologicalEntanglementEntropy(), Infinite()
    )
    # Z₂ sibling cross-check: ToricCode also has γ = log 2 (same Z₂ topological
    # order, total quantum dimension 𝒟 = 2 ⇒ γ = log 𝒟 = log 2).
    @test γ == QAtlas.fetch(ToricCode(), TopologicalEntanglementEntropy(), Infinite())
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "KagomeHeisenbergAFM — verification cards" begin
    # The 2D kagome AFM is a frustrated quantum spin liquid with no
    # closed form; all reference numbers are DMRG/iDMRG literature values.
    verify(
        KagomeHeisenbergAFM(; J=1.0),
        Energy(:per_site),
        Infinite();
        route=:literature_value,
        independent=-0.4386,
        agree_within=5e-3,
        refs=[
            "Yan-Huse-White 2011; Depenbrock-McCulloch-Schollwöck 2012 DMRG: e ≈ -0.4386 J"
        ],
    )

    verify(
        KagomeHeisenbergAFM(; J=1.0),
        MassGap(),
        Infinite();
        route=:literature_value,
        independent=0.13,
        agree_within=5e-2,
        refs=["Depenbrock et al. 2012 DMRG: spin gap ≈ 0.13 J"],
    )

    verify(
        KagomeHeisenbergAFM(; J=1.0),
        TopologicalEntanglementEntropy(),
        Infinite();
        route=:literature_value,
        independent=log(2),
        agree_within=1e-6,
        refs=["Z2 spin liquid: gamma = log 2 (Jiang-Wang-Balents 2012)"],
    )
end
# ── additional verification cards (#381 batch 4) ─────────────────────────
@testset "KagomeHeisenbergAFM — DMRG Energy + MassGap (#381 batch 4)" begin
    # Kagome Heisenberg AF DMRG: e₀ ≈ -0.4386 per spin and singlet-triplet
    # gap Δ ≈ 0.13 J (Yan-Huse-White 2011 Science 332, 1173;
    # Depenbrock-McCulloch-Schollwöck 2012 PRL 109, 067201).
    verify(
        KagomeHeisenbergAFM(),
        Energy(:per_site),
        Infinite();
        route=:literature_value,
        independent=-0.4386,
        agree_within=1e-3,
        refs=["Yan-Huse-White 2011 Science 332 1173: Kagome HAFM DMRG e₀ ≈ -0.4386 per spin"],
    )
    verify(
        KagomeHeisenbergAFM(),
        MassGap(),
        Infinite();
        route=:literature_value,
        independent=0.13,
        agree_within=1e-2,
        refs=["Yan-Huse-White 2011: Kagome HAFM singlet-triplet gap Δ ≈ 0.13 J (DMRG)"],
    )
end
# ── additional verification cards (#381 batch 3) ─────────────────────────
@testset "KagomeHeisenbergAFM — TEE Z2 spin liquid (#381 batch 3)" begin
    # KagomeHeisenberg AF is widely believed to realise a Z2 topological
    # spin liquid (Yan-Huse-White 2011 DMRG; Depenbrock-McCulloch-Schollwöck
    # 2012). Z2 topological order ⇒ total quantum dimension D = 2 ⇒
    # TEE γ = log D = log 2.
    verify(
        KagomeHeisenbergAFM(),
        TopologicalEntanglementEntropy(),
        Infinite();
        route=:literature_value,
        independent=log(2),
        agree_within=1e-12,
        refs=["Yan-Huse-White 2011; Depenbrock et al. 2012: Kagome HAFM Z2 spin liquid ⇒ TEE = log 2"],
    )
end

