# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: ToricCode — Kitaev (2003) Z₂ surface code
#
# Verifies all five closed-form quantities exposed by the model and the
# `DomainError` boundary cases for negative J / negative genus.
#
# Run targeted:
#   julia --project=test test/standalone/test_toric_code.jl
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "ToricCode — Kitaev 2003 closed-form quantities" begin
    # ── GroundStateEnergyDensity = −(J_e + J_m) ──────────────────────────
    @testset "ε₀ = −(J_e + J_m)" begin
        for (Je, Jm) in ((1.0, 1.0), (2.0, 1.0), (1.0, 3.0), (0.5, 2.5), (3.7, 0.0))
            m = ToricCode(; J_e=Je, J_m=Jm)
            ε = QAtlas.fetch(m, GroundStateEnergyDensity(), Infinite())
            @test ε == -(Je + Jm)
            @test ε isa Float64
        end
    end

    # ── MassGap = 2 min(J_e, J_m) ─────────────────────────────────────────
    @testset "Δ = 2·min(J_e, J_m)" begin
        @test QAtlas.fetch(ToricCode(; J_e=1, J_m=1), MassGap(), Infinite()) == 2.0
        @test QAtlas.fetch(ToricCode(; J_e=2, J_m=1), MassGap(), Infinite()) == 2.0
        @test QAtlas.fetch(ToricCode(; J_e=1, J_m=3), MassGap(), Infinite()) == 2.0
        @test QAtlas.fetch(ToricCode(; J_e=0.7, J_m=2.5), MassGap(), Infinite()) ≈ 1.4
        @test QAtlas.fetch(ToricCode(; J_e=0.0, J_m=5.0), MassGap(), Infinite()) == 0.0
    end

    # ── GroundStateDegeneracy on closed surface of genus g = 4^g ─────────
    @testset "GSD = 4^genus" begin
        m = ToricCode()
        @test QAtlas.fetch(m, GroundStateDegeneracy(), PBC(0); genus=0) == 1
        @test QAtlas.fetch(m, GroundStateDegeneracy(), PBC(0); genus=1) == 4
        @test QAtlas.fetch(m, GroundStateDegeneracy(), PBC(0); genus=2) == 16
        @test QAtlas.fetch(m, GroundStateDegeneracy(), PBC(0); genus=3) == 64
        @test QAtlas.fetch(m, GroundStateDegeneracy(), PBC(0); genus=4) == 256

        # default genus = 1 (torus)
        @test QAtlas.fetch(m, GroundStateDegeneracy(), PBC(0)) == 4

        # GSD is independent of (J_e, J_m)
        m2 = ToricCode(; J_e=2.5, J_m=0.3)
        @test QAtlas.fetch(m2, GroundStateDegeneracy(), PBC(0); genus=2) == 16
    end

    # ── Topological entanglement entropy γ = log 2 ───────────────────────
    @testset "γ = log 2 (Kitaev–Preskill / Levin–Wen)" begin
        m = ToricCode()
        γ = QAtlas.fetch(m, TopologicalEntanglementEntropy(), Infinite())
        @test γ ≈ log(2.0) atol = 1e-14
        @test γ ≈ 0.6931471805599453 atol = 1e-14

        # Independent of J
        γ2 = QAtlas.fetch(
            ToricCode(; J_e=2.7, J_m=0.13), TopologicalEntanglementEntropy(), Infinite()
        )
        @test γ2 == γ
    end

    # ── AnyonStatistics — default :em returns mutual phase π ─────────────
    @testset "AnyonStatistics(:em) → mutual phase π" begin
        m = ToricCode()
        em = QAtlas.fetch(m, AnyonStatistics(); type=:em)
        @test em isa NamedTuple
        @test em.label === :em_braiding
        @test em.mutual_phase ≈ π atol = 1e-14
        @test em.anyons === (:e, :m)

        # default kwarg
        @test QAtlas.fetch(m, AnyonStatistics()) == em

        # :em_braiding alias matches :em
        @test QAtlas.fetch(m, AnyonStatistics(); type=:em_braiding) == em

        # Individual anyons
        e = QAtlas.fetch(m, AnyonStatistics(); type=:e)
        @test e.label === :e
        @test e.statistics === :boson
        @test e.self_phase == 0.0
        @test e.quantum_dim == 1.0
        @test e.fusion == ((:e, :e) => Symbol("1"))

        mAnyon = QAtlas.fetch(m, AnyonStatistics(); type=:m)
        @test mAnyon.label === :m
        @test mAnyon.statistics === :boson

        ε = QAtlas.fetch(m, AnyonStatistics(); type=:ε)
        @test ε.label === :ε
        @test ε.statistics === :fermion
        @test ε.self_phase ≈ π atol = 1e-14
        @test QAtlas.fetch(m, AnyonStatistics(); type=:epsilon) == ε

        vac = QAtlas.fetch(m, AnyonStatistics(); type=:vacuum)
        @test vac.label === Symbol("1")
        @test vac.statistics === :boson
        @test QAtlas.fetch(m, AnyonStatistics(); type=Symbol("1")) == vac

        # Unknown type
        @test_throws ArgumentError QAtlas.fetch(m, AnyonStatistics(); type=:bogus)
    end

    # ── Asymmetry: J_e ≠ J_m gives different ε but same degeneracy ───────
    @testset "J_e ≠ J_m: different ε, same GSD/γ" begin
        ma = ToricCode(; J_e=1.0, J_m=2.0)
        mb = ToricCode(; J_e=2.0, J_m=1.0)

        εa = QAtlas.fetch(ma, GroundStateEnergyDensity(), Infinite())
        εb = QAtlas.fetch(mb, GroundStateEnergyDensity(), Infinite())
        @test εa == εb  # symmetric in (J_e, J_m) — both give -3.0

        mc = ToricCode(; J_e=1.0, J_m=3.0)
        εc = QAtlas.fetch(mc, GroundStateEnergyDensity(), Infinite())
        @test εc == -4.0
        @test εc != εa

        # Degeneracy / γ unchanged
        @test QAtlas.fetch(ma, GroundStateDegeneracy(), PBC(0); genus=1) ==
            QAtlas.fetch(mc, GroundStateDegeneracy(), PBC(0); genus=1) ==
            4
        @test QAtlas.fetch(ma, TopologicalEntanglementEntropy(), Infinite()) ==
            QAtlas.fetch(mc, TopologicalEntanglementEntropy(), Infinite())

        # MassGap differs across asymmetric models
        Δa = QAtlas.fetch(ma, MassGap(), Infinite())  # 2 min(1,2) = 2
        Δc = QAtlas.fetch(mc, MassGap(), Infinite())  # 2 min(1,3) = 2
        @test Δa == Δc == 2.0
        Δd = QAtlas.fetch(ToricCode(; J_e=0.5, J_m=3.0), MassGap(), Infinite())
        @test Δd == 1.0
        @test Δd != Δa
    end

    # ── Domain errors ────────────────────────────────────────────────────
    @testset "DomainError on negative J or genus" begin
        @test_throws DomainError ToricCode(; J_e=-0.1, J_m=1.0)
        @test_throws DomainError ToricCode(; J_e=1.0, J_m=-2.0)
        @test_throws DomainError ToricCode(; J_e=-1.0, J_m=-1.0)

        m = ToricCode()
        @test_throws DomainError QAtlas.fetch(m, GroundStateDegeneracy(), PBC(0); genus=-1)
        @test_throws DomainError QAtlas.fetch(m, GroundStateDegeneracy(), PBC(0); genus=-3)
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "ToricCode — verification cards" begin
    verify(
        ToricCode(; J_e=1.0, J_m=1.0),
        GroundStateEnergyDensity(),
        Infinite();
        route=:second_closed_form,
        independent=-2.0,
        agree_within=1e-10,
        refs=["Kitaev 2003: toric-code e0 = -(J_e + J_m)"],
    )
    verify(
        ToricCode(; J_e=1.0, J_m=1.0),
        MassGap(),
        Infinite();
        route=:second_closed_form,
        independent=2.0,
        agree_within=1e-10,
        refs=["Toric-code gap = 2 min(J_e, J_m)"],
    )
    verify(
        ToricCode(; J_e=1.0, J_m=1.0),
        TopologicalEntanglementEntropy(),
        Infinite();
        route=:second_closed_form,
        independent=log(2),
        agree_within=1e-10,
        refs=["Z2 topological order: gamma = log 2 (Kitaev-Preskill)"],
    )
end

# ── additional verification card (#381 batch) ─────────────────────────────
@testset "ToricCode — GroundStateDegeneracy/PBC (#381 batch)" begin
    # Kitaev 2003 §4.1: GSD(g) = 4^g on a closed orientable surface of
    # genus g. Independent of J_e, J_m — purely topological
    # (dim H₁(Σ_g; Z₂) = 2g, so 2^{2g} = 4^g logical states; logical
    # operators are labelled by first homology classes). Sphere g=0 ⇒
    # unique GS; torus g=1 ⇒ canonical 4-fold; double torus g=2 ⇒ 16.
    # PBC(0) matches the call style used by the existing testsets above
    # (PBC() and PBC(0) dispatch identically for this hub).
    for genus in (0, 1, 2, 3)
        verify(
            ToricCode(; J_e=1.0, J_m=1.0),
            GroundStateDegeneracy(),
            PBC(0);
            route=:second_closed_form,
            independent=4^genus,
            agree_within=0,
            refs=["Kitaev 2003 §4.1 (Ann. Phys. 303): ToricCode GSD = 4^g on a genus-g closed orientable surface (purely topological)"],
            fetch_kw=(; genus=genus),
        )
    end
    # Same model with asymmetric couplings — degeneracy must be independent.
    for (J_e, J_m) in ((0.5, 2.0), (3.0, 1.0))
        verify(
            ToricCode(; J_e=J_e, J_m=J_m),
            GroundStateDegeneracy(),
            PBC(0);
            route=:second_closed_form,
            independent=4,
            agree_within=0,
            refs=["Kitaev 2003 §4.1: GSD is purely topological — independent of J_e, J_m"],
            fetch_kw=(; genus=1),
        )
    end
end

