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
    # ── Energy{:per_site} = −(J_e + J_m) ──────────────────────────
    @testset "ε₀ = −(J_e + J_m)" begin
        for (Je, Jm) in ((1.0, 1.0), (2.0, 1.0), (1.0, 3.0), (0.5, 2.5), (3.7, 0.0))
            m = ToricCode(; J_e=Je, J_m=Jm)
            ε = QAtlas.fetch(m, Energy{:per_site}(), Infinite())
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

    # ── #819: one quantity per return SCHEMA ─────────────────────────────
    #
    # These were one `AnyonStatistics` whose returned NamedTuple changed shape with
    # its `type` kwarg — `:em_braiding` gave (label, mutual_phase, anyons) while
    # `:e`/`:m`/`:ε` gave (label, statistics, self_phase, quantum_dim, fusion).  A
    # caller could not use the result without branching on the argument it had just
    # passed.  The label is still a kwarg, because WHICH anyon is an instance
    # selector like a momentum; what is fixed is that the fields no longer move.
    @testset "AnyonSelfStatistics: every label returns the same fields" begin
        m = ToricCode()
        fields = (:label, :statistics, :self_phase, :quantum_dim, :fusion)
        for l in (:vacuum, Symbol("1"), :e, :m, :ε, :epsilon, :eps)
            r = QAtlas.fetch(m, AnyonSelfStatistics(); label=l)
            # the point of the split: the SCHEMA is label-independent
            @test keys(r) === fields
            @test r.quantum_dim == 1.0                 # Abelian theory throughout
            @test r.fusion == ((r.label, r.label) => Symbol("1"))
        end

        e = QAtlas.fetch(m, AnyonSelfStatistics(); label=:e)
        @test e.label === :e
        @test e.statistics === :boson
        @test e.self_phase == 0.0

        mAnyon = QAtlas.fetch(m, AnyonSelfStatistics(); label=:m)
        @test mAnyon.label === :m
        @test mAnyon.statistics === :boson
        @test mAnyon.self_phase == 0.0

        ε = QAtlas.fetch(m, AnyonSelfStatistics(); label=:ε)
        @test ε.label === :ε
        @test ε.statistics === :fermion
        @test ε.self_phase ≈ π atol = 1e-14
        # the aliases are the SAME anyon, not merely a passing call
        @test QAtlas.fetch(m, AnyonSelfStatistics(); label=:epsilon) == ε
        @test QAtlas.fetch(m, AnyonSelfStatistics(); label=:eps) == ε

        vac = QAtlas.fetch(m, AnyonSelfStatistics(); label=:vacuum)
        @test vac.label === Symbol("1")
        @test vac.statistics === :boson
        @test QAtlas.fetch(m, AnyonSelfStatistics(); label=Symbol("1")) == vac

        @test_throws ArgumentError QAtlas.fetch(m, AnyonSelfStatistics(); label=:bogus)
    end

    @testset "AnyonMutualStatistics: e/m braid is π, everything else trivial" begin
        m = ToricCode()
        em = QAtlas.fetch(m, AnyonMutualStatistics())        # default pair
        @test keys(em) === (:anyons, :mutual_phase)
        @test em.mutual_phase ≈ π atol = 1e-14
        @test em.anyons === (:e, :m)
        # braiding is symmetric in the pair
        @test QAtlas.fetch(m, AnyonMutualStatistics(); anyons=(:m, :e)).mutual_phase ≈ π atol =
            1e-14
        # ...and trivial for every other pair, including an anyon with itself
        for pair in ((:e, :e), (:m, :m), (:ε, :ε), (:e, :vacuum), (:vacuum, :m))
            @test QAtlas.fetch(m, AnyonMutualStatistics(); anyons=pair).mutual_phase == 0.0
        end
        @test_throws ArgumentError QAtlas.fetch(
            m, AnyonMutualStatistics(); anyons=(:e, :bogus)
        )
    end

    @testset "the two quantities are not interchangeable" begin
        # The defect the split removes, pinned: no argument to one produces the
        # other'"'"'s fields, so a caller can never get a shape it did not ask for.
        m = ToricCode()
        self_fields = keys(QAtlas.fetch(m, AnyonSelfStatistics(); label=:ε))
        mutual_fields = keys(QAtlas.fetch(m, AnyonMutualStatistics()))
        @test isempty(intersect(self_fields, mutual_fields))
        # ε is a fermion BECAUSE of the e/m mutual phase — the physics that used to
        # justify keeping them together, now stated across the two quantities
        εphase = QAtlas.fetch(m, AnyonSelfStatistics(); label=:ε).self_phase
        ephase = QAtlas.fetch(m, AnyonSelfStatistics(); label=:e).self_phase
        mphase = QAtlas.fetch(m, AnyonSelfStatistics(); label=:m).self_phase
        braid = QAtlas.fetch(m, AnyonMutualStatistics(); anyons=(:e, :m)).mutual_phase
        @test εphase ≈ ephase + mphase + braid atol = 1e-14
    end

    # ── Asymmetry: J_e ≠ J_m gives different ε but same degeneracy ───────
    @testset "J_e ≠ J_m: different ε, same GSD/γ" begin
        ma = ToricCode(; J_e=1.0, J_m=2.0)
        mb = ToricCode(; J_e=2.0, J_m=1.0)

        εa = QAtlas.fetch(ma, Energy{:per_site}(), Infinite())
        εb = QAtlas.fetch(mb, Energy{:per_site}(), Infinite())
        @test εa == εb  # symmetric in (J_e, J_m) — both give -3.0

        mc = ToricCode(; J_e=1.0, J_m=3.0)
        εc = QAtlas.fetch(mc, Energy{:per_site}(), Infinite())
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
        Energy{:per_site}(),
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
            refs=[
                "Kitaev 2003 §4.1 (Ann. Phys. 303): ToricCode GSD = 4^g on a genus-g closed orientable surface (purely topological)",
            ],
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
