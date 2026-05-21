# ─────────────────────────────────────────────────────────────────────────────
# Test: Kitaev (2001) 1D p-wave wire.
#
# Values are verified by the verify() cards below. This file retains
# only the structural / error / identity / relational guards that
# verify() architecturally cannot express (DomainError on the gapless
# line, TFIM↔Kitaev BdG-spectrum cross-model identity, OBC ExactSpectrum
# shape, EdgeModeEnergy == MassGap@OBC identity, gapless-metal guard,
# trivial-phase bulk-gap relational bounds).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test
using QAtlas:
    Kitaev1D,
    TFIM,
    ExactSpectrum,
    Energy,
    MassGap,
    CorrelationLength,
    TopologicalInvariant,
    EdgeModeEnergy,
    OBC,
    Infinite,
    fetch

@testset "Kitaev1D — structural / error / identity guards" begin
    @testset "TopologicalInvariant: gapless-line DomainError" begin
        # On |μ|=2t one Pfaffian vanishes ⇒ invariant ill-defined.
        @test_throws ErrorException fetch(
            Kitaev1D(; μ=2.0, t=1.0, Δ=1.0), TopologicalInvariant(), Infinite()
        )
        @test_throws ErrorException fetch(
            Kitaev1D(; μ=-2.0, t=1.0, Δ=1.0), TopologicalInvariant(), Infinite()
        )
    end

    @testset "TFIM correspondence: Kitaev1D BdG spectrum == TFIM BdG (cross-model)" begin
        # Cross-model identity: Kitaev1D(μ=-2h, t=J, Δ=J) reproduces the
        # TFIM(J,h) BdG spectrum. Not a single-value verify card.
        N = 20
        for (J, h) in [(1.0, 0.5), (1.0, 1.5), (1.0, 1.0), (0.7, 0.3)]
            μ_eq = -2h
            spec_kitaev = QAtlas._kitaev1d_bdg_spectrum(N, μ_eq, J, J)
            spec_tfim = QAtlas._tfim_bdg_spectrum(N, J, h)
            K = length(spec_tfim)
            @test K > 0
            top_kitaev = spec_kitaev[(end - K + 1):end]
            @test isapprox(sort(top_kitaev), sort(spec_tfim); atol=1e-9)
        end
    end

    @testset "OBC ExactSpectrum shape (length N, sorted, non-negative)" begin
        m = Kitaev1D(; μ=0.0, t=1.0, Δ=1.0)
        for N in (10, 20, 30)
            spec = fetch(m, ExactSpectrum(), OBC(N))
            @test length(spec) == N
            @test issorted(spec)
            @test all(spec .>= -1e-10)
        end
    end

    @testset "EdgeModeEnergy/OBC == MassGap/OBC (definitional identity)" begin
        # Registry: same value; named for the Majorana boundary-mode
        # interpretation. Identity between two fetches — not a card.
        m = Kitaev1D(; μ=0.0, t=1.0, Δ=1.0)
        N = 40
        @test fetch(m, EdgeModeEnergy(), OBC(N)) == fetch(m, MassGap(), OBC(N))
    end

    @testset "OBC trivial-phase EdgeModeEnergy is order of bulk gap" begin
        # Relational: at trivial μ=3, t=Δ=1 the bulk gap is |μ|-2t = 1;
        # the OBC lowest energy is of that order and approaches it
        # exponentially in N (no Majorana edge mode in trivial phase).
        m = Kitaev1D(; μ=3.0, t=1.0, Δ=1.0)
        Λmin = fetch(m, EdgeModeEnergy(), OBC(40))
        bulk_gap = abs(abs(3.0) - 2 * 1.0)              # = 1.0
        @test 0.5 * bulk_gap <= Λmin <= 1.5 * bulk_gap
        # Cross-check against the analytic infinite-chain gap.
        @test isapprox(Λmin, fetch(m, MassGap(), Infinite()); atol=5e-2)
    end

    @testset "CorrelationLength gapless-line: ξ = Inf (structural)" begin
        # ξ on |μ|=2t diverges (no verify card can encode Inf cleanly).
        @test fetch(Kitaev1D(; μ=2.0, t=1.0, Δ=1.0), CorrelationLength(), Infinite()) == Inf
    end

    @testset "Energy(:per_site) Infinite — finite/negative + :natural delegation" begin
        # Relational/sanity: ε < 0 and finite; :natural router resolves
        # to :per_site at Infinite() per native_energy_granularity.
        for (μ, t, Δ) in [(0.0, 1.0, 1.0), (1.0, 1.0, 1.0), (3.0, 1.0, 1.0)]
            ε = fetch(Kitaev1D(; μ=μ, t=t, Δ=Δ), Energy(:per_site), Infinite())
            @test isfinite(ε) && ε < 0
        end
        @test fetch(Kitaev1D(), Energy(), Infinite()) ==
            fetch(Kitaev1D(), Energy(:per_site), Infinite())
    end
end

@testset "Kitaev1D gapless-metal guard (OBC)" begin
    m_metal = Kitaev1D(; μ=1.0, t=1.0, Δ=0.0)
    @test_throws ErrorException QAtlas.fetch(m_metal, MassGap(), OBC(20))
    @test_throws ErrorException QAtlas.fetch(m_metal, EdgeModeEnergy(), OBC(20))

    m_atomic = Kitaev1D(; μ=3.0, t=1.0, Δ=0.0)
    Δgap = QAtlas.fetch(m_atomic, MassGap(), OBC(20))
    @test isfinite(Δgap)
    @test Δgap > 0
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "Kitaev1D — verification cards" begin
    # MassGap Infinite at characteristic points (sweet/trivial/critical).
    verify(
        Kitaev1D(; μ=0.0, t=1.0, Δ=1.0),
        MassGap(),
        Infinite();
        route=:second_closed_form,
        independent=2.0,
        agree_within=1e-9,
        refs=["Kitaev chain sweet spot: gap = 2|Δ|"],
    )
    verify(
        Kitaev1D(; μ=3.0, t=1.0, Δ=1.0),
        MassGap(),
        Infinite();
        route=:second_closed_form,
        independent=1.0,
        agree_within=1e-9,
        refs=["Kitaev chain trivial phase: gap = ||μ| - 2|t||"],
    )
    verify(
        Kitaev1D(; μ=2.0, t=1.0, Δ=1.0),
        MassGap(),
        Infinite();
        route=:second_closed_form,
        independent=0.0,
        agree_within=1e-9,
        refs=["Kitaev chain topological transition at |μ|=2|t|: gap = 0"],
    )

    # TopologicalInvariant Infinite: closed form ν = sgn(μ²−4t²).
    # |μ|<2|t| topological ν=-1; |μ|>2|t| trivial ν=+1.
    for (μ, t, Δ) in (
        (0.0, 1.0, 1.0),
        (1.5, 1.0, 1.0),
        (-1.5, 1.0, 0.5),
        (3.0, 1.0, 1.0),
        (-2.5, 1.0, 1.0),
    )
        ν_closed = sign(μ^2 - 4 * t^2)              # -1 topological, +1 trivial
        verify(
            Kitaev1D(; μ=μ, t=t, Δ=Δ),
            TopologicalInvariant(),
            Infinite();
            route=:second_closed_form,
            independent=ν_closed,
            agree_within=1e-12,
            refs=[
                "Kitaev 2001; Asboth-Oroszlany-Palyi 2016: ν = sgn(μ²-4t²) " *
                "(-1 topological, +1 trivial)",
            ],
        )
    end

    # CorrelationLength Infinite in the trivial phase μ=3, t=Δ=1: ξ = 1/Δ_gap = 1.
    verify(
        Kitaev1D(; μ=3.0, t=1.0, Δ=1.0),
        CorrelationLength(),
        Infinite();
        route=:second_closed_form,
        independent=1.0,
        agree_within=1e-9,
        refs=["Kitaev chain trivial μ=3, t=Δ=1: bulk gap = 1 ⇒ ξ = 1/Δ_gap = 1"],
    )

    # EdgeModeEnergy OBC at the sweet spot (μ=0, t=Δ=1): the two
    # Majoranas decouple exactly to the chain ends ⇒ splitting ≈ 0.
    verify(
        Kitaev1D(; μ=0.0, t=1.0, Δ=1.0),
        EdgeModeEnergy(),
        OBC(40);
        route=:second_closed_form,
        independent=0.0,
        agree_within=1e-9,
        refs=["Kitaev sweet spot μ=0, t=Δ=1: exact Majorana boundary, splitting ~ 0"],
    )

    # MassGap OBC at the sweet spot: same value as EdgeModeEnergy (defn).
    verify(
        Kitaev1D(; μ=0.0, t=1.0, Δ=1.0),
        MassGap(),
        OBC(40);
        route=:second_closed_form,
        independent=0.0,
        agree_within=1e-9,
        refs=["Kitaev sweet spot μ=0, t=Δ=1: OBC gap ≈ 0 (Majorana edge mode)"],
    )
end

# ── additional verification cards (#381 batch) ─────────────────────────────
@testset "Kitaev1D — additional closed-form cards (#381 batch)" begin
    # TopologicalInvariant/Infinite (Kitaev 2001 Pfaffian Z2): ν = sgn(μ² − 4t²).
    # |μ| < 2|t| ⇒ topological (ν = −1); |μ| > 2|t| ⇒ trivial (ν = +1).
    for (μ, t, Δ, ν_expected) in (
        (0.0, 1.0, 1.0, -1), (0.5, 1.0, 1.0, -1), (3.0, 1.0, 1.0, +1), (-3.0, 1.0, 0.5, +1)
    )
        verify(
            Kitaev1D(; μ=μ, t=t, Δ=Δ),
            TopologicalInvariant(),
            Infinite();
            route=:second_closed_form,
            independent=ν_expected,
            agree_within=0,
            refs=["Kitaev 2001 Pfaffian Z2: ν = sgn(μ² − 4t²) on the gapped phases"],
        )
    end

    # Energy/Infinite sweet-spot closed form: μ = 0, t = Δ ⇒ E(k) = 2t,
    # so ε₀ = −(1/2π) ∫₀^π E(k) dk = −(1/2π)·2t·π = −t. Independent of integral.
    for t in (0.5, 1.0, 2.0)
        verify(
            Kitaev1D(; μ=0.0, t=t, Δ=t),
            Energy(:per_site),
            Infinite();
            route=:second_closed_form,
            independent=(-t),
            agree_within=1e-9,
            refs=[
                "Kitaev 2001 sweet spot μ=0, t=Δ: dispersion is flat E(k)=2t, so ε₀ = −t"
            ],
        )
    end

    # CorrelationLength/Infinite = 1 / bulk gap (gapped phases).
    # Sweet spot (μ=0, t=Δ=1): gap = 2 ⇒ ξ = 1/2.
    verify(
        Kitaev1D(; μ=0.0, t=1.0, Δ=1.0),
        CorrelationLength(),
        Infinite();
        route=:second_closed_form,
        independent=0.5,
        agree_within=1e-12,
        refs=[
            "Kitaev 2001 sweet spot (μ=0, t=Δ): bulk gap = 2|Δ|, general ξ = 1/Δ_gap ⇒ ξ = 1/(2|Δ|)",
        ],
    )
    # Trivial phase (μ=3, t=1, Δ=1): gap = |μ|−2|t| = 1 ⇒ ξ = 1.
    verify(
        Kitaev1D(; μ=3.0, t=1.0, Δ=1.0),
        CorrelationLength(),
        Infinite();
        route=:second_closed_form,
        independent=1.0,
        agree_within=1e-12,
        refs=["Kitaev 2001 trivial phase: bulk gap = |μ|−2|t| ⇒ ξ = 1/(|μ|−2|t|)"],
    )

    # EdgeModeEnergy/OBC at the exact sweet spot (μ=0, t=Δ): the two
    # end-localised Majorana modes do NOT hybridise — the lowest BdG
    # eigenvalue is exactly 0 for any N (Kitaev 2001, original example).
    for N in (6, 8, 16, 32)
        verify(
            Kitaev1D(; μ=0.0, t=1.0, Δ=1.0),
            EdgeModeEnergy(),
            OBC(N);
            route=:second_closed_form,
            independent=0.0,
            agree_within=1e-10,
            refs=[
                "Kitaev 2001 sweet spot OBC: Majorana zero modes are exact (E_edge = 0 for any N)",
            ],
        )
    end
end
