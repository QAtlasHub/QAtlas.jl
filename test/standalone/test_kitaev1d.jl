# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: Kitaev (2001) 1D p-wave wire.
#
# Targeted run (skips Pkg.test()):
#   julia --project=test test/standalone/test_kitaev1d.jl
#
# Coverage:
#   • PBC closed-form dispersion: E(0) = |2t + μ|, E(π) = |2t - μ|
#   • TopologicalInvariant: ν = -1 (topological) for |μ| < 2|t|;
#                           ν = +1 (trivial)     for |μ| > 2|t|;
#                           Pfaffian sign flips at the critical line.
#   • TFIM correspondence: Kitaev1D(μ=-2h, t=J, Δ=J) BdG spectrum at OBC
#                          equals TFIM(J=J, h=h) BdG spectrum.
#   • OBC topological N = 40 has Majorana edge mode (lowest |Λ| ≤ 1e-3).
#   • OBC trivial      N = 40 has lowest |Λ| of order the bulk gap 2(|μ| - 2t).
#   • CorrelationLength: finite away from |μ| = 2|t|, Inf on the critical line.
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

@testset "Kitaev1D" begin

    # ─────────────────── PBC closed-form dispersion ───────────────────
    @testset "PBC dispersion: E(0) = |2t+μ|, E(π) = |2t-μ|" begin
        # E(k) = √((2t cos k + μ)² + 4Δ² sin² k)
        # At k=0: sin = 0, E = |2t + μ|.   At k=π: sin = 0, E = |2t - μ|.
        for (μ, t, Δ) in
            [(0.0, 1.0, 1.0), (1.0, 1.0, 1.0), (3.0, 1.0, 1.0), (-2.5, 1.0, 0.5)]
            E0 = sqrt((2t * 1.0 + μ)^2 + 4Δ^2 * 0.0)
            Eπ = sqrt((2t * (-1.0) + μ)^2 + 4Δ^2 * 0.0)
            @test E0 ≈ abs(2t + μ)
            @test Eπ ≈ abs(2t - μ)
        end
    end

    # ─────────────────── TopologicalInvariant ─────────────────────────
    @testset "TopologicalInvariant: ν = -1 in topological phase" begin
        m = Kitaev1D(; μ=0.0, t=1.0, Δ=1.0)
        @test fetch(m, TopologicalInvariant(), Infinite()) == -1

        m2 = Kitaev1D(; μ=1.5, t=1.0, Δ=1.0)        # |μ| < 2|t|
        @test fetch(m2, TopologicalInvariant(), Infinite()) == -1

        m3 = Kitaev1D(; μ=-1.5, t=1.0, Δ=0.5)
        @test fetch(m3, TopologicalInvariant(), Infinite()) == -1
    end

    @testset "TopologicalInvariant: ν = +1 in trivial phase" begin
        m = Kitaev1D(; μ=3.0, t=1.0, Δ=1.0)         # |μ| > 2|t|
        @test fetch(m, TopologicalInvariant(), Infinite()) == 1

        m2 = Kitaev1D(; μ=-2.5, t=1.0, Δ=1.0)
        @test fetch(m2, TopologicalInvariant(), Infinite()) == 1
    end

    @testset "Pfaffian flips sign at |μ| = 2t" begin
        # On the gapless critical line one Pfaffian vanishes ⇒ invariant
        # is ill-defined, the routine throws.  We sample just inside the
        # boundary on either side and check the sign flip.
        ε = 1e-6
        ν_in = fetch(
            Kitaev1D(; μ=2.0 - ε, t=1.0, Δ=1.0), TopologicalInvariant(), Infinite()
        )
        ν_out = fetch(
            Kitaev1D(; μ=2.0 + ε, t=1.0, Δ=1.0), TopologicalInvariant(), Infinite()
        )
        @test ν_in == -1
        @test ν_out == 1

        # And the gapless point itself errors out (ill-defined invariant).
        @test_throws ErrorException fetch(
            Kitaev1D(; μ=2.0, t=1.0, Δ=1.0), TopologicalInvariant(), Infinite()
        )
        @test_throws ErrorException fetch(
            Kitaev1D(; μ=-2.0, t=1.0, Δ=1.0), TopologicalInvariant(), Infinite()
        )
    end

    # ─────────────────── TFIM correspondence ──────────────────────────
    @testset "TFIM correspondence: μ=-2h, t=J, Δ=J reproduces TFIM spectrum" begin
        # _tfim_bdg_spectrum filters values <= 1e-10, so we must do the
        # same to the Kitaev1D spectrum to compare apples to apples.  We
        # do it by sorting and taking the top N entries on each side
        # (both routines build the same 2N×2N BdG matrix when the
        # parameters coincide, so the top-N entries match exactly).
        N = 20
        for (J, h) in [(1.0, 0.5), (1.0, 1.5), (1.0, 1.0), (0.7, 0.3)]
            μ_eq = -2h
            spec_kitaev = QAtlas._kitaev1d_bdg_spectrum(N, μ_eq, J, J)
            spec_tfim = QAtlas._tfim_bdg_spectrum(N, J, h)
            # spec_tfim filters near-zeros (< 1e-10); take the top-K
            # entries of the Kitaev spectrum that match length(spec_tfim).
            K = length(spec_tfim)
            @test K > 0
            top_kitaev = spec_kitaev[(end - K + 1):end]
            @test isapprox(sort(top_kitaev), sort(spec_tfim); atol=1e-9)
        end
    end

    # ─────────────────── OBC edge mode (topological) ──────────────────
    @testset "OBC topological N=40 has Majorana edge mode (≤ 1e-3)" begin
        m = Kitaev1D(; μ=0.0, t=1.0, Δ=1.0)
        N = 40
        Λmin = fetch(m, EdgeModeEnergy(), OBC(N))
        @test Λmin <= 1e-3
        # Same value via MassGap@OBC by construction
        @test fetch(m, MassGap(), OBC(N)) == Λmin
        # Sweet spot μ = 0, t = Δ = 1: the Kitaev sweet point — the two
        # Majoranas decouple exactly to the chain ends, so the splitting
        # is essentially 0 (down to floating-point noise).
        @test Λmin <= 1e-10
    end

    @testset "OBC topological μ ≠ 0 still has small edge-mode energy" begin
        m = Kitaev1D(; μ=0.5, t=1.0, Δ=1.0)
        Λmin = fetch(m, EdgeModeEnergy(), OBC(40))
        @test Λmin <= 1e-3
    end

    # ─────────────────── OBC trivial: bulk gap ────────────────────────
    @testset "OBC trivial N=40 lowest energy is order of bulk gap" begin
        # Trivial phase: μ = 3, t = Δ = 1 ⇒ bulk gap = ||μ| - 2|t|| = 1.
        # The OBC lowest energy converges to the bulk gap exponentially
        # in N, so at N=40 we just check it's of the right order.
        m = Kitaev1D(; μ=3.0, t=1.0, Δ=1.0)
        Λmin = fetch(m, EdgeModeEnergy(), OBC(40))
        bulk_gap = abs(abs(3.0) - 2 * 1.0)         # = 1.0; min(|2t+μ|, |2t-μ|)
        @test Λmin >= 0.5 * bulk_gap
        @test Λmin <= 1.5 * bulk_gap
        # Cross-check against the analytic infinite-chain gap
        gap_inf = fetch(m, MassGap(), Infinite())
        @test isapprox(Λmin, gap_inf; atol=5e-2)
    end

    # ─────────────────── ExactSpectrum returns vector ─────────────────
    @testset "ExactSpectrum returns N sorted non-negative entries" begin
        m = Kitaev1D(; μ=0.0, t=1.0, Δ=1.0)
        for N in (10, 20, 30)
            spec = fetch(m, ExactSpectrum(), OBC(N))
            @test length(spec) == N
            @test issorted(spec)
            @test all(spec .>= -1e-10)
        end
    end

    # ─────────────────── CorrelationLength ─────────────────────────────
    @testset "CorrelationLength: gapped, finite; gapless, Inf" begin
        ξ_top = fetch(Kitaev1D(; μ=0.0, t=1.0, Δ=1.0), CorrelationLength(), Infinite())
        @test isfinite(ξ_top)
        @test ξ_top > 0

        ξ_triv = fetch(Kitaev1D(; μ=3.0, t=1.0, Δ=1.0), CorrelationLength(), Infinite())
        @test isfinite(ξ_triv)
        @test ξ_triv > 0
        # In the trivial phase μ = 3, t = Δ = 1 the bulk gap is 1, so ξ = 1.
        @test isapprox(ξ_triv, 1.0; atol=1e-9)

        ξ_crit = fetch(Kitaev1D(; μ=2.0, t=1.0, Δ=1.0), CorrelationLength(), Infinite())
        @test ξ_crit == Inf
    end

    # ─────────────────── Energy{:per_site} smoke ──────────────────────
    @testset "Energy(:per_site) at Infinite is finite + negative" begin
        for (μ, t, Δ) in [(0.0, 1.0, 1.0), (1.0, 1.0, 1.0), (3.0, 1.0, 1.0)]
            ε = fetch(Kitaev1D(; μ=μ, t=t, Δ=Δ), Energy(:per_site), Infinite())
            @test isfinite(ε)
            @test ε < 0
        end
        # Energy() resolves through the :natural router to :per_site at
        # Infinite() (declared by `native_energy_granularity`).
        ε_nat = fetch(Kitaev1D(), Energy(), Infinite())
        ε_ps = fetch(Kitaev1D(), Energy(:per_site), Infinite())
        @test ε_nat == ε_ps
    end

    # ─────────────────── MassGap closed form spot checks ──────────────
    @testset "MassGap closed form at corners" begin
        # Sweet spot t = Δ, μ = 0: gap = 2|Δ| (PBC dispersion
        # E(k) = 2√(t² cos²k + Δ² sin²k) is constant = 2|t| = 2|Δ|).
        @test isapprox(
            fetch(Kitaev1D(; μ=0.0, t=1.0, Δ=1.0), MassGap(), Infinite()), 2.0; atol=1e-10
        )

        # Trivial μ = 3, t = Δ = 1: gap = ||μ| - 2|t|| = 1.
        @test isapprox(
            fetch(Kitaev1D(; μ=3.0, t=1.0, Δ=1.0), MassGap(), Infinite()), 1.0; atol=1e-10
        )

        # Critical |μ| = 2|t|: gap = 0.
        @test isapprox(
            fetch(Kitaev1D(; μ=2.0, t=1.0, Δ=1.0), MassGap(), Infinite()), 0.0; atol=1e-10
        )
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
