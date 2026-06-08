# ─────────────────────────────────────────────────────────────────────────────
# Test: Su-Schrieffer-Heeger (1979) 1D dimerised tight-binding chain.
#
# Values are pinned by the verify() cards (closed forms with independent
# witnesses). This file also keeps the structural / error / identity guards
# verify() cannot express, plus the load-bearing INDEPENDENT cross-check: the
# Infinite per-site energy (Gauss-Kronrod over the band) and the bulk gap must
# agree with a direct dense-ED diagonalisation of the OBC chain — two
# genuinely different computations.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test
using QAtlas:
    SSH,
    ExactSpectrum,
    Energy,
    MassGap,
    CorrelationLength,
    TopologicalInvariant,
    EdgeModeEnergy,
    OBC,
    Infinite,
    fetch

@testset "SSH — structural / error / identity guards" begin
    @testset "TopologicalInvariant: gapless-line (|v|=|w|) error" begin
        @test_throws ErrorException fetch(
            SSH(; v=1.0, w=1.0), TopologicalInvariant(), Infinite()
        )
        @test_throws ErrorException fetch(
            SSH(; v=0.7, w=0.7), TopologicalInvariant(), Infinite()
        )
        # v = −w is ALSO gapless (q(0) = 0), not only v = w
        @test_throws ErrorException fetch(
            SSH(; v=-1.0, w=1.0), TopologicalInvariant(), Infinite()
        )
    end

    @testset "OBC ExactSpectrum shape (length N, sorted, non-negative)" begin
        m = SSH(; v=0.6, w=1.0)
        for N in (8, 16, 24)
            spec = fetch(m, ExactSpectrum(), OBC(N))
            @test length(spec) == N
            @test issorted(spec)
            @test all(spec .>= -1e-10)
        end
    end

    @testset "EdgeModeEnergy/OBC == MassGap/OBC (definitional identity)" begin
        m = SSH(; v=0.4, w=1.0)
        N = 32
        @test fetch(m, EdgeModeEnergy(), OBC(N)) == fetch(m, MassGap(), OBC(N))
        # MassGap@OBC is the smallest ExactSpectrum entry.
        @test fetch(m, MassGap(), OBC(N)) == fetch(m, ExactSpectrum(), OBC(N))[1]
    end

    @testset "topological vs trivial OBC edge mode (w>v vs v>w)" begin
        # Topological (w>v): edge mode ≪ bulk gap and shrinks with N.
        topo = SSH(; v=0.4, w=1.0)
        e20 = fetch(topo, EdgeModeEnergy(), OBC(20))
        e40 = fetch(topo, EdgeModeEnergy(), OBC(40))
        @test e40 < 1e-3                         # ≪ single-particle gap |v-w| = 0.6
        @test e40 < e20                          # exponential decay in N
        # Trivial (v>w): lowest OBC level is of order the bulk gap (no edge mode).
        triv = SSH(; v=1.0, w=0.4)
        @test fetch(triv, EdgeModeEnergy(), OBC(40)) >
            0.5 * fetch(triv, MassGap(), Infinite())
    end

    @testset "CorrelationLength gapless-line: ξ = Inf (structural)" begin
        @test fetch(SSH(; v=1.0, w=1.0), CorrelationLength(), Infinite()) == Inf
    end

    @testset "Energy(:per_site) Infinite — finite/negative + :natural delegation" begin
        for (v, w) in ((1.0, 0.4), (0.4, 1.0), (1.0, 1.0))
            ε = fetch(SSH(; v=v, w=w), Energy(:per_site), Infinite())
            @test isfinite(ε) && ε < 0
        end
        @test fetch(SSH(), Energy(), Infinite()) ==
            fetch(SSH(), Energy(:per_site), Infinite())
    end
end

# ── INDEPENDENT cross-check: closed forms vs direct dense-ED ───────────────────
@testset "SSH — Infinite closed forms agree with OBC dense-ED" begin
    @testset "per-site energy: Gauss-Kronrod integral == OBC ED average" begin
        for (v, w) in ((1.0, 0.4), (0.4, 1.0), (0.7, 1.3), (-0.5, 0.7), (0.6, -1.0))
            m = SSH(; v=v, w=w)
            ε_inf = fetch(m, Energy(:per_site), Infinite())
            N = 200
            # half-filled OBC ground energy per site = (filled negative band) / 2N
            #   = -Σ(non-negative single-particle energies) / 2N — a dense-ED
            #     computation independent of the Gauss-Kronrod band integral.
            ε_obc = -sum(fetch(m, ExactSpectrum(), OBC(N))) / (2N)
            @test isapprox(ε_inf, ε_obc; atol=1e-2)
        end
    end

    @testset "single-particle gap ||v|-|w|| == OBC gap in trivial phase (no edge mode)" begin
        # trivial |v|>|w| (no edge mode), incl. OPPOSITE-SIGN hoppings — this is the
        # independent ED witness that catches the ||v|-|w|| vs |v-w| convention.
        for (v, w) in ((1.0, 0.4), (-1.0, 0.4), (1.0, -0.4))
            m = SSH(; v=v, w=w)
            @test isapprox(
                fetch(m, MassGap(), OBC(60)), fetch(m, MassGap(), Infinite()); atol=1e-2
            )
        end
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────────
@testset "SSH — verification cards" begin
    # MassGap Infinite = single-particle gap ||v|-|w|| = min_k|q(k)| (closed form).
    # Includes OPPOSITE-SIGN hoppings (vw<0), where the minimum sits at k=0 and the
    # naive |v-w| would be wrong (e.g. (-0.5,0.7): ||v|-|w||=0.2, not |v-w|=1.2).
    for (v, w, gap) in (
        (1.0, 0.4, 0.6),
        (0.4, 1.0, 0.6),
        (0.0, 1.0, 1.0),
        (1.0, 0.0, 1.0),
        (-0.5, 0.7, 0.2),
        (0.6, -1.0, 0.4),
        (-1.0, -0.4, 0.6),
    )
        verify(
            SSH(; v=v, w=w),
            MassGap(),
            Infinite();
            route=:second_closed_form,
            independent=gap,
            agree_within=1e-12,
            refs=["SSH 1979: single-particle gap = min_k|q(k)| = ||v|−|w|| (band gap 2×)"],
        )
    end

    # TopologicalInvariant winding W: 1 (|w|>|v|) / 0 (|w|<|v|).
    # Fetch integrates Im(q'/q); the independent witness is the |w|≷|v| threshold.
    for (v, w) in (
        (0.4, 1.0),
        (0.0, 1.0),
        (1.0, 0.4),
        (1.0, 0.0),
        (-0.5, 1.2),
        (1.3, -0.5),
        (0.3, -1.5),
    )
        W_expected = abs(w) > abs(v) ? 1.0 : 0.0   # |W| ∈ {0,1}; topological iff |w|>|v|
        verify(
            SSH(; v=v, w=w),
            TopologicalInvariant(),
            Infinite();
            route=:second_closed_form,
            independent=W_expected,
            agree_within=1e-9,
            refs=["SSH 1979; Asbóth-Oroszlány-Pályi 2016: W = 1 (|w|>|v|) / 0 (|w|<|v|)"],
        )
    end

    # Energy(:per_site) Infinite at the fully dimerised sweet spots: flat band,
    # ε₀ = −max(|v|,|w|)/2.  Independent of the integral.
    for t in (0.5, 1.0, 2.0)
        verify(
            SSH(; v=0.0, w=t),
            Energy(:per_site),
            Infinite();
            route=:second_closed_form,
            independent=(-t / 2),
            agree_within=1e-9,
            refs=["SSH 1979 sweet spot v=0: flat band |q|=|w| ⇒ ε₀ = −|w|/2"],
        )
        verify(
            SSH(; v=t, w=0.0),
            Energy(:per_site),
            Infinite();
            route=:second_closed_form,
            independent=(-t / 2),
            agree_within=1e-9,
            refs=["SSH 1979 sweet spot w=0: flat band |q|=|v| ⇒ ε₀ = −|v|/2"],
        )
    end

    # CorrelationLength Infinite = 1/||v|-|w||.
    verify(
        SSH(; v=1.0, w=0.4),
        CorrelationLength(),
        Infinite();
        route=:second_closed_form,
        independent=(1 / 0.6),
        agree_within=1e-9,
        refs=["SSH 1979: ξ = 1/Δ_gap = 1/||v|−|w||; v=1,w=0.4 ⇒ ξ = 1/0.6"],
    )
    verify(
        SSH(; v=-0.5, w=0.7),                          # opposite-sign: gap min at k=0
        CorrelationLength(),
        Infinite();
        route=:second_closed_form,
        independent=(1 / 0.2),
        agree_within=1e-9,
        refs=["SSH 1979 opposite-sign: ξ = 1/||v|−|w|| = 1/0.2 = 5"],
    )

    # EdgeModeEnergy OBC at the topological sweet spot (v=0): the end sites
    # decouple, so the boundary zero modes are EXACT for any N.
    for N in (6, 8, 16, 32)
        verify(
            SSH(; v=0.0, w=1.0),
            EdgeModeEnergy(),
            OBC(N);
            route=:second_closed_form,
            independent=0.0,
            agree_within=1e-12,
            refs=[
                "SSH 1979 topological sweet spot v=0: exact end zero modes (E_edge = 0 ∀N)"
            ],
        )
    end
end
