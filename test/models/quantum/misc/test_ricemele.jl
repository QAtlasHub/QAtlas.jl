# ─────────────────────────────────────────────────────────────────────────────
# Test: Rice-Mele (1982) dimerised chain with a staggered on-site potential.
#
# The oracle is SSH: at Δ = 0 this model IS SSH, and SSH is an independent
# implementation here that shares no code with it. The rest is closed forms with
# a witness, plus the guards a copy of SSH would fail.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test
using Random: seed!
using QAtlas:
    RiceMele, SSH, ExactSpectrum, Energy, MassGap, CorrelationLength, OBC, Infinite, fetch

@testset "RiceMele — Δ = 0 is SSH, on an independent implementation" begin
    # Not one fixture: sweep the dimerisation, both signs, and the critical line.
    for (v, w) in ((0.6, 1.0), (1.0, 0.6), (1.0, 1.0), (0.0, 1.0), (1.0, 0.0), (-0.7, 1.3))
        rm = RiceMele(; v, w, Δ=0.0)
        ssh = SSH(; v, w)
        @test fetch(rm, Energy{:per_site}(), Infinite()) ≈
            fetch(ssh, Energy{:per_site}(), Infinite()) atol = 1.0e-12
        @test fetch(rm, MassGap(), Infinite()) ≈ fetch(ssh, MassGap(), Infinite()) atol =
            1.0e-14
        for N in (4, 9, 20)
            full = fetch(rm, ExactSpectrum(), OBC(N))
            @test length(full) == 2N
            # Chiral symmetry is back at Δ = 0, and this is the unambiguous statement of it.
            @test full ≈ -reverse(full) atol = 1.0e-10
            @test full[(N + 1):(2N)] ≈ fetch(ssh, ExactSpectrum(), OBC(N)) atol = 1.0e-10
            # And the two OBC gap conventions coincide exactly here, which is the
            # whole reason the Δ ≠ 0 one is written with a half in it.
            @test fetch(rm, MassGap(), OBC(N)) ≈ fetch(ssh, MassGap(), OBC(N)) atol =
                1.0e-10
        end
    end
end

@testset "RiceMele — the zero-mode degeneracy, where the non-negative half is N+1" begin
    # v = 0 gives a two-fold exact zero, so the non-negative energies number N+1, not N.
    # Returning all 2N leaves nothing to choose between them.
    N = 4
    full = fetch(RiceMele(; v=0.0, w=1.0, Δ=0.0), ExactSpectrum(), OBC(N))
    @test length(full) == 2N
    @test full ≈ [-1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0, 1.0] atol = 1.0e-12
    @test count(x -> abs(x) < 1.0e-12, full) == 2               # the two edge modes
    @test count(>=(-1.0e-12), full) == N + 1                    # N+1, not N
    @test length(fetch(SSH(; v=0.0, w=1.0), ExactSpectrum(), OBC(N))) == N
end

@testset "RiceMele — what Δ ≠ 0 changes, which a copy of SSH would not" begin
    m = RiceMele(; v=0.75, w=0.25, Δ=0.1)

    @testset "the OBC spectrum stays ±-symmetric, and why" begin
        # Δ breaks the BULK chiral symmetry, not this one: S = ΓR anticommutes with H
        # for a complete-cell chain. Tested because the first draft claimed the opposite.
        for N in (5, 12)
            ε = fetch(m, ExactSpectrum(), OBC(N))
            @test length(ε) == 2N
            @test issorted(ε)
            @test ε ≈ -reverse(ε) atol = 1.0e-12
            @test count(<(0), ε) == N
        end
        # The control: one site short and the antisymmetry is gone.
        odd = QAtlas._rice_mele_obc_spectrum_sites(9, 0.75, 0.25, 0.1)
        @test !isapprox(odd, -reverse(odd); atol=1.0e-3)
        @test maximum(abs.(odd + reverse(odd))) > 0.1     # measured 2Δ = 0.2 here
    end

    @testset "the gap never closes" begin
        # |v| = |w| is SSH's Dirac point. Δ gaps it out, and the value is exactly |Δ|.
        for Δ in (0.05, 0.3, -0.2)
            @test fetch(RiceMele(; v=1.0, w=1.0, Δ), MassGap(), Infinite()) ≈ abs(Δ) atol =
                1.0e-14
            @test isfinite(
                fetch(RiceMele(; v=1.0, w=1.0, Δ), CorrelationLength(), Infinite())
            )
        end
        @test fetch(SSH(; v=1.0, w=1.0), MassGap(), Infinite()) == 0.0   # the contrast
    end

    @testset "no TopologicalInvariant row exists" begin
        # Not quantised once Δ breaks the chiral symmetry, so the model must refuse.
        @test_throws ErrorException fetch(m, QAtlas.TopologicalInvariant(), Infinite())
    end
end

@testset "RiceMele — the ± symmetry over a random sweep, and where it stops" begin
    # The claim is that a complete-cell chain is ±-symmetric for EVERY (v,w,Δ), and that
    # the site count — not Δ — is what carries it. A single fixture cannot say either.
    seed!(846)
    worst_even, worst_odd = 0.0, 0.0
    for _ in 1:300
        v, w, Δ = 2rand() - 1, 2rand() - 1, 2rand() - 1
        for N in (3, 4, 7, 12)
            ε = QAtlas._rice_mele_obc_spectrum(N, v, w, Δ)
            worst_even = max(worst_even, maximum(abs.(ε + reverse(ε))))
        end
        for n in (5, 9, 15)
            o = QAtlas._rice_mele_obc_spectrum_sites(n, v, w, Δ)
            worst_odd = max(worst_odd, maximum(abs.(o + reverse(o))))
        end
    end
    @test worst_even < 1.0e-12          # measured 2.2e-15
    @test worst_odd > 0.1               # measured 2.0

    # And the (hx, hy, hz) mapping, over the same kind of sweep rather than three triples.
    seed!(1982)
    worst_map = 0.0
    for _ in 1:200
        hx, hy, hz = 2rand() - 1, 2rand() - 1, 2rand() - 1
        m = RiceMele(; v=(hx + hy) / 2, w=(hx - hy) / 2, Δ=hz)
        for k in range(-2π, 2π; length=97)
            paper = sqrt((hx * cos(k / 2))^2 + (hy * sin(k / 2))^2 + hz^2)
            worst_map = max(
                worst_map, abs(QAtlas._rice_mele_dispersion(k, m.v, m.w, m.Δ) - paper)
            )
        end
    end
    @test worst_map < 1.0e-13           # measured 2.2e-16
end

@testset "RiceMele — closed forms with a witness" begin
    @testset "flat-band sweet spot w = 0" begin
        # Decoupled dimers: E± = ±√(v²+Δ²) for every k, so ε₀ = −√(v²+Δ²)/2 exactly.
        for (v, Δ) in ((1.0, 0.0), (1.0, 0.3), (0.4, 0.9), (2.0, -0.5))
            m = RiceMele(; v, w=0.0, Δ)
            @test fetch(m, Energy{:per_site}(), Infinite()) ≈ -sqrt(v^2 + Δ^2) / 2 atol =
                1.0e-12
            @test fetch(m, MassGap(), Infinite()) ≈ sqrt(v^2 + Δ^2) atol = 1.0e-14
        end
    end

    @testset "MassGap@OBC does not converge to Infinite when |w| > |v|" begin
        # The edge-mode scale, not the bulk gap — SSH has the same caveat, and the docstring
        # says so. Flat in N, so this is not a slow transient.
        m = RiceMele(; v=0.3, w=0.9, Δ=-0.4)
        obc = [fetch(m, MassGap(), OBC(N)) for N in (10, 40, 200)]
        @test all(x -> isapprox(x, obc[1]; atol=1.0e-8), obc)
        @test !isapprox(obc[1], fetch(m, MassGap(), Infinite()); atol=1.0e-2)
        # Trivial side |v| > |w|: there it does converge, polynomially — residual 4.6e-5 at
        # N = 200, against the 0.32 the topological side never loses.
        t = RiceMele(; v=0.9, w=0.3, Δ=-0.4)
        @test fetch(t, MassGap(), OBC(200)) ≈ fetch(t, MassGap(), Infinite()) atol = 1.0e-3
    end

    @testset "OBC energy density converges to the Infinite one as O(1/N)" begin
        m = RiceMele(; v=0.75, w=0.25, Δ=0.1)
        e∞ = fetch(m, Energy{:per_site}(), Infinite())
        errs = [abs(fetch(m, Energy{:per_site}(), OBC(N)) - e∞) for N in (20, 40, 80, 160)]
        @test issorted(errs; rev=true)                 # it converges at all
        # Halving 1/N should halve the error; allow a factor of two of slack either way.
        ratios = errs[1:(end - 1)] ./ errs[2:end]
        @test all(1.0 .< ratios .< 4.0)
        # The ratio bound alone passes for a sequence that plateaus at the WRONG limit
        # (ratios → 1⁺). Anchor the value too: measured 8.4e-6 at N = 1280 for the correct
        # implementation, 2.6e-3 for a Δ/2 confined to the OBC diagonal.
        @test abs(fetch(m, Energy{:per_site}(), OBC(1280)) - e∞) < 1.0e-4
    end

    @testset "the (hx, hy, hz) parametrisation of the response literature" begin
        # The bridge to the nonlinear-response packages, which use this form.
        for (hx, hy, hz) in ((1.0, 0.5, 0.1), (1.0, 0.0, 0.4), (0.8, 0.8, 0.2))
            m = RiceMele(; v=(hx + hy) / 2, w=(hx - hy) / 2, Δ=hz)
            @test fetch(m, MassGap(), Infinite()) ≈ sqrt(hy^2 + hz^2) atol = 1.0e-14
        end
        # The paper's own numbers: hx=1, hy=0.5, hz=0.1 has a band gap of 1.0198.
        m = RiceMele(; v=0.75, w=0.25, Δ=0.1)
        @test 2 * fetch(m, MassGap(), Infinite()) ≈ 1.019803902718557 atol = 1.0e-12
        # Every Infinite observable here is symmetric under v ↔ w, so the mapping needs an
        # OBC witness to be sensitive to which hopping sits on the boundary.
        @test fetch(m, MassGap(), OBC(40)) ≉
            fetch(RiceMele(; v=0.25, w=0.75, Δ=0.1), MassGap(), OBC(40))
    end
end

@testset "RiceMele — structural and error guards" begin
    @test_throws ArgumentError RiceMele(; v=NaN, w=1.0, Δ=0.0)
    @test_throws ArgumentError RiceMele(; v=1.0, w=Inf, Δ=0.0)
    @test_throws ArgumentError RiceMele(; v=1.0, w=1.0, Δ=NaN)
    # `OBC(0)` is caught by the boundary layer first, so the model's own `N ≥ 1` guard is
    # unreachable through `fetch` and is tested where it can fire.
    @test_throws ErrorException fetch(RiceMele(), ExactSpectrum(), OBC(0))
    @test_throws ArgumentError QAtlas._rice_mele_obc_spectrum(0, 1.0, 1.0, 0.0)
    @test_throws ArgumentError QAtlas._rice_mele_obc_spectrum_sites(0, 1.0, 1.0, 0.0)

    m = RiceMele(; v=0.6, w=1.0, Δ=0.25)
    @test RiceMele() === RiceMele(; v=1.0, w=1.0, Δ=0.0)      # Δ defaults to the SSH point
    @test fetch(m, CorrelationLength(), Infinite()) ≈ 1 / fetch(m, MassGap(), Infinite())
    # Sign conventions: the answer cannot depend on the sign of Δ, nor on swapping the
    # two hoppings, since both are relabellings of the same chain.
    @test fetch(m, Energy{:per_site}(), Infinite()) ≈
        fetch(RiceMele(; v=0.6, w=1.0, Δ=-0.25), Energy{:per_site}(), Infinite())
    @test fetch(m, Energy{:per_site}(), Infinite()) ≈
        fetch(RiceMele(; v=1.0, w=0.6, Δ=0.25), Energy{:per_site}(), Infinite())
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────────
@testset "RiceMele — verification cards" begin
    # MassGap@Infinite against a brute-force minimum of E₊(k) over a fine grid — no
    # reference to the √((|v|−|w|)²+Δ²) expression, so a wrong closed form is visible.
    for (v, w, Δ) in (
        (0.75, 0.25, 0.1),
        (0.25, 0.75, 0.1),
        (1.0, 1.0, 0.3),
        (-0.5, 0.7, 0.2),
        (0.6, -1.0, -0.4),
        (1.0, 0.0, 0.5),
    )
        ks = range(-π, π; length=400_001)
        gap_num = minimum(sqrt(v^2 + w^2 + 2 * v * w * cos(k) + Δ^2) for k in ks)
        verify(
            RiceMele(; v, w, Δ),
            MassGap(),
            Infinite();
            route=:second_closed_form,
            independent=gap_num,
            agree_within=1.0e-9,
            refs=["Rice & Mele 1982: E₊(k) = √(v²+w²+2vw cos k + Δ²), gap = min_k E₊"],
        )
        verify(
            RiceMele(; v, w, Δ),
            CorrelationLength(),
            Infinite();
            route=:second_closed_form,
            independent=1 / gap_num,
            agree_within=1.0e-8,
            refs=["ξ = 1/gap, with the gap taken from the same grid minimum"],
        )
    end

    # Energy@Infinite against a hand-written Simpson rule — the source uses adaptive
    # Gauss-Kronrod, so the quadrature is a different one, not the same call re-typed.
    for (v, w, Δ) in ((0.75, 0.25, 0.1), (1.0, 1.0, 0.0), (0.4, 1.1, -0.6))
        n = 200_000
        h = 2π / n
        f(k) = sqrt(v^2 + w^2 + 2 * v * w * cos(k) + Δ^2)
        simpson =
            f(-π) +
            f(π) +
            4 * sum(f(-π + (2j - 1) * h) for j in 1:(n ÷ 2)) +
            2 * sum(f(-π + 2j * h) for j in 1:(n ÷ 2 - 1))
        verify(
            RiceMele(; v, w, Δ),
            Energy{:per_site}(),
            Infinite();
            route=:second_closed_form,
            independent=(-(simpson * h / 3) / (4π)),
            agree_within=1.0e-9,
            refs=["Rice & Mele 1982: ε₀ = −(1/4π)∫E₊(k)dk at half filling"],
        )
    end

    # The fully dimerised limits are exact: w = 0 decouples the chain into N dimers of
    # energy ±√(v²+Δ²), and v = 0 into N−1 dimers plus two isolated ±Δ end sites.
    for (v, Δ) in ((1.0, 0.0), (1.0, 0.3), (0.4, 0.9), (2.0, -0.5))
        verify(
            RiceMele(; v, w=0.0, Δ),
            Energy{:per_site}(),
            Infinite();
            route=:limiting_case,
            independent=(-sqrt(v^2 + Δ^2) / 2),
            agree_within=1.0e-10,
            refs=["w = 0: flat band E± = ±√(v²+Δ²), so ε₀ = −√(v²+Δ²)/2"],
        )
        for N in (4, 9)
            verify(
                RiceMele(; v, w=0.0, Δ),
                Energy{:per_site}(),
                OBC(N);
                route=:limiting_case,
                independent=(-sqrt(v^2 + Δ^2) / 2),
                agree_within=1.0e-10,
                refs=["w = 0: N decoupled dimers, no boundary term at any N"],
            )
            # `verify` compares scalars; the full array is pinned by the plain @test above.
            verify(
                RiceMele(; v, w=0.0, Δ),
                ExactSpectrum(),
                OBC(N);
                route=:limiting_case,
                independent=2N * sqrt(v^2 + Δ^2),
                agree_within=1.0e-10,
                subject_extract=s -> sum(abs, s),
                refs=["w = 0: the 2N levels are ±√(v²+Δ²), so Σ|ε| = 2N√(v²+Δ²)"],
            )
            verify(
                RiceMele(; v, w=0.0, Δ),
                MassGap(),
                OBC(N);
                route=:limiting_case,
                independent=sqrt(v^2 + Δ^2),
                agree_within=1.0e-10,
                refs=["w = 0: (ε_{N+1} − ε_N)/2 = √(v²+Δ²), N-independent"],
            )
        end
    end
end
