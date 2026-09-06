# ─────────────────────────────────────────────────────────────────────────────
# Test: Rice-Mele (1982) dimerised chain with a staggered on-site potential.
#
# The oracle is SSH: at Δ = 0 this model IS SSH, and SSH is an independent
# implementation here that shares no code with it. The rest is closed forms with
# a witness, plus the guards a copy of SSH would fail.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test
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
        @test maximum(abs.(odd + reverse(odd))) > 0.1
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
        @test_throws Exception fetch(m, QAtlas.TopologicalInvariant(), Infinite())
    end
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

    @testset "OBC energy density converges to the Infinite one as O(1/N)" begin
        m = RiceMele(; v=0.75, w=0.25, Δ=0.1)
        e∞ = fetch(m, Energy{:per_site}(), Infinite())
        errs = [abs(fetch(m, Energy{:per_site}(), OBC(N)) - e∞) for N in (20, 40, 80, 160)]
        @test issorted(errs; rev=true)                 # it converges at all
        # Halving 1/N should halve the error; allow a factor of two of slack either way.
        ratios = errs[1:(end - 1)] ./ errs[2:end]
        @test all(1.0 .< ratios .< 4.0)
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
