# The finite-L Calabrese-Cardy block entropy (core/cft_entanglement.jl).
#
# Every assertion here is anchored to an EXACT route the atlas already has, not
# to a literature number: TFIM/OBC is Peschel's correlation-matrix method at any
# L, XXZ1D/OBC is dense ED.  That is the legitimate use of an exponential row —
# the instrument a closed form is validated against.

using Test
using QAtlas
using QAtlas: chord_distance, cft_block_entropy, fetch, OBC, Infinite
using QAtlas: TFIM, XXZ1D, VonNeumannEntropy, LuttingerParameter
using QAtlas: cft_region_entropy, FermionicEntanglementEntropy, Region

# least-squares offset only (the non-universal constant this file does not return)
function _fit_offset(exact, pred)
    c1 = sum(exact .- pred) / length(exact)
    r = exact .- pred .- c1
    return c1, sqrt(sum(abs2, r) / length(r))
end

@testset "chord_distance refuses a block with no complement" begin
    @test chord_distance(4, 8) ≈ (16 / π) * sin(π / 2)
    @test chord_distance(4, 8; periodic=true) ≈ (8 / π) * sin(π / 2)
    # the open form is exactly twice the periodic one at the same (ℓ, L)
    @test chord_distance(3, 10) ≈ 2 * chord_distance(3, 10; periodic=true)
    @test_throws ArgumentError chord_distance(8, 8)     # whole chain: pure state
    @test_throws ArgumentError chord_distance(0, 8)
end

@testset "the open prefactor is c/6 and the periodic one c/3" begin
    # ratio of the two leading terms at the SAME chord distance isolates the
    # prefactor, which is the one thing the formula must not get wrong
    d_open = chord_distance(4, 16)
    open_val = cft_block_entropy(0.5, 4, 16)
    @test open_val ≈ (0.5 / 6) * log(d_open)
    per_val = cft_block_entropy(0.5, 4, 16; periodic=true)
    @test per_val ≈ (0.5 / 3) * log(chord_distance(4, 16; periodic=true))
end

@testset "reproduces TFIM's EXACT correlation-matrix entropies (c = 1/2)" begin
    # No alternating term for the free-fermion Ising chain — measured rms ~0.002
    # at both sizes, which is what "the plain chord formula is enough here" means.
    for L in (32, 64)
        ℓs = 4:2:(L ÷ 2)
        exact = [
            fetch(TFIM(; J=1.0, h=1.0), VonNeumannEntropy(), OBC(L); ℓ=ℓ, beta=Inf) for
            ℓ in ℓs
        ]
        pred = [cft_block_entropy(0.5, ℓ, L) for ℓ in ℓs]
        c1, rms = _fit_offset(exact, pred)
        @test rms < 0.005
        # the non-universal constant is a property of model+boundary, NOT of L:
        # its stability across L is the evidence the leading term is right
        @test 0.20 < c1 < 0.27
    end
end

@testset "the alternating term is needed at c = 1, and its SIGN is fixed" begin
    # L = 10, not 12, and the difference is the whole point of this file: the
    # free-fermion route does L = 64 with 15 blocks in 0.09 s, while dense ED
    # takes 5.05 s at L = 10 and 108 s at L = 12 for five.  MEASURED that the fit
    # is unchanged — A = 0.180 vs 0.174 at Δ = 0, 0.621 vs 0.650 at Δ = 1 — so
    # the larger chain buys nothing here except twenty times the runtime.
    L, ℓs = 10, 2:5
    for (Δ, Amin, Amax) in ((0.0, 0.10, 0.25), (1.0, 0.45, 0.85))
        exact = [
            fetch(XXZ1D(; Δ=Δ), VonNeumannEntropy(), OBC(L); ℓ=ℓ, beta=Inf) for ℓ in ℓs
        ]
        K = fetch(XXZ1D(; Δ=Δ), LuttingerParameter(), Infinite())
        _, rms_plain = _fit_offset(exact, [cft_block_entropy(1.0, ℓ, L) for ℓ in ℓs])
        # joint least squares over (constant, amplitude)
        d = [chord_distance(ℓ, L) for ℓ in ℓs]
        lead = [(1 / 6) * log(x) for x in d]
        basis = [(-1)^(ℓ + 1) * x^(-1 / (2K)) for (ℓ, x) in zip(ℓs, d)]
        coef = hcat(ones(length(ℓs)), basis) \ (exact .- lead)
        A = coef[2]
        _, rms_alt = _fit_offset(
            exact, [cft_block_entropy(1.0, ℓ, L; K=K, alt_amplitude=A) for ℓ in ℓs]
        )
        # the plain formula is BAD here — 20% of the leading log — which is the
        # reason the term exists at all; assert that rather than only the fix
        @test rms_plain > 0.05
        @test rms_alt < rms_plain / 3
        # POSITIVE amplitude, i.e. odd blocks sit ABOVE the smooth curve.  Writing
        # (-1)^ℓ instead of (-1)^(ℓ+1) flips this and makes the fit worse; that is
        # how the sign was caught, so it is pinned here.
        @test Amin < A < Amax
    end
end

@testset "the non-universal amplitude is not silently invented" begin
    # no K, or a periodic chain, with a requested oscillation is an ERROR: the
    # exponent needs K, and the oscillation is a property of the open boundary
    @test_throws ArgumentError cft_block_entropy(1.0, 3, 12; alt_amplitude=0.5)
    @test_throws ArgumentError cft_block_entropy(
        1.0, 3, 12; periodic=true, K=1.0, alt_amplitude=0.5
    )
    @test_throws ArgumentError cft_block_entropy(1.0, 3, 12; K=-1.0, alt_amplitude=0.5)
    # and the default really is "leading term only", matching the convention the
    # infinite-chain forms already follow
    @test cft_block_entropy(1.0, 3, 12; K=1.0) == cft_block_entropy(1.0, 3, 12)
end

# ─────────────────────────────────────────────────────────────────────────────
# Arbitrary regions — the signed pair sum (cft_region_entropy)
# ─────────────────────────────────────────────────────────────────────────────

# I(A:B) = S(A) + S(B) - S(A∪B).  Every non-universal constant and the cutoff
# cancel here, so the comparison below has NOTHING fitted — unlike the block
# tests above, which fit an offset.
_mi(f, A, B) = f(A) + f(B) - f(sort(vcat(A, B)))

@testset "cft_region_entropy reduces to cft_block_entropy on ONE interval" begin
    # The whole claim of the pair sum is that the region enters only through its
    # endpoints, so a single interval must give back the bipartite formula
    # identically — not approximately.
    for L in (16, 33, 64), ℓ in (1, 2, 5, L - 1)
        for periodic in (false, true)
            @test cft_region_entropy(1.0, ℓ, L; periodic=periodic) ≈
                cft_block_entropy(1.0, ℓ, L; periodic=periodic) rtol = 1e-12
            @test cft_region_entropy(0.5, Region(1:ℓ...), L; periodic=periodic) ≈
                cft_block_entropy(0.5, ℓ, L; periodic=periodic) rtol = 1e-12
        end
    end
    # linear in c, like the kernel it is built from
    @test cft_region_entropy(2.0, [2, 3, 7, 8], 20) ≈
        2 * cft_region_entropy(1.0, [2, 3, 7, 8], 20)
end

@testset "cft_region_entropy vs the EXACT free-fermion route (nothing fitted)" begin
    # TFIM at criticality: c = 1/2, and the Peschel covariance restricted to any
    # site set is the FERMIONIC entropy — the quantity the pair sum predicts.
    m = TFIM(1.0, 1.0)
    L = 96
    Sf(sites) = fetch(m, FermionicEntanglementEntropy(), OBC(L); region=Region(sites...))
    Sp(sites) = cft_region_entropy(0.5, sites, L)

    # two blocks, several separations
    for (ℓ, g) in ((6, 3), (6, 6), (6, 12), (10, 5), (10, 20))
        p = (L - (2ℓ + g)) ÷ 2 + 1
        A = collect(p:(p + ℓ - 1))
        B = collect((p + ℓ + g):(p + ℓ + g + ℓ - 1))
        exact, pred = _mi(Sf, A, B), _mi(Sp, A, B)
        @test exact > 0                                  # subadditivity, as a sanity floor
        @test isapprox(exact, pred; rtol=0.05)
        # …and the agreement is much better than that bound on the well-separated
        # configurations, where the 1/ℓ² lattice correction is smallest
        g ≥ 2ℓ && @test isapprox(exact, pred; rtol=0.02)
    end

    # three blocks: S(A)+S(B)+S(C)-S(A∪B∪C) is constant-free the same way
    ℓ, g = 6, 6
    p = (L - (3ℓ + 2g)) ÷ 2 + 1
    A = collect(p:(p + ℓ - 1))
    B = collect((p + ℓ + g):(p + 2ℓ + g - 1))
    C = collect((p + 2ℓ + 2g):(p + 3ℓ + 2g - 1))
    all3 = sort(vcat(A, B, C))
    ex3 = Sf(A) + Sf(B) + Sf(C) - Sf(all3)
    pr3 = Sp(A) + Sp(B) + Sp(C) - Sp(all3)
    @test isapprox(ex3, pr3; rtol=0.05)

    # THE RESIDUAL IS A LATTICE ARTIFACT, not a theory-specific F(x): at fixed
    # cross ratio it must SHRINK as the blocks grow.  A cross-ratio-dependent
    # term would sit at a constant relative size instead.  (Measured decay is
    # 1/ℓ², a factor of 4 per doubling; asserted loosely as ">2x better".)
    function _relerr(Lx)
        ℓx = Lx ÷ 16
        px = (Lx - 3ℓx) ÷ 2 + 1
        Ax = collect(px:(px + ℓx - 1))
        Bx = collect((px + 2ℓx):(px + 3ℓx - 1))
        f(s) = fetch(m, FermionicEntanglementEntropy(), OBC(Lx); region=Region(s...))
        e, q = _mi(f, Ax, Bx), _mi(s -> cft_region_entropy(0.5, s, Lx), Ax, Bx)
        return abs(e - q) / abs(q)
    end
    e64, e128, e256 = _relerr(64), _relerr(128), _relerr(256)
    @test e64 > 2 * e128 > 4 * e256
    @test e256 < 2e-3
end

@testset "the spin and fermionic routes agree ONLY on one interval" begin
    m = TFIM(1.0, 1.0)
    bc = OBC(12)
    for ℓ in 1:6
        r = Region((1:ℓ)...)
        @test fetch(m, VonNeumannEntropy(), bc; region=r) ≈
            fetch(m, FermionicEntanglementEntropy(), bc; region=r) rtol = 1e-12
    end
    # a block in the BULK is still one interval, so still both
    r = Region(4, 5, 6)
    @test fetch(m, VonNeumannEntropy(), bc; region=r) ≈
        fetch(m, FermionicEntanglementEntropy(), bc; region=r) rtol = 1e-12
    # a disconnected region: the spin route REFUSES rather than return the
    # fermionic number under the spin name, and the fermionic route answers
    @test_throws ArgumentError fetch(m, VonNeumannEntropy(), bc; region=Region(1, 2, 5, 6))
    @test fetch(m, FermionicEntanglementEntropy(), bc; region=Region(1, 2, 5, 6)) > 0
end

@testset "cft_region_entropy: geometry the pair sum must get right" begin
    L = 40
    # A boundary-touching block merges with its own MIRROR image, so it has one
    # endpoint pair while a bulk block of the same length has two.  If the merge
    # were skipped the boundary block would pick up a log(0).
    @test isfinite(cft_region_entropy(1.0, 1:5, L))
    @test isfinite(cft_region_entropy(1.0, 36:40, L))          # the far boundary
    @test cft_region_entropy(1.0, 1:5, L) ≈ cft_region_entropy(1.0, 36:40, L) rtol = 1e-12
    # reflection symmetry of the whole construction, on a two-block region
    @test cft_region_entropy(1.0, [3, 4, 5, 12, 13, 14], L) ≈
        cft_region_entropy(1.0, [L + 1 - s for s in [3, 4, 5, 12, 13, 14]], L) rtol = 1e-12
    # a bulk block sits BELOW the same-length boundary block: two endpoints cost
    # more than one
    @test cft_region_entropy(1.0, 16:20, L) > cft_region_entropy(1.0, 1:5, L)

    # PERIODIC: sites L and 1 are neighbours, so a region holding both ends is
    # ONE arc across the seam.  Rotating a block around the ring must not change
    # its entropy — which is exactly what a missed seam merge would break.
    ref = cft_region_entropy(1.0, 5:9, L; periodic=true)
    @test cft_region_entropy(1.0, [39, 40, 1, 2, 3], L; periodic=true) ≈ ref rtol = 1e-12
    @test cft_region_entropy(1.0, vcat(38:40, 1:2), L; periodic=true) ≈ ref rtol = 1e-12
    @test cft_region_entropy(1.0, 1:5, L; periodic=true) ≈ ref rtol = 1e-12
end

@testset "cft_region_entropy: refusals rather than silent answers" begin
    @test_throws ArgumentError cft_region_entropy(1.0, Int[], 12)          # empty
    @test_throws ArgumentError cft_region_entropy(1.0, 1:12, 12)           # no complement
    @test_throws ArgumentError cft_region_entropy(1.0, [0, 1], 12)         # outside 1:L
    @test_throws ArgumentError cft_region_entropy(1.0, [11, 13], 12)       # outside 1:L
    @test_throws ArgumentError cft_region_entropy(1.0, 1:1, 1)             # L must exceed 1
end
