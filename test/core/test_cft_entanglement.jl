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
