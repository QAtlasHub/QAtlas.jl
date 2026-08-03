# ─────────────────────────────────────────────────────────────────────────────
# core/cft_entanglement.jl — the finite-chain Calabrese-Cardy entanglement forms.
#
# The infinite-chain forms live per-model (e.g. TFIM_cft_entanglement.jl); this
# file carries the FINITE-L ones, which are model-independent given `c` and, for
# the alternating term, the Luttinger parameter `K`.
#
# For a critical chain of L sites, a block of the first ℓ, the universal leading
# term is a log of the CHORD (conformal) distance
#
#     d(ℓ, L) = (2L/π) sin(πℓ/L)        open boundaries
#     d(ℓ, L) = ( L/π) sin(πℓ/L)        periodic boundaries
#
# with prefactor c/6 for OPEN and c/3 for PERIODIC (one boundary point versus
# two), plus a non-universal constant this file does NOT return — the same
# convention the infinite-chain forms already follow: "QAtlas returns the leading
# log term only; downstream fits should account for the offset."
#
# THE ALTERNATING TERM.  For a Luttinger liquid with open boundaries the block
# entropy carries a parity oscillation on top of the log,
#
#     S_alt(ℓ) = (-1)^(ℓ+1) · A · d(ℓ, L)^(-1/(2K)),
#
# SIGN, fixed by measurement rather than convention: a POSITIVE amplitude raises
# ODD blocks and lowers even ones.  That is the direction the data goes — the
# odd-ℓ residual is positive on both chains measured (XX at L = 32: +0.070 at
# ℓ = 5; Heisenberg at L = 12: +0.159 at ℓ = 3) — and it is the physically
# expected one for an antiferromagnetic open chain, where an odd block leaves an
# unpaired spin and carries MORE entanglement.  Writing `(-1)^ℓ` instead makes
# the fit worse, not better: MEASURED rms 0.116 -> 0.271 at Δ = 1, which is how
# the sign error was caught.
#
# whose EXPONENT is universal (fixed by K) and whose AMPLITUDE A is not.  The
# exponent is derived here from the model's `LuttingerParameter`; `A` is a
# caller-supplied kwarg defaulting to 0, exactly as the additive constant is
# left to the caller.
#
# MEASURED, so the caller has a starting point rather than a bare parameter:
#
#   * the EXPONENT is right.  On the open XX chain (c = 1, K = 1) the odd-ℓ
#     residual times d^(1/2K) converges to a CONSTANT — 0.220 at d = 10, 0.157 at
#     d = 40, 0.1470 at d = 81 — verified on a free-fermion correlation-matrix
#     calculation up to L = 128.  If the exponent were wrong that product would
#     drift without limit instead of settling.
#   * the AMPLITUDE, by joint least squares over (constant, A) against the atlas's
#     own exact routes at L = 12:
#
#         XXZ Δ = 0 (K = 1)    A ≈ 0.174    rms 0.0707 -> 0.0051   (14x better)
#         XXZ Δ = 1 (K = 1/2)  A ≈ 0.650    rms 0.1164 -> 0.0275   (4.2x better)
#
#     The K = 1 amplitude drifts from 0.174 (fit at L = 12) to ≈ 0.147 (large-L
#     free-fermion data), which is the L = 12 fit absorbing subleading corrections
#     — so treat these as starting points, not constants of nature.  The K = 1/2
#     value rests on ED at N ≤ 12 alone, since there is no cheaper exact route for
#     the interacting chain; it is the least determined number here.
#
# WHY THE ALTERNATING TERM MATTERS.  Without it the plain chord formula is wrong
# by ~20% of the leading log at accessible sizes for c = 1: MEASURED residuals of
# 0.09-0.16 at L = 12, alternating in sign with the parity of ℓ.  For the free-
# fermion Ising chain (c = 1/2) there is no such term to speak of — the plain
# formula reproduces the exact BdG entropies with rms 0.002 at L = 32 and 64.
#
# References: Calabrese & Cardy, [CalabreseCardy2004](@cite) (the log and the
# open/periodic prefactors); the parity oscillation in open Luttinger chains is
# standard in the entanglement literature and its amplitude is fitted, not
# derived, which is why this file does not ship one.
# ─────────────────────────────────────────────────────────────────────────────

"""
    chord_distance(ℓ, L; periodic=false) -> Float64

Conformal (chord) distance of a block of `ℓ` sites in a critical chain of `L`:
`(2L/π) sin(πℓ/L)` open, `(L/π) sin(πℓ/L)` periodic.

Throws unless `0 < ℓ < L`, so a block that is the whole chain — where the
entropy is that of a pure state and the formula does not apply — is refused
rather than returning `log(0)`.
"""
function chord_distance(ℓ::Integer, L::Integer; periodic::Bool=false)
    0 < ℓ < L || throw(
        ArgumentError(
            "chord_distance: need 0 < ℓ < L (got ℓ = $ℓ, L = $L); a block equal " *
            "to the whole chain has no complement to be entangled with",
        ),
    )
    return (periodic ? L : 2L) / π * sin(π * ℓ / L)
end

"""
    cft_block_entropy(c, ℓ, L; periodic=false, K=nothing, alt_amplitude=0.0) -> Float64

Universal Calabrese-Cardy block entanglement entropy of the first `ℓ` sites of a
critical chain of `L` sites with central charge `c`:

    S(ℓ) = (c/6) log d(ℓ, L)  +  (-1)^(ℓ+1) · A · d(ℓ, L)^(-1/(2K))  open
    S(ℓ) = (c/3) log d(ℓ, L)  +  …                                    periodic

The **non-universal additive constant is not included** — the convention the
infinite-chain forms already use.

`K` is the Luttinger parameter; supply it together with `alt_amplitude` to get
the parity oscillation of an open Luttinger chain.  The exponent `1/(2K)` is
universal; the amplitude is not, which is why it has no default (see the file
header for measured values).  A POSITIVE amplitude raises ODD blocks.

Passing `alt_amplitude` without `K`, or with `periodic = true`, is an error
rather than a silent no-op: the oscillation is a property of the open boundary
and its exponent needs `K`.
"""
function cft_block_entropy(
    c::Real,
    ℓ::Integer,
    L::Integer;
    periodic::Bool=false,
    K::Union{Nothing,Real}=nothing,
    alt_amplitude::Real=0.0,
)
    d = chord_distance(ℓ, L; periodic=periodic)
    leading = (periodic ? c / 3 : c / 6) * log(d)
    iszero(alt_amplitude) && return leading
    periodic && throw(
        ArgumentError(
            "cft_block_entropy: the parity oscillation is a property of the OPEN " *
            "boundary; drop `alt_amplitude` for periodic, or set periodic = false",
        ),
    )
    K === nothing && throw(
        ArgumentError(
            "cft_block_entropy: `alt_amplitude` needs `K` — the oscillation's " *
            "exponent is 1/(2K) and only the AMPLITUDE is non-universal",
        ),
    )
    K > 0 || throw(ArgumentError("cft_block_entropy: K must be positive; got $K"))
    return leading + (-1)^(ℓ + 1) * alt_amplitude * d^(-1 / (2K))
end

# ─────────────────────────────────────────────────────────────────────────────
# ARBITRARY REGIONS — the same kernel, summed over endpoint PAIRS.
#
# For the free massless Dirac field in two dimensions the mutual information is
# EXTENSIVE (Casini & Huerta, [CasiniHuerta2009](@cite); the modular Hamiltonian
# that underlies it is [CasiniHuerta2009Modular](@cite)), and the entanglement
# entropy of a union of intervals is then a SIGNED SUM of the bipartite chord
# kernel over the endpoint set: with arcs `[aᵢ, bᵢ]` on a circle of circumference
# `C` and chord `D(x,y) = (C/π)|sin(π(x−y)/C)|`,
#
#     pairsum = Σ_{i,j} ln D(aᵢ,bⱼ) − Σ_{i<j} ln D(aᵢ,aⱼ) − Σ_{i<j} ln D(bᵢ,bⱼ).
#
# The region enters ONLY through its endpoints: neither the number of intervals
# nor their arrangement appears anywhere else, which is what makes an arbitrary
# region derivable from the bipartite result rather than a separate calculation.
# One interval gives exactly one term, `ln D(a,b)`, so `cft_region_entropy`
# reduces identically to `cft_block_entropy` — asserted in the tests, not
# asserted here.
#
# THE OPEN CHAIN IS THE CHIRAL HALF of a periodic system of circumference 2L, in
# which the region is accompanied by its MIRROR IMAGE −A.  Hence `C = 2L`, the
# arcs are `A ∪ (−A)`, and the prefactor is c/6 = ½ · c/3.  A block touching a
# boundary merges with its own image (`[0,ℓ]` and `[−ℓ,0]` become `[−ℓ,ℓ]`),
# which is why the boundary-touching block has ONE endpoint pair and recovers
# `(c/6) ln[(2L/π) sin(πℓ/L)]`.
#
# MEASURED against the exact free-fermion route (mutual information, so every
# non-universal constant and the cutoff cancel and there is nothing to fit):
#
#   * the residual is a LATTICE artifact, not a theory-specific term.  At fixed
#     cross ratio it falls by 4x per doubling of the block length — 1/ℓ²:
#
#         ℓ = 4    XX −6.13%    TFIM −1.41%
#         ℓ = 16   XX −0.34%    TFIM −0.09%
#         ℓ = 64   XX −0.03%    TFIM −0.01%
#
#   * across cross ratios 0.04 … 0.94 the ratio I(TFIM)/I(XX) sits at 0.4996 …
#     0.5002, i.e. exactly c_Ising/c_Dirac with no drift.  A theory-specific
#     F(x) would bend it.  `c` is the only model input.
#
#     Worth flagging, since it is broader than what is cited: the references
#     above are the DIRAC field, and the critical TFIM is a MAJORANA one.  The
#     agreement there is measured here, not taken from them.
#
# WHAT IT PREDICTS is the `FermionicEntanglementEntropy`, NOT the spin entropy —
# see that quantity's docstring, and `core/regions.jl` for the Jordan-Wigner
# reason.  The two agree on a single contiguous interval and nowhere else.
#
# NO ALTERNATING TERM here.  The parity oscillation of `cft_block_entropy` is
# indexed by the parity of one block length; a region with several intervals has
# no single such parity, and the multi-interval amplitude is not the bipartite
# one.  `alt_amplitude` is refused rather than applied to some chosen interval.
# ─────────────────────────────────────────────────────────────────────────────

# maximal runs of consecutive sites, as (first, last) integer pairs
function _site_runs(sites::AbstractVector{<:Integer})
    runs = Tuple{Int,Int}[]
    isempty(sites) && return runs
    p = sites[1]
    for k in 2:length(sites)
        if sites[k] != sites[k - 1] + 1
            push!(runs, (p, sites[k - 1]))
            p = sites[k]
        end
    end
    push!(runs, (p, sites[end]))
    return runs
end

# Site run `p:q` occupies the continuum interval `[p−1, q]` of a chain laid out
# on `[0, L]`, so the leading block `1:ℓ` is `[0, ℓ]` — the convention that makes
# the boundary-touching case agree with `chord_distance`.
#
# Returned arcs are `(left, right)` traversed in ONE direction; the pair sum
# distinguishes left from right endpoints, so an arc that crosses the
# identification point is written with a right endpoint beyond the period (e.g.
# `(a, 2L−a)` for a block touching the far boundary) rather than wrapped.  The
# chord is periodic, so out-of-range representatives are exact, not approximate.
function _region_arcs(runs::Vector{Tuple{Int,Int}}, L::Int, periodic::Bool)
    if periodic
        arcs = [(float(p - 1), float(q)) for (p, q) in runs]
        # on a ring, site L and site 1 are neighbours: a region holding both ends
        # is ONE arc across the seam, not two
        if length(arcs) > 1 && first(runs[1]) == 1 && last(runs[end]) == L
            a_first, b_first = arcs[1]
            a_last, _ = arcs[end]
            arcs = vcat(arcs[2:(end - 1)], [(a_last, b_first + L)])
        end
        return arcs
    end
    arcs = Tuple{Float64,Float64}[]
    for (p, q) in runs
        a, b = float(p - 1), float(q)
        if a == 0.0                        # touches the near boundary: merges with −A
            push!(arcs, (-b, b))
        elseif b == float(L)               # touches the far boundary: merges across ±L
            push!(arcs, (a, 2L - a))
        else
            push!(arcs, (a, b))
            push!(arcs, (-b, -a))
        end
    end
    return arcs
end

function _chord_pairsum(arcs::Vector{Tuple{Float64,Float64}}, C::Real)
    as = first.(arcs)
    bs = last.(arcs)
    n = length(arcs)
    d(x, y) = (C / π) * abs(sin(π * (x - y) / C))
    lg = function (x, y, what)
        v = d(x, y)
        v > 0 || throw(
            ArgumentError(
                "cft_region_entropy: two $what coincide on the conformal circle " *
                "(x = $x, y = $y), so the pair sum diverges; this means the " *
                "region's intervals were not maximal or the region fills the chain",
            ),
        )
        return log(v)
    end
    s = 0.0
    for x in as, y in bs
        s += lg(x, y, "endpoints")
    end
    for i in 1:(n - 1), j in (i + 1):n
        s -= lg(as[i], as[j], "left endpoints")
        s -= lg(bs[i], bs[j], "right endpoints")
    end
    return s
end

"""
    cft_region_entropy(c, region, L; periodic=false) -> Float64

Universal Calabrese-Cardy entanglement entropy of an **arbitrary** `region` of a
critical chain of `L` sites with central charge `c` — the multi-interval
generalisation of [`cft_block_entropy`](@ref), obtained by summing the same
chord kernel over the region's endpoint PAIRS with alternating signs
(Casini & Huerta, [CasiniHuerta2009](@cite)).

`region` is an AbstractQAtlas `Region`, any iterable of sites, or an `Integer`
block length `ℓ` (meaning `1:ℓ`), matching the `region`/`ℓ` sugar of the
entanglement `fetch` methods.  For a single interval this returns exactly
`cft_block_entropy(c, ℓ, L; periodic)`.

As there, the **non-universal additive constant is not included** — here it is
one constant per endpoint, so it cancels outright in a mutual information.

This predicts the [`FermionicEntanglementEntropy`](@ref) of the region.  For a
disconnected region that is **not** the spin `VonNeumannEntropy`: the
Jordan-Wigner string leaves the region, and the two differ by O(0.1) nats
without violating any entropy inequality (see `core/regions.jl`).  They coincide
on a single contiguous interval.

There is no `alt_amplitude` here — see the section header.
"""
function cft_region_entropy(c::Real, region, L::Integer; periodic::Bool=false)
    L > 1 || throw(ArgumentError("cft_region_entropy: need L > 1; got $L"))
    sites = _region_sites(region)
    isempty(sites) && throw(ArgumentError("cft_region_entropy: the region is empty."))
    (first(sites) ≥ 1 && last(sites) ≤ L) || throw(
        ArgumentError(
            "cft_region_entropy: region sites must lie in 1:L; got $sites with L = $L."
        ),
    )
    length(sites) ≤ L - 1 || throw(
        ArgumentError(
            "cft_region_entropy: the region must leave a non-empty complement " *
            "(the entropy of the whole pure chain is 0, and the formula does not " *
            "apply); got |region| = $(length(sites)) with L = $L.",
        ),
    )
    arcs = _region_arcs(_site_runs(sites), Int(L), periodic)
    circumference = periodic ? L : 2L
    return (periodic ? c / 3 : c / 6) * _chord_pairsum(arcs, circumference)
end
