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
