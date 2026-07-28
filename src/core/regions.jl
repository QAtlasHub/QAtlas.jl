# ─────────────────────────────────────────────────────────────────────────────
# core/regions.jl — a REGION as a fetch argument for the entanglement entropies.
#
# AbstractQAtlas keys a region entropy by VALUE:
#
#     entanglement_entropy(A) == VariableKey(VonNeumannEntropy, RegionSupport(A))
#
# so `region_report` can auto-discover subadditivity, Araki–Lieb, strong
# subadditivity and weak monotonicity over a bag of entropies without any
# A/B/AB hand-labelling.  What was missing on this side is the ability to ASK a
# hub for `S(A)` at all: every implementation here took a block length `ℓ`,
# which names the one region `1:ℓ` and nothing else.
#
# This file supplies the translation.  `ℓ` stays exactly the sugar it always
# was — `ℓ = k` means `region = Region(1:k)` — so no existing call changes.
#
# THE JORDAN–WIGNER CAVEAT, and why it is a hard error rather than a footnote.
# For a free-fermion (BdG / Peschel) implementation the reduced state on a set
# of sites is obtained by restricting the Majorana covariance to that set's
# Majorana indices.  That is exact for the FERMIONS.  It equals the SPIN
# entropy only when the region is a single contiguous interval, because only
# then does the Jordan–Wigner string factorise across the boundary (the parity
# operator left on A commutes with ρ_A — Peschel 2003; Fagotti–Calabrese 2010).
# For two or more disjoint intervals the two entropies genuinely differ, and
# they differ in a way no inequality would catch: both are honest von Neumann
# entropies of honest states, so subadditivity and friends hold for either one.
# A silently-fermionic answer would therefore pass every check in the atlas
# while answering a different question than `VonNeumannEntropy` on a SPIN model
# asks.  So the free-fermion path refuses multi-interval regions outright and
# says why; the dense-ED path, which traces out the complement of the actual
# spin sites, has no such restriction and takes any region.
#
# This is not a limitation in practice for the region inequalities: choosing the
# regions ADJACENT — A = 1:2, B = 3:4, C = 5:6 — makes every entropy that
# `region_report` needs (S(A), S(B), S(C), S(A∪B), S(B∪C), S(A∪B∪C)) a single
# contiguous interval, so the free-fermion hubs materialise all four relations
# with no caveat at all.
# ─────────────────────────────────────────────────────────────────────────────

"""
    _region_sites(region) -> Vector{Int}

The sites of `region` as a sorted `Vector{Int}`.  Accepts a `Region`, any
iterable of site labels, or an `Integer` block length (the legacy `ℓ` sugar,
meaning `1:ℓ`).
"""
_region_sites(r::Region) = sort!(collect(Int, r.sites))
_region_sites(ℓ::Integer) = collect(1:Int(ℓ))
_region_sites(sites) = sort!(collect(Int, sites))

"""
    _is_contiguous(sites::AbstractVector{<:Integer}) -> Bool

Whether the sorted site list is a single unbroken interval `a:b`.
"""
function _is_contiguous(sites::AbstractVector{<:Integer})
    isempty(sites) && return true
    return last(sites) - first(sites) + 1 == length(sites)
end

"""
    _entanglement_sites(N::Int, region, ℓ) -> Vector{Int}

Resolve the `region` / `ℓ` kwarg pair into a validated, sorted site list on an
`N`-site chain.  Exactly one of the two must be given; `ℓ = k` is sugar for
`Region(1:k)`.

The region must be non-empty, lie inside `1:N`, and leave a non-empty
complement — `S(A)` with `A` the whole system is the entropy of a pure state,
which is 0 by construction and carries no information about the implementation.
"""
function _entanglement_sites(N::Int, region, ℓ)
    if region === nothing && ℓ === nothing
        throw(
            ArgumentError(
                "VonNeumannEntropy: supply either `region` (an AbstractQAtlas `Region`) " *
                "or `ℓ` (a block length, meaning the region `1:ℓ`).",
            ),
        )
    elseif region !== nothing && ℓ !== nothing
        throw(
            ArgumentError(
                "VonNeumannEntropy: `region` and `ℓ` are alternatives — `ℓ = k` IS " *
                "`region = Region(1:k)`; got both (region = $region, ℓ = $ℓ).",
            ),
        )
    end
    # A bare Integer in the `region` slot reads two ways — `region = 3` could be
    # site 3 or the block `1:3` — and both are plausible enough that guessing
    # would be a silent wrong answer for half the callers.  Refuse it and name
    # the two spellings instead.
    region isa Integer && throw(
        ArgumentError(
            "VonNeumannEntropy: `region = $region` is ambiguous — write " *
            "`region = Region($region)` for the single site, or `ℓ = $region` " *
            "(equivalently `region = Region(1:$region...)`) for the leading block.",
        ),
    )
    sites = region === nothing ? _region_sites(ℓ) : _region_sites(region)
    isempty(sites) && throw(ArgumentError("VonNeumannEntropy: the region is empty."))
    (first(sites) ≥ 1 && last(sites) ≤ N) || throw(
        ArgumentError(
            "VonNeumannEntropy: region sites must lie in 1:N; got $sites with N = $N."
        ),
    )
    length(sites) ≤ N - 1 || throw(
        ArgumentError(
            "VonNeumannEntropy: the region must leave a non-empty complement " *
            "(S of the whole pure system is 0 by construction); got |region| = " *
            "$(length(sites)) with N = $N.",
        ),
    )
    return sites
end

"""
    _require_single_interval(sites, who::AbstractString)

Throw unless `sites` is a single contiguous interval.  Called by the free-fermion
entanglement implementations, where a multi-interval region would return the
FERMIONIC entropy rather than the spin one — see this file's header.
"""
function _require_single_interval(sites::AbstractVector{<:Integer}, who::AbstractString)
    _is_contiguous(sites) && return nothing
    return throw(
        ArgumentError(
            "$who: the free-fermion (Peschel covariance) route is exact only for a " *
            "SINGLE contiguous interval, because the Jordan–Wigner string factorises " *
            "across the boundary only then; for the multi-interval region $sites it " *
            "would return the FERMIONIC entanglement entropy, which differs from the " *
            "spin one and which no entropy inequality would flag (both are honest von " *
            "Neumann entropies).  Use adjacent intervals — A = 1:2, B = 3:4 makes " *
            "A ∪ B = 1:4 contiguous — or a dense-ED hub, which traces out the spin " *
            "complement and takes any region.",
        ),
    )
end

"""
    _partial_trace_sites(ρ, sites, N; d=2) -> Matrix

Reduced density matrix of `ρ` (a `d^N × d^N` matrix) on `sites`, tracing out the
complement.  `sites` may be any subset — contiguous or not.

THE CONVENTION, which is normative for this atlas: **site 1 is the slowest-
varying tensor factor**, i.e. the leftmost in `reduce(kron, ops)`.  That is what
`_pauli_string` and `_spin1_string` build, so it is what the density matrices
handed to this function carry.  Under it, site `i` is tensor dimension `N+1-i`
of the ket group and `2N+1-i` of the bra group.

The convention is not readable off a `reshape`: column-major makes the FIRST
reshape slot the FASTEST-varying one, so `reshape(ρ, (dA, dB, dA, dB))` keeps
the LAST sites while `reshape(ρ, (dB, dA, dB, dA))` keeps the first.  The atlas
had one of each (#785), and neither produced a wrong number, because purity
(`S(A) = S(Aᶜ)`) and reflection symmetry make `S(first ℓ) = S(last ℓ)` on these
chains.  A region API removes both alibis: ask for `Region(2,3)` under the wrong
convention and you get `Region(N-2, N-1)` — and the entropy inequalities cannot
catch it, since a consistent relabelling of regions satisfies all of them.  So
this is settled by construction in the tests, on a product state where the
entropy is 0 everywhere and only the reduced state itself can decide.
"""
function _partial_trace_sites(
    ρ::AbstractMatrix, sites::AbstractVector{<:Integer}, N::Int; d::Int=2
)
    rest = setdiff(1:N, sites)
    isempty(rest) && return Matrix(ρ)
    T = reshape(Array(ρ), ntuple(_ -> d, 2N))
    ket(i) = N + 1 - i
    bra(i) = 2N + 1 - i
    # `reverse` so that within each group the LOWEST site number ends up
    # slowest-varying again, preserving the convention inside the block.
    perm = (
        reverse(ket.(sites))...,
        reverse(ket.(rest))...,
        reverse(bra.(sites))...,
        reverse(bra.(rest))...,
    )
    dA = d^length(sites)
    dB = d^length(rest)
    # `permutedims` moves the A dimensions to the FRONT, which in column-major
    # order makes them the FAST group -- the opposite of the unpermuted layout,
    # where site 1 is slowest overall and A is the SLOW group.  So the reshape
    # here is (dA, dB, dA, dB), not the (dB, dA, dB, dA) that is correct for a
    # helper reading the natural layout (`_s1_partial_trace_A`).  Within the A
    # composite the `reverse` above leaves site 1 slowest, matching `kron`.
    R = reshape(permutedims(T, perm), (dA, dB, dA, dB))
    ρA = zeros(eltype(ρ), dA, dA)
    @inbounds for a1 in 1:dA, a2 in 1:dA
        s = zero(eltype(ρ))
        for b in 1:dB
            s += R[a1, b, a2, b]
        end
        ρA[a1, a2] = s
    end
    return ρA
end

"""
    _reduced_from_pure(ψ, sites, N; d=2) -> Matrix

Reduced density matrix of the pure state `ψ` on `sites`, without forming the
`d^N × d^N` outer product.  Same convention as [`_partial_trace_sites`](@ref).
"""
function _reduced_from_pure(
    ψ::AbstractVector, sites::AbstractVector{<:Integer}, N::Int; d::Int=2
)
    rest = setdiff(1:N, sites)
    isempty(rest) && return ψ * ψ'
    T = reshape(Array(ψ), ntuple(_ -> d, N))
    perm = (reverse((N + 1) .- sites)..., reverse((N + 1) .- rest)...)
    dA = d^length(sites)
    dB = d^length(rest)
    # A was permuted to the front ⇒ it is the FIRST slot of the reshape
    M = reshape(permutedims(T, perm), (dA, dB))
    return M * M'
end

"""
    _majorana_indices(sites) -> Vector{Int}

The Majorana indices of `sites` under the convention that the pair
`(γ_{2i-1}, γ_{2i})` is local to spin site `i` (the σˣ-string Jordan–Wigner
convention used by `_majorana_ham`).  Sorted whenever `sites` is.
"""
function _majorana_indices(sites::AbstractVector{<:Integer})
    idx = Vector{Int}(undef, 2 * length(sites))
    @inbounds for (a, s) in enumerate(sites)
        idx[2a - 1] = 2s - 1
        idx[2a] = 2s
    end
    return idx
end
