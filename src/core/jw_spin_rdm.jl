# ─────────────────────────────────────────────────────────────────────────────
# core/jw_spin_rdm.jl — the SPIN reduced density matrix of a DISCONNECTED region
# of a free-fermion chain, from the Majorana covariance alone.
#
# `core/regions.jl` explains why the free-fermion routes refuse a multi-interval
# region: restricting the covariance to a site set gives the FERMIONIC entropy,
# and on a disconnected region that is a different number from the spin one,
# which no entropy inequality would flag.  This file removes the restriction
# rather than working around it.
#
# THE CONSTRUCTION.  Write `γ̃` for the Jordan-Wigner fermions built *within* the
# region and `γ` for the chain's own.  They differ by the string over the
# non-region sites to the left, so for a Majorana index in run `r`,
#
#     γ̃_a = P₀ · G₁ ⋯ G_{r-1} · γ_a ,      G_g = ∏_{m ∈ gap g} σˣ_m
#
# with `P₀` the prefix before the first run.  A Majorana monomial `γ̃_T` therefore
# equals `∏_g G_g^{n_g}` times `γ_T`, where `n_g` counts the `T`-indices lying in
# runs AFTER gap `g`; `P₀` appears `|T|` times and cancels, since only even `|T|`
# has non-zero expectation.  So
#
#     ⟨γ̃_T⟩ = ⟨ γ_{S} γ_T ⟩ ,   S = the Majorana indices of the gaps with n_g ODD,
#
# and the spin reduced density matrix is the monomial sum
#
#     ρ_s = 2^{-n} Σ_T (-1)^{|T|(|T|-1)/2} ⟨γ̃_T⟩ γ̃_T .
#
# WHY THIS NEEDS NO GAUSSIAN-OPERATOR MACHINERY.  `G_g = ∏_{m} σˣ_m` is itself a
# Majorana MONOMIAL — `σˣ_m = -i γ_{2m-1} γ_{2m}` in the σˣ-string convention
# `_majorana_ham` uses — so `⟨γ_S γ_T⟩` is one more Pfaffian of the SAME plain
# covariance, over an enlarged index set.  Every ingredient is `Pf(M[·])` with
# `M_ab = ⟨γ_a γ_b⟩ = -i Σ_ab` (verified against exact diagonalisation).
#
# That matters for more than tidiness.  The natural alternative — normalise by
# `⟨G⟩` and build a Gaussian operator with the string-inserted covariance — is
# singular exactly where it is most needed: `⟨G⟩` vanishes identically for an
# odd-length gap at half filling, by particle-hole symmetry, which is half of all
# configurations.  Nothing here divides by `⟨G⟩`, so those cases are ordinary.
#
# COST is `4^{|A|}` Pfaffians and a `2^{|A|}` accumulation each — exponential in
# the REGION, polynomial in the CHAIN.  That is the whole point: dense ED costs
# `2^N` and caps at `N = 12`, while this handles `N = 200` provided the region
# stays small.  It is not a replacement for the single-interval route, which is
# `O(ℓ³)` and should keep serving contiguous blocks.
# ─────────────────────────────────────────────────────────────────────────────

"""
    _pfaffian(A::AbstractMatrix) -> Number

Pfaffian of a skew-symmetric matrix by the Parlett-Reid algorithm, `O(m³)`.

Zero for odd dimension (the Pfaffian of an odd skew-symmetric matrix), one for
the empty matrix.  Pivoting is on magnitude, which is what keeps it usable when
a whole block of the covariance is near zero — the half-filling case that makes
the normalised route singular.
"""
function _pfaffian(A0::AbstractMatrix)
    n = size(A0, 1)
    n == size(A0, 2) || throw(ArgumentError("_pfaffian: matrix must be square"))
    isodd(n) && return zero(eltype(A0))
    n == 0 && return one(eltype(A0))
    A = Matrix(copy(A0))
    pf = one(eltype(A))
    for k in 1:2:(n - 1)
        # pivot the largest |A[k+1:n, k]| into row/col k+1
        kp = k + argmax(abs.(@view A[(k + 1):n, k]))
        if kp != k + 1
            A[[k + 1, kp], :] = A[[kp, k + 1], :]
            A[:, [k + 1, kp]] = A[:, [kp, k + 1]]
            pf = -pf
        end
        iszero(A[k, k + 1]) && return zero(eltype(A))
        pf *= A[k, k + 1]
        if k + 2 <= n
            τ = A[k, (k + 2):n] ./ A[k, k + 1]
            v = A[(k + 2):n, k + 1]
            A[(k + 2):n, (k + 2):n] .+= τ * transpose(v) .- v * transpose(τ)
        end
    end
    return pf
end

# A Majorana monomial in the Jordan-Wigner representation IS a Pauli string, so
# it is stored as (i-power, X-mask, Z-mask) and written into a dense matrix in
# O(2ⁿ) — one non-zero per column.  Building the monomials as dense matrices
# instead would cost O(8ⁿ) each and put the whole routine out of reach at |A| = 6.
#
# Convention (internal, and it does not escape): qubit j occupies bit n-j, and
# γ̃_{2j-1} = (∏_{m<j} X_m) Z_j,  γ̃_{2j} = (∏_{m<j} X_m) Y_j  with Y = i·X·Z.
# Only the SPECTRUM of the result is used, so this need not match the site
# ordering of `_partial_trace_sites`.
function _majorana_pauli(a::Int, n::Int)
    j = (a + 1) ÷ 2
    xmask = 0
    for m in 1:(j - 1)
        xmask |= 1 << (n - m)
    end
    bit = 1 << (n - j)
    return isodd(a) ? (0, xmask, bit) : (1, xmask | bit, bit)
end

# (p₁,x₁,z₁)·(p₂,x₂,z₂):  Z₁ past X₂ gives (-1)^{popcount(z₁ & x₂)}
function _pauli_mul(P, Q)
    p1, x1, z1 = P
    p2, x2, z2 = Q
    return (p1 + p2 + 2 * count_ones(z1 & x2), x1 ⊻ x2, z1 ⊻ z2)
end

# accumulate  coeff · (i^p ∏ X^x Z^z)  into the dense ρ
function _add_pauli!(ρ::AbstractMatrix, coeff, P, n::Int)
    p, xmask, zmask = P
    ip = (im)^mod(p, 4)
    c = coeff * ip
    @inbounds for b in 0:(2 ^ n - 1)
        s = isodd(count_ones(b & zmask)) ? -1 : 1
        ρ[(b ⊻ xmask) + 1, b + 1] += s * c
    end
    return ρ
end

"""
    MAX_JW_SPIN_REGION :: Int

Largest region size the [`spin_rdm_from_covariance`](@ref) monomial sum will
attempt.  The work is `4^{|A|}` Pfaffians, so this is a wall, not a preference:
`|A| = 8` runs in about a second and each further site costs roughly 8×.
"""
const MAX_JW_SPIN_REGION = 10

"""
    spin_rdm_from_covariance(Σ, sites, N) -> Matrix{ComplexF64}

The **spin** reduced density matrix on `sites` of an `N`-site free-fermion chain
whose Majorana covariance is `Σ` (the real antisymmetric matrix
`_majorana_covariance_gs` returns, related to the correlator by
`⟨γ_a γ_b⟩ = -i Σ_ab`).

`sites` may be **any** subset — this is the multi-interval case the free-fermion
routes otherwise refuse.  For a single contiguous interval the result is the
fermionic reduced state as well, and the two agree exactly; for a disconnected
one they differ, and this returns the spin answer.

Cost is exponential in `|sites|` and polynomial in `N`, which is the opposite
trade to dense ED — see the file header.
"""
function spin_rdm_from_covariance(
    Σ::AbstractMatrix, sites::AbstractVector{<:Integer}, N::Int
)
    A = sort!(collect(Int, sites))
    n = length(A)
    n ≥ 1 || throw(ArgumentError("spin_rdm_from_covariance: the region is empty."))
    (first(A) ≥ 1 && last(A) ≤ N) || throw(
        ArgumentError(
            "spin_rdm_from_covariance: region sites must lie in 1:N; got $A with N = $N.",
        ),
    )
    n ≤ MAX_JW_SPIN_REGION || throw(
        ArgumentError(
            "spin_rdm_from_covariance: |region| = $n exceeds MAX_JW_SPIN_REGION = " *
            "$(MAX_JW_SPIN_REGION); the monomial sum is 4^|region| Pfaffians. Use a " *
            "smaller region, or the single-interval route if the region is contiguous.",
        ),
    )
    size(Σ, 1) == 2N || throw(
        ArgumentError(
            "spin_rdm_from_covariance: Σ must be 2N × 2N; got $(size(Σ,1)) with N = $N."
        ),
    )
    M = (-im) .* Matrix{ComplexF64}(Σ)
    for a in axes(M, 1)
        M[a, a] = 0                       # the Pfaffian sees only the skew part
    end

    # runs of consecutive sites, and the gaps between them
    runs = Vector{UnitRange{Int}}()
    p = A[1]
    for k in 2:n
        if A[k] != A[k - 1] + 1
            push!(runs, p:A[k - 1])
            p = A[k]
        end
    end
    push!(runs, p:A[end])
    gaps = [(last(runs[g]) + 1):(first(runs[g + 1]) - 1) for g in 1:(length(runs) - 1)]
    gap_majoranas = [vcat([[2m - 1, 2m] for m in g]...) for g in gaps]
    # run index of each local Majorana, so n_g can be counted per monomial
    run_of_local = Vector{Int}(undef, 2n)
    for (loc, site) in enumerate(A)
        r = findfirst(rr -> site in rr, runs)
        run_of_local[2loc - 1] = r
        run_of_local[2loc] = r
    end
    global_majorana = vcat([[2s - 1, 2s] for s in A]...)

    ρ = zeros(ComplexF64, 2^n, 2^n)
    monomial = Vector{Tuple{Int,Int,Int}}(undef, 2n + 1)   # prefix products
    T = Int[]
    invtwo = 1 / 2^n

    # depth-first over subsets in increasing order, so each extension is one
    # Pauli multiply rather than a rebuild
    monomial[1] = (0, 0, 0)
    function walk(start::Int)
        k = length(T)
        if iseven(k)
            ngs = [count(t -> run_of_local[t] > g, T) for g in 1:length(gaps)]
            S = Int[]
            for g in 1:length(gaps)
                isodd(ngs[g]) && append!(S, gap_majoranas[g])
            end
            Tg = [global_majorana[t] for t in T]
            U = vcat(S, Tg)
            val = if isempty(U)
                one(ComplexF64)
            else
                # No reordering sign: sorting `U = [S; T_global]` is always an EVEN
                # permutation, so the factor is identically +1. Each gap block
                # contributes 2·|gap_g| indices and each moves past the same count
                # m_g of T-indices, i.e. 2·|gap_g|·m_g transpositions. MEASURED over
                # 1311 configurations, and mutation-checked: multiplying by the sign
                # changed no test, which is what flagged it as dead. The parity claim
                # is asserted in test_jw_spin_rdm.jl rather than left implicit here.
                (-im)^(length(S) ÷ 2) * _pfaffian(M[sort(U), sort(U)])
            end
            if val != 0
                ε = iseven(k * (k - 1) ÷ 2) ? 1 : -1
                _add_pauli!(ρ, ε * val * invtwo, monomial[k + 1], n)
            end
        end
        for a in start:(2n)
            push!(T, a)
            monomial[length(T) + 1] = _pauli_mul(monomial[length(T)], _majorana_pauli(a, n))
            walk(a + 1)
            pop!(T)
        end
        return nothing
    end
    walk(1)
    return (ρ + ρ') ./ 2      # Hermitian by construction; symmetrise off round-off
end
