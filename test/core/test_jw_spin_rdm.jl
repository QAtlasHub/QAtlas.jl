# The Jordan-Wigner string reinstated: the SPIN reduced density matrix of a
# DISCONNECTED region from the Majorana covariance (core/jw_spin_rdm.jl, #832).
#
# The oracle is dense ED, which is exactly the right use of an `:exponential`
# row — the instrument a cheaper route is checked against.  The point of the new
# route is that it is exponential in the REGION and polynomial in the CHAIN, so
# the assertions that matter are (a) it agrees with ED where ED can run, and
# (b) it still runs where ED cannot.

using Test
using LinearAlgebra
using QAtlas
using QAtlas: TFIM, OBC, fetch, VonNeumannEntropy, FermionicEntanglementEntropy
using QAtlas: Region, RenyiEntropy
using QAtlas: spin_rdm_from_covariance, MAX_JW_SPIN_REGION
using QAtlas: _majorana_ham, _majorana_covariance_gs, _pfaffian

@testset "_pfaffian against the recursive definition" begin
    # Parlett-Reid is O(m³) and not obviously right; the recursive expansion IS
    # the definition, so it is the oracle here rather than a second fast method.
    function pf_recursive(A)
        n = size(A, 1)
        n == 0 && return one(eltype(A))
        isodd(n) && return zero(eltype(A))
        n == 2 && return A[1, 2]
        s = zero(eltype(A))
        for j in 2:n
            rest = [k for k in 1:n if k != 1 && k != j]
            s += (-1)^j * A[1, j] * pf_recursive(A[rest, rest])
        end
        return s
    end
    for m in (0, 2, 4, 6, 8), trial in 1:3
        B = randn(m, m)
        A = B - B'
        @test _pfaffian(A) ≈ pf_recursive(A) rtol = 1e-9
        # Pf(A)² = det(A) — an independent structural law, not a second algorithm
        m > 0 && @test _pfaffian(A)^2 ≈ det(A) rtol = 1e-9
    end
    # odd dimension is zero by definition, not by round-off
    @test _pfaffian(zeros(3, 3)) == 0
    @test _pfaffian(zeros(0, 0)) == 1
    # complex entries: the string-inserted blocks are imaginary
    for trial in 1:3
        B = randn(ComplexF64, 6, 6)
        A = B - transpose(B)
        @test _pfaffian(A) ≈ pf_recursive(A) rtol = 1e-9
    end
end

@testset "the string/monomial reordering is always an EVEN permutation" begin
    # `spin_rdm_from_covariance` sorts `[gap-Majoranas; region-Majoranas]` before
    # taking the Pfaffian and applies NO sign for it.  That is a claim, not an
    # omission: each gap block contributes 2·|gap| indices which all move past the
    # same count of region indices, so the transposition count is even.  A
    # mutation that inserted the sign changed no other test — this is what makes
    # the claim checked rather than merely true.
    function permsign(v)
        w = collect(v)
        s = 1
        for i in 1:length(w), j in 1:(length(w) - 1)
            if w[j] > w[j + 1]
                w[j], w[j + 1] = w[j + 1], w[j]
                s = -s
            end
        end
        return s
    end
    nchecked = 0
    for A in (
        [1, 3],
        [1, 2, 5, 6],
        [1, 2, 4, 5],
        [2, 3, 6, 7],
        [1, 3, 5],
        [2, 4, 6, 8],
        [2, 3, 4, 8],
        [1, 4, 7, 10],
        [1, 2, 3, 7, 8],
    )
        n = length(A)
        runs = UnitRange{Int}[]
        p = A[1]
        for k in 2:n
            A[k] == A[k - 1] + 1 || (push!(runs, p:A[k - 1]); p=A[k])
        end
        push!(runs, p:A[end])
        gaps = [(last(runs[g]) + 1):(first(runs[g + 1]) - 1) for g in 1:(length(runs) - 1)]
        gapm = [vcat([[2m - 1, 2m] for m in g]...) for g in gaps]
        runof = zeros(Int, 2n)
        for (loc, site) in enumerate(A)
            r = findfirst(rr -> site in rr, runs)
            runof[2loc - 1] = r
            runof[2loc] = r
        end
        gm = vcat([[2s - 1, 2s] for s in A]...)
        for mask in 0:(2 ^ (2n) - 1)
            iseven(count_ones(mask)) || continue
            T = [t for t in 1:(2n) if (mask >> (t - 1)) & 1 == 1]
            S = Int[]
            for g in 1:length(gaps)
                isodd(count(t -> runof[t] > g, T)) && append!(S, gapm[g])
            end
            U = vcat(S, [gm[t] for t in T])
            isempty(U) && continue
            nchecked += 1
            @test permsign(U) == 1
        end
    end
    @test nchecked > 1000          # the sweep really ran, rather than matching nothing
end

# dense-ED oracle: spin reduced density matrix of the OBC TFIM ground state.
# The ground state is CACHED per (N, J, h) — rebuilding and diagonalising it per
# region cost 3m43s for seven regions, and the oracle is the same state each time.
const _ED_CACHE = Dict{Tuple{Int,Float64,Float64},Vector{ComplexF64}}()
function _ed_ground_state(N, J, h)
    return get!(_ED_CACHE, (N, J, h)) do
        sx, sz = [0.0 1; 1 0], [1.0 0; 0 -1]
        function o(m, i)
            return kron(
                Matrix(1.0I, 2^(i - 1), 2^(i - 1)),
                ComplexF64.(m),
                Matrix(1.0I, 2^(N - i), 2^(N - i)),
            )
        end
        H = zeros(ComplexF64, 2^N, 2^N)
        for i in 1:(N - 1)
            H -= J * o(sz, i) * o(sz, i + 1)
        end
        for i in 1:N
            H -= h * o(sx, i)
        end
        return eigen(Hermitian(H)).vectors[:, 1]
    end
end
function _ed_spin_entropy(N, J, h, sites)
    ψ = _ed_ground_state(N, J, h)
    rest = setdiff(1:N, sites)
    perm = (reverse((N + 1) .- sites)..., reverse((N + 1) .- rest)...)
    M = reshape(
        permutedims(reshape(Array(ψ), ntuple(_ -> 2, N)), perm),
        (2^length(sites), 2^length(rest)),
    )
    ρ = M * M'
    return -sum(λ -> λ > 1e-13 ? λ * log(λ) : 0.0, real(eigvals(Hermitian(ρ))))
end

@testset "disconnected SPIN entropy matches dense ED" begin
    J, h = 1.0, 1.0
    # Two chain lengths, one EVEN and one ODD.  The agreement is exact at any N,
    # so this is not a convergence check — it is the "vary an axis that should
    # not matter" one, and odd N in particular has no half-filling symmetry to
    # accidentally rescue a sign.
    for N in (10, 11)
        m = TFIM(J, h)
        bc = OBC(N)
        for sites in (
            [1, 3],                 # single sites, gap of 1 — the ⟨G⟩ = 0 case
            [1, 2, 5, 6],           # two blocks, EVEN gap
            [1, 2, 4, 5],           # two blocks, ODD gap
            [2, 3, 6, 7],           # neither block touches the boundary
            [2, 3, 4, 8],           # unequal blocks, long gap
            [1, 3, 5],              # THREE runs — two independent strings
            [2, 4, 6, 8],           # four runs
        )
            last(sites) ≤ N - 1 || continue
            got = fetch(m, VonNeumannEntropy(), bc; region=Region(sites...))
            @test got ≈ _ed_spin_entropy(N, J, h, sites) atol = 1e-10
        end
    end
    # off criticality too — the construction must not lean on h = J
    for (J, h) in ((1.0, 0.5), (1.0, 2.0))
        got = fetch(TFIM(J, h), VonNeumannEntropy(), OBC(10); region=Region(1, 2, 5, 6))
        @test got ≈ _ed_spin_entropy(10, J, h, [1, 2, 5, 6]) atol = 1e-10
    end
end

@testset "the two routes agree on a contiguous region and differ off it" begin
    N = 12
    m = TFIM(1.0, 1.0)
    bc = OBC(N)
    Σ = _majorana_covariance_gs(_majorana_ham(N, 1.0, 1.0))
    # contiguous: the Jordan-Wigner string factorises, so the reconstructed spin
    # state must reproduce the O(ℓ³) covariance route EXACTLY.  This is the
    # assertion that would catch a sign or ordering error in the monomial sum,
    # because it compares against a route that shares none of its code.
    for sites in ([1, 2], [3, 4, 5], [4, 5, 6, 7], [9, 10])
        ρ = spin_rdm_from_covariance(Σ, sites, N)
        S = -sum(λ -> λ > 1e-13 ? λ * log(λ) : 0.0, real(eigvals(Hermitian(ρ))))
        @test S ≈ fetch(m, VonNeumannEntropy(), bc; region=Region(sites...)) atol = 1e-10
        @test tr(ρ) ≈ 1 atol = 1e-10
        @test minimum(real(eigvals(Hermitian(ρ)))) > -1e-10   # a state, not just Hermitian
    end
    # disconnected: spin and fermionic are DIFFERENT, and by a margin no
    # inequality would notice
    for sites in ([1, 3], [1, 2, 5, 6], [2, 3, 6, 7])
        r = Region(sites...)
        Ss = fetch(m, VonNeumannEntropy(), bc; region=r)
        Sf = fetch(m, FermionicEntanglementEntropy(), bc; region=r)
        @test !isapprox(Ss, Sf; atol=1e-3)
        @test abs(Ss - Sf) > 0.05
        @test Ss > 0 && Sf > 0
    end
end

@testset "Renyi follows the same split, and matches ED" begin
    # Leaving Rényi guarded while von Neumann answered would let the SHAPE of the
    # region decide which quantities exist. Checked against the same dense-ED
    # oracle, on the reduced state rather than through any shared code path.
    N, J, h = 10, 1.0, 1.0
    m = TFIM(J, h)
    ψ = _ed_ground_state(N, J, h)
    for sites in ([1, 3], [1, 2, 5, 6], [2, 4, 6, 8]), α in (0.5, 2.0, 3.0)
        rest = setdiff(1:N, sites)
        perm = (reverse((N + 1) .- sites)..., reverse((N + 1) .- rest)...)
        M = reshape(
            permutedims(reshape(Array(ψ), ntuple(_ -> 2, N)), perm),
            (2^length(sites), 2^length(rest)),
        )
        λ = real(eigvals(Hermitian(M * M')))
        want = log(sum(x -> x > 1e-13 ? x^α : 0.0, λ)) / (1 - α)
        got = fetch(m, RenyiEntropy(α), OBC(N); region=Region(sites...))
        @test got ≈ want atol = 1e-9
    end
    # contiguous still takes the O(ℓ³) route and still agrees with it
    @test fetch(m, RenyiEntropy(2.0), OBC(N); region=Region(1, 2, 3)) ≈
        fetch(m, RenyiEntropy(2.0), OBC(N); ℓ=3) atol = 1e-12
end

@testset "it runs where dense ED cannot" begin
    # The whole trade: exponential in the REGION, polynomial in the CHAIN.  Dense
    # ED caps at N = 12; these are N = 64 and N = 160.
    m = TFIM(1.0, 1.0)
    for N in (64, 160)
        r = Region(10, 11, 30, 31)
        S = fetch(m, VonNeumannEntropy(), OBC(N); region=r)
        @test isfinite(S) && 0 < S < 4 * log(2)
        # purity: S(A) of a disconnected A is bounded by |A| log 2
        @test S < 4 * log(2) + 1e-9
        # subadditivity against its own parts, which the region generator relies on
        S1 = fetch(m, VonNeumannEntropy(), OBC(N); region=Region(10, 11))
        S2 = fetch(m, VonNeumannEntropy(), OBC(N); region=Region(30, 31))
        @test S1 + S2 - S > -1e-10                       # I(A:B) ≥ 0
        @test S ≥ abs(S1 - S2) - 1e-10                   # Araki–Lieb
    end
    # far-separated blocks decorrelate: the mutual information must DECAY, which
    # a route that silently dropped the string would not reproduce (it would sit
    # at the fermionic value instead)
    N = 160
    mi = Float64[]
    for d in (4, 8, 16, 32, 64)
        S1 = fetch(m, VonNeumannEntropy(), OBC(N); region=Region(20, 21))
        S2 = fetch(m, VonNeumannEntropy(), OBC(N); region=Region(20 + d, 21 + d))
        S12 = fetch(m, VonNeumannEntropy(), OBC(N); region=Region(20, 21, 20 + d, 21 + d))
        push!(mi, S1 + S2 - S12)
    end
    @test all(>(−1e-10), mi)
    @test issorted(mi; rev=true)                         # monotone decay in separation
    @test mi[end] < mi[1] / 10
end

@testset "refusals rather than a hang" begin
    N = 40
    Σ = _majorana_covariance_gs(_majorana_ham(N, 1.0, 1.0))
    big = collect(1:(MAX_JW_SPIN_REGION + 1))
    @test_throws ArgumentError spin_rdm_from_covariance(Σ, big, N)
    @test_throws ArgumentError spin_rdm_from_covariance(Σ, Int[], N)
    @test_throws ArgumentError spin_rdm_from_covariance(Σ, [0, 2], N)
    @test_throws ArgumentError spin_rdm_from_covariance(Σ, [1, N + 1], N)
    @test_throws ArgumentError spin_rdm_from_covariance(Σ, [1, 3], N + 1)   # Σ size
end
