using QAtlas, LinearAlgebra, Test

# Region as a fetch argument for the TFIM free-fermion entanglement entropies
# (#780).  Two facts decide this design, and both are measured here rather than
# assumed:
#
#   (a) a contiguous interval ANYWHERE in the chain — not just `1:ℓ` — has
#       S_fermion = S_spin.  This is what lets the region entropy inequalities
#       be instantiated on ADJACENT blocks A = 1:2, B = 3:4, C = 5:6, where
#       every union `region_report` needs (A∪B, B∪C, A∪B∪C) is again a single
#       interval.
#
#   (b) a MULTI-interval region does NOT.  The free-fermion submatrix answers
#       a different question there, and the gap is large.  The implementation
#       throws rather than returning it; the test below records how wrong the
#       silent answer would have been, so that "just allow multi-interval"
#       cannot look like a harmless relaxation later.
#
# `_build_tfim_dense` comes from test/util/tfim_dense_ed.jl (loaded by
# runtests.jl) and uses site 1 as the OUTERMOST kron factor.

# Reduced density matrix of the SPIN state on an arbitrary site subset.
# (test/util/spinhalf_ed.jl only traces out the complement of a leading block.)
function _rdm_spin_sites(ψ::AbstractVector, N::Int, sites::AbstractVector{Int})
    rest = setdiff(1:N, sites)
    T = reshape(ψ, ntuple(_ -> 2, N))
    # site i is dimension N+1-i under the site-1-outermost kron convention
    perm = (reverse(N .+ 1 .- sites)..., reverse(N .+ 1 .- rest)...)
    M = reshape(permutedims(T, perm), 2^length(sites), 2^length(rest))
    return M * M'
end

function _vn_entropy(ρ::AbstractMatrix)
    return -sum(λ -> λ > 1e-14 ? λ * log(λ) : 0.0, eigvals(Hermitian(ρ)))
end

# The free-fermion answer on an ARBITRARY site set — i.e. exactly what the
# implementation would return if `_require_single_interval` were removed.
function _free_fermion_entropy(m::TFIM, N::Int, sites::Vector{Int}; beta=Inf)
    Σ = QAtlas._majorana_thermal_covariance(QAtlas._majorana_ham(N, m.J, m.h), beta)
    idx = QAtlas._majorana_indices(sites)
    λ = eigvals(Hermitian(im .* Σ[idx, idx]))
    L = length(sites)
    return sum(QAtlas._peschel_mode_entropy(λ[k]) for k in (L + 1):(2L))
end

@testset "TFIM entanglement — `region` reproduces `ℓ`" begin
    N = 12
    for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0))
        m = TFIM(; J=J, h=h)
        for ℓ in 1:(N - 1)
            S_ℓ = QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); ℓ=ℓ)
            S_r = QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); region=Region(1:ℓ...))
            @test S_r ≈ S_ℓ atol = 1e-12
        end
        # Rényi carries the same sugar
        for α in (0.5, 2.0)
            S_ℓ = QAtlas.fetch(m, RenyiEntropy(α), OBC(N); ℓ=4)
            S_r = QAtlas.fetch(m, RenyiEntropy(α), OBC(N); region=Region(1, 2, 3, 4))
            @test S_r ≈ S_ℓ atol = 1e-12
        end
    end
end

@testset "TFIM entanglement — (a) a contiguous block ANYWHERE matches full ED" begin
    # This is the fact the adjacent-block instantiation of the region
    # inequalities rests on.  If it ever breaks, `region_report` over a
    # free-fermion hub silently starts comparing fermionic to spin entropies.
    N = 8
    for (J, h) in ((1.0, 0.5), (1.0, 1.0))
        ψ = eigen(Symmetric(Matrix(real(_build_tfim_dense(N, J, h))))).vectors[:, 1]
        m = TFIM(; J=J, h=h)
        for sites in ([1], [3], [1, 2], [3, 4], [5, 6], [2, 3, 4], [3, 4, 5, 6], [5, 6, 7])
            S_ed = _vn_entropy(_rdm_spin_sites(ψ, N, sites))
            S_ff = QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); region=Region(sites...))
            @test S_ff ≈ S_ed atol = 1e-10
        end
    end
end

@testset "TFIM entanglement — OBC reflection symmetry of a shifted block" begin
    # An OBC chain is symmetric under i → N+1-i, so a block and its mirror
    # carry the same entropy.  Independent of the ED reference above: this one
    # is a structural law, not a second computation of the same number.
    N = 10
    m = TFIM(; J=1.0, h=1.0)
    for (a, b) in ((1, 2), (2, 4), (3, 5))
        left = QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); region=Region((a:b)...))
        right = QAtlas.fetch(
            m, VonNeumannEntropy(), OBC(N); region=Region(((N + 1 - b):(N + 1 - a))...)
        )
        @test left ≈ right atol = 1e-12
    end
end

@testset "TFIM entanglement — (b) multi-interval regions answer the SPIN question" begin
    N = 8
    m = TFIM(; J=1.0, h=0.5)
    ψ = eigen(Symmetric(Matrix(real(_build_tfim_dense(N, 1.0, 0.5))))).vectors[:, 1]

    # This used to assert a THROW.  #832 reinstates the Jordan-Wigner string
    # explicitly, so the route now answers — and the assertion that matters is
    # that it answers the SPIN question, checked against this file's own dense-ED
    # helper rather than against anything the new code shares.
    for sites in ([1, 3], [1, 4], [2, 4], [1, 2, 4], [1, 2, 5, 6])
        got = QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); region=Region(sites...))
        S_spin = _vn_entropy(_rdm_spin_sites(ψ, N, sites))
        S_ferm = _free_fermion_entropy(m, N, sites)
        @test got ≈ S_spin atol = 1e-10

        # and the free-fermion submatrix, which is what a naive route would have
        # returned, differs by ~0.5 nats — not a rounding error.  Both are honest
        # von Neumann entropies, so no entropy inequality downstream would have
        # caught the substitution; only this comparison does.
        @test abs(S_spin - S_ferm) > 0.4
        @test !isapprox(got, S_ferm; atol=0.1)
        # …and asking for the fermionic one explicitly gives that other number
        @test QAtlas.fetch(
            m, FermionicEntanglementEntropy(), OBC(N); region=Region(sites...)
        ) ≈ S_ferm atol = 1e-10
    end

    # the Rényi method takes the same split, so the SHAPE of a region never
    # decides which quantities exist.  Checked against Tr ρ^α of this file's own
    # dense-ED reduced state, not against the von Neumann route.
    for sites in ([1, 3], [1, 2, 5, 6]), α in (0.5, 2.0)
        λ = eigvals(Symmetric(Matrix(real(_rdm_spin_sites(ψ, N, sites)))))
        want = log(sum(x -> x > 1e-13 ? x^α : 0.0, real(λ))) / (1 - α)
        @test QAtlas.fetch(m, RenyiEntropy(α), OBC(N); region=Region(sites...)) ≈ want atol =
            1e-9
    end
end

@testset "TFIM entanglement — region argument validation" begin
    N = 8
    m = TFIM(; J=1.0, h=1.0)
    # neither given
    @test_throws ArgumentError QAtlas.fetch(m, VonNeumannEntropy(), OBC(N))
    # both given — `ℓ = k` IS `region = Region(1:k)`, so asking for both is a bug
    @test_throws ArgumentError QAtlas.fetch(
        m, VonNeumannEntropy(), OBC(N); region=Region(1, 2), ℓ=2
    )
    # outside the chain
    @test_throws ArgumentError QAtlas.fetch(
        m, VonNeumannEntropy(), OBC(N); region=Region(7, 8, 9)
    )
    @test_throws ArgumentError QAtlas.fetch(
        m, VonNeumannEntropy(), OBC(N); region=Region(0)
    )
    # the whole system: S of a pure state is 0 by construction, never evidence
    @test_throws ArgumentError QAtlas.fetch(
        m, VonNeumannEntropy(), OBC(N); region=Region((1:N)...)
    )
    @test_throws ArgumentError QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); ℓ=N)
    # a bare Integer in the `region` slot is ambiguous (site 3 vs the block 1:3)
    # and is refused rather than guessed
    @test_throws ArgumentError QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); region=3)
    # ... while both disambiguated spellings work and mean different things
    @test QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); region=Region(3)) !=
        QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); ℓ=3)
end
