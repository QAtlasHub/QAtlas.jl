using QAtlas, LinearAlgebra, Test

# Regions on the dense-ED entanglement hubs (#780 step 2), and the site-index
# convention they turn on (#785).
#
# Unlike the free-fermion hubs, these trace out the SPIN complement directly, so
# there is no Jordan-Wigner restriction: any subset of sites is fair game,
# contiguous or not.  That is what makes genuinely disjoint regions -- and so
# subadditivity and Araki-Lieb with real content -- reachable at all.
#
# THE CONVENTION CANNOT BE READ OFF A `reshape`.  Column-major makes the first
# reshape slot the FASTEST-varying one, so `(dA, dB, dA, dB)` keeps the LAST
# sites and `(dB, dA, dB, dA)` keeps the first.  The atlas held one of each:
# `_s1_partial_trace_A` kept the first `ℓ` (matching its docstring and
# `_pauli_string`), `_xxz1d_partial_trace_B` kept the last `ℓ` while claiming the
# first.  Neither ever produced a wrong number, because purity (S(A) = S(Aᶜ)) and
# the chains' reflection symmetry make S(first ℓ) = S(last ℓ).  A region argument
# removes both alibis, so the convention is pinned here BY CONSTRUCTION: on a
# product state the entropy is 0 on every region and only the reduced state
# itself can say which site came back.

@testset "dense-ED partial trace — the site convention, by construction" begin
    # spin-1/2: site k carries cos θ_k|0⟩ + sin θ_k|1⟩ with distinct θ_k, so
    # ρ_A[1,1] identifies the site unambiguously.
    N = 5
    θ = [0.1, 0.7, 1.2, 1.5, 0.4]
    kets = [[cos(t), sin(t)] for t in θ]
    ψ = reduce(kron, kets)      # same kron order as `_pauli_string`
    for k in 1:N
        ρk = QAtlas._partial_trace_sites(ψ * ψ', [k], N; d=2)
        @test real(ρk[1, 1]) ≈ cos(θ[k])^2 atol = 1e-12
    end
    # a NON-contiguous pair comes back as the product of exactly those two sites
    ρ13 = QAtlas._partial_trace_sites(ψ * ψ', [1, 3], N; d=2)
    want = kron(kets[1], kets[3]) * kron(kets[1], kets[3])'
    @test ρ13 ≈ want atol = 1e-12

    # the pure-state shortcut must agree with the full outer product
    for sites in ([2], [2, 4], [1, 3, 5])
        @test QAtlas._reduced_from_pure(ψ, sites, N; d=2) ≈
            QAtlas._partial_trace_sites(ψ * ψ', sites, N; d=2) atol = 1e-12
    end
end

@testset "dense-ED partial trace — agrees with the independent S1 implementation" begin
    # `_s1_partial_trace_A` reads the NATURAL layout (no permutedims) and was
    # already on the normative convention, so it is an independent oracle for
    # the shared helper's arithmetic -- a different code path, not a restatement.
    N = 3
    φ = [0.2, 0.9, 1.4]
    k3 = [normalize([cos(t), sin(t) * 0.6, sin(t) * 0.8]) for t in φ]
    ψ3 = normalize(reduce(kron, k3) + 0.7 * reduce(kron, reverse(k3)))  # entangled
    ρ3 = ψ3 * ψ3'
    for ℓ in 1:(N - 1)
        @test QAtlas._partial_trace_sites(ρ3, collect(1:ℓ), N; d=3) ≈
            QAtlas._s1_partial_trace_A(ρ3, ℓ, N) atol = 1e-12
    end
end

@testset "dense-ED partial trace — purity S(A) = S(Aᶜ), including disjoint A" begin
    # A structural law, independent of any reference implementation: for a PURE
    # state a region and its complement carry equal entropy.  It fails loudly
    # under a scrambled site convention, which is how the first draft of the
    # shared helper was caught.
    N = 5
    vn(ρ) = -sum(λ -> λ > 1e-14 ? λ * log(λ) : 0.0, real.(eigvals(Hermitian((ρ + ρ') / 2))))
    m = XXZ1D(; J=1.0, Δ=0.5)
    H = QAtlas._xxz1d_hamiltonian_matrix(m, N)
    ψ = eigen(Hermitian(H)).vectors[:, 1]
    for sites in ([1], [1, 2], [1, 3], [2, 4], [1, 3, 5], [2, 3])
        comp = setdiff(1:N, sites)
        @test vn(QAtlas._reduced_from_pure(ψ, sites, N; d=2)) ≈
            vn(QAtlas._reduced_from_pure(ψ, comp, N; d=2)) atol = 1e-10
    end
end

@testset "XXZ1D / S1Heisenberg1D — `region` reproduces `ℓ`" begin
    N = 6
    for m in (XXZ1D(; J=1.0, Δ=0.5), XXZ1D(; J=1.0, Δ=1.0))
        for ℓ in 1:(N - 1), beta in (Inf, 1.0)
            a = QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); ℓ=ℓ, beta=beta)
            b = QAtlas.fetch(
                m, VonNeumannEntropy(), OBC(N); region=Region(1:ℓ...), beta=beta
            )
            @test a ≈ b atol = 1e-12
        end
        @test QAtlas.fetch(m, RenyiEntropy(2.0), OBC(N); ℓ=3) ≈
            QAtlas.fetch(m, RenyiEntropy(2.0), OBC(N); region=Region(1, 2, 3)) atol = 1e-12
    end
    N1 = 4
    m1 = S1Heisenberg1D(; J=1.0)
    for ℓ in 1:(N1 - 1)
        @test QAtlas.fetch(m1, VonNeumannEntropy(), OBC(N1); ℓ=ℓ) ≈
            QAtlas.fetch(m1, VonNeumannEntropy(), OBC(N1); region=Region(1:ℓ...)) atol =
            1e-12
    end
    # Heisenberg1D forwards the pair through to XXZ1D
    mh = Heisenberg1D()
    @test QAtlas.fetch(mh, VonNeumannEntropy(), OBC(N); ℓ=2) ≈
        QAtlas.fetch(mh, VonNeumannEntropy(), OBC(N); region=Region(1, 2)) atol = 1e-12
end

@testset "XXZ1D — non-contiguous regions are ALLOWED here, and have content" begin
    # The free-fermion hubs refuse these (Jordan-Wigner); a dense-ED hub must
    # not, since it traces out the spin complement directly.
    N = 6
    m = XXZ1D(; J=1.0, Δ=0.5)
    S(sites; beta=Inf) =
        QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); region=Region(sites...), beta=beta)

    for sites in ([1, 3], [2, 4], [1, 4, 6], [1, 2, 5])
        s = S(sites)
        @test isfinite(s)
        @test s > 0                       # a disjoint region of a correlated GS is mixed
        @test s ≤ length(sites) * log(2) + 1e-10   # maximal-mixing bound
    end

    # Subadditivity and Araki-Lieb on GENUINELY disjoint, non-adjacent regions --
    # the instances the free-fermion route cannot supply at all.  These are the
    # relations #780 exists to reach.
    for (A, B) in (([1], [3]), ([1], [4]), ([1, 2], [4, 5]), ([2], [5]))
        S_A, S_B, S_AB = S(A), S(B), S(vcat(A, B))
        @test S_A + S_B - S_AB ≥ -1e-10              # I(A:B) ≥ 0
        @test S_AB - abs(S_A - S_B) ≥ -1e-10         # Araki-Lieb
    end

    # ... and at finite temperature, where the state is genuinely mixed and
    # purity no longer relates a region to its complement.
    for (A, B) in (([1], [3]), ([1, 2], [4, 5]))
        S_A, S_B, S_AB = S(A; beta=1.5), S(B; beta=1.5), S(vcat(A, B); beta=1.5)
        @test S_A + S_B - S_AB ≥ -1e-10
        @test S_AB - abs(S_A - S_B) ≥ -1e-10
    end
end

@testset "dense-ED regions — argument validation matches the free-fermion hubs" begin
    N = 5
    m = XXZ1D(; J=1.0, Δ=0.5)
    @test_throws ArgumentError QAtlas.fetch(m, VonNeumannEntropy(), OBC(N))
    @test_throws ArgumentError QAtlas.fetch(
        m, VonNeumannEntropy(), OBC(N); region=Region(1, 2), ℓ=2
    )
    @test_throws ArgumentError QAtlas.fetch(
        m, VonNeumannEntropy(), OBC(N); region=Region(4, 5, 6)
    )
    @test_throws ArgumentError QAtlas.fetch(
        m, VonNeumannEntropy(), OBC(N); region=Region((1:N)...)
    )
    @test_throws ArgumentError QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); region=3)
end
