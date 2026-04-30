using QAtlas, Lattice2D, LinearAlgebra, Test

# `_build_tfim_dense` is loaded once from `runtests.jl` via
# `test/util/tfim_dense_ed.jl`, so it is already in scope when this file
# is included by the test harness.

@testset "TFIM RenyiEntropy(α) OBC via correlation matrix" begin
    @testset "α → 1 limit matches Von Neumann" begin
        # α = 1.0001 (close to 1 but valid; the constructor rejects α = 1).
        # The Taylor expansion of the per-mode Rényi entropy around α = 1
        # gives s_α(ν) = s_vN(ν) + (1 - α)/2 · Var_p(log p) + O((1-α)²),
        # so the relative error scales linearly in |1 - α| = 1e-4.
        N, ℓ, h = 12, 6, 0.7
        model = TFIM(; J=1.0, h=h)
        S_vn = QAtlas.fetch(model, VonNeumannEntropy(), OBC(N); ℓ=ℓ, beta=Inf)
        S_renyi = QAtlas.fetch(model, RenyiEntropy(1.0001), OBC(N); ℓ=ℓ, beta=Inf)
        @test S_renyi ≈ S_vn rtol = 1e-3
    end

    @testset "α = 2 (collision entropy)" begin
        # Rényi is monotonically decreasing in α, so for α > 1
        # we must have 0 ≤ S_α ≤ S_vN.
        for h in (0.5, 1.0, 1.5), ℓ in (3, 6), N in (12,)
            model = TFIM(; J=1.0, h=h)
            S_vn = QAtlas.fetch(model, VonNeumannEntropy(), OBC(N); ℓ=ℓ)
            S_2 = QAtlas.fetch(model, RenyiEntropy(2.0), OBC(N); ℓ=ℓ)
            @test 0 ≤ S_2
            @test S_2 ≤ S_vn + 1e-10
        end
    end

    @testset "α = 0.5 (≥ S_vn)" begin
        # For α < 1 the inequality flips: S_α ≥ S_vN.
        for h in (0.5, 1.0, 1.5), ℓ in (3, 6), N in (12,)
            model = TFIM(; J=1.0, h=h)
            S_vn = QAtlas.fetch(model, VonNeumannEntropy(), OBC(N); ℓ=ℓ)
            S_h = QAtlas.fetch(model, RenyiEntropy(0.5), OBC(N); ℓ=ℓ)
            @test S_h + 1e-10 ≥ S_vn
        end
    end

    @testset "Compare against full ED at small N" begin
        # Build the GS, partial-trace, get RDM eigenvalues, compute
        # Rényi directly.  Tolerance is α-dependent: at α < 1 the sum
        # `Σ p^α` is dominated by *small* Schmidt eigenvalues, which the
        # dense `eigvals(Hermitian(ρA))` resolves only down to its
        # round-off floor `~ eps · max(p) ≈ 1e-16`.  Filtering those at
        # `1e-15` drops mass of order `2^ℓ · (1e-16)^α`, which is
        # `~1e-7` for α=0.5 and only `~1e-32` for α=2.  The Peschel
        # mode-product formula is exact, so the residual is the ED
        # reference's precision limit.
        for h in (0.5, 1.0, 1.5)
            N = 8
            ℓ = 4
            H = _build_tfim_dense(N, 1.0, h)
            E, V = eigen(H)
            ψ = V[:, 1]
            dA = 2^ℓ
            dB = 2^(N - ℓ)
            Ψ = reshape(ψ, (dA, dB))
            ρA = Ψ * Ψ'
            evals = real.(eigvals(Hermitian(ρA)))
            evals = filter(x -> x > 1e-15, evals)
            for α in (0.5, 2.0, 3.0)
                S_ed = log(sum(p^α for p in evals)) / (1 - α)
                S_qa = QAtlas.fetch(TFIM(; J=1.0, h=h), RenyiEntropy(α), OBC(N); ℓ=ℓ)
                # α < 1: ED's round-off floor at ~1e-16 propagates to
                # `~2^ℓ · (eps)^α ≈ 1e-7` in the entropy via the
                # `(1-α)^{-1}` Rényi prefactor.  α ≥ 1: round-off
                # contribution scales as `eps^α ≪ eps`, so the strict
                # 1e-10 atol holds.
                atol = α < 1 ? 1e-6 : 1e-10
                @test S_qa ≈ S_ed atol = atol
            end
        end
    end

    @testset "Renyi increases monotonically with ℓ (ℓ ≤ N/2)" begin
        # For the OBC ground state, S_α(ℓ) is non-decreasing on
        # 1 ≤ ℓ ≤ N/2; this is a sanity check (CFT log scaling at h = J,
        # bounded area-law growth otherwise).
        N, h = 16, 0.5
        model = TFIM(; J=1.0, h=h)
        prev = -1.0
        for ℓ in 1:(N ÷ 2)
            S = QAtlas.fetch(model, RenyiEntropy(2.0), OBC(N); ℓ=ℓ)
            @test S > prev - 1e-10
            prev = S
        end
    end

    @testset "thermal beta argument" begin
        # Same beta-handshake as the Von Neumann path: β = Inf must
        # match the omitted-beta call to machine precision.
        N = 12
        model = TFIM(; J=1.0, h=1.0)
        ℓ = N ÷ 2
        S_gs = QAtlas.fetch(model, RenyiEntropy(2.0), OBC(N); ℓ=ℓ)
        S_inf = QAtlas.fetch(model, RenyiEntropy(2.0), OBC(N); ℓ=ℓ, beta=Inf)
        @test S_gs ≈ S_inf atol = 1e-12

        # Finite β adds thermal mixing → larger entropy than the GS.
        S_thermal = QAtlas.fetch(model, RenyiEntropy(2.0), OBC(N); ℓ=ℓ, beta=1.0)
        @test S_thermal > S_gs
    end

    @testset "input validation" begin
        @test_throws ArgumentError QAtlas.fetch(
            TFIM(; J=1.0, h=1.0), RenyiEntropy(2.0), OBC(10); ℓ=0
        )
        @test_throws ArgumentError QAtlas.fetch(
            TFIM(; J=1.0, h=1.0), RenyiEntropy(2.0), OBC(10); ℓ=10
        )
        # The α = 1 rejection lives in the RenyiEntropy constructor itself.
        @test_throws ArgumentError RenyiEntropy(1.0)
    end
end
