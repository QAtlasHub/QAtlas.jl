using QAtlas, Lattice2D, LinearAlgebra, Test

# The full-ED reference needs the spin-ED utilities loaded by runtests.jl;
# when this file is included from the test harness they are already in scope.

@testset "TFIM VonNeumannEntropy — Peschel vs full ED at small N" begin
    # N = 10 → 2^N = 1024 full state vector; the SVD reference is
    # cheap and the Peschel covariance is a 20 × 20 matrix.
    N = 10
    for (J, h) in ((1.0, 0.3), (1.0, 1.0), (1.0, 3.0))
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())
        H = build_tfim(lat, J, h)
        ψ0 = eigen(Symmetric(H)).vectors[:, 1]
        model = TFIM(; J=J, h=h)
        for ℓ in 1:(N - 1)
            S_ed = entanglement_entropy(ψ0, ℓ, N)
            S_peschel = QAtlas.fetch(model, VonNeumannEntropy(), OBC(N); ℓ=ℓ)
            @test S_peschel ≈ S_ed atol = 1e-8
        end
    end
end

@testset "TFIM VonNeumannEntropy — symmetry and structure" begin
    # S(ℓ) = S(N − ℓ) for pure states (ground state); equivalent to
    # saying the Peschel covariance of the two halves gives the same
    # spectrum.
    N = 16
    model = TFIM(; J=1.0, h=1.0)  # critical, where entanglement is richest
    Ss = [QAtlas.fetch(model, VonNeumannEntropy(), OBC(N); ℓ=ℓ) for ℓ in 1:(N - 1)]
    for ℓ in 1:(N ÷ 2 - 1)
        @test Ss[ℓ] ≈ Ss[N - ℓ] rtol = 1e-10
    end
    @test argmax(Ss) == N ÷ 2 || argmax(Ss) == N ÷ 2 - 1

    # Disordered phase (h ≫ J): S ≪ 1 for any ℓ (area law; here exactly 0
    # in the h → ∞ product state, near-zero at h = 10 J).
    model_dis = TFIM(; J=1.0, h=10.0)
    S_dis = QAtlas.fetch(model_dis, VonNeumannEntropy(), OBC(N); ℓ=N ÷ 2)
    @test S_dis < 0.1

    # Ordered phase far from criticality: also area-law (ignoring the
    # near-degenerate Z_2 tunnel splitting at this finite N).
    model_ord = TFIM(; J=1.0, h=0.1)
    S_ord = QAtlas.fetch(model_ord, VonNeumannEntropy(), OBC(N); ℓ=N ÷ 2)
    @test S_ord < 1.1  # upper-bounded by ln 2 (Z_2 cat state) + finite-size noise
end

@testset "TFIM VonNeumannEntropy — Calabrese–Cardy log scaling at h = J" begin
    # S(ℓ, N) = (c/6) log( (2N/π) sin(πℓ/N) ) + s₁ for an OBC chain at
    # the Ising CFT critical point.  Reference value comes from
    # `Universality(:Ising)` so the test breaks if either the registry
    # entry or the extracted central charge drift away from the Ising
    # CFT value c = 1/2.
    c_ising = Float64(QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2).c)
    @test c_ising ≈ 0.5 atol = 1e-12

    function _extract_c_obc(N::Int, model)
        ℓs = collect(10:10:(N - 10))
        Ss = [QAtlas.fetch(model, VonNeumannEntropy(), OBC(N); ℓ=ℓ) for ℓ in ℓs]
        ξs = [log((2N / π) * sin(π * ℓ / N)) for ℓ in ℓs]
        n = length(ℓs)
        ξ̄ = sum(ξs) / n
        S̄ = sum(Ss) / n
        slope =
            sum((ξs[i] - ξ̄) * (Ss[i] - S̄) for i in 1:n) / sum((ξs[i] - ξ̄)^2 for i in 1:n)
        return 6 * slope
    end

    function _linear_fit(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real})
        n = length(xs)
        x̄ = sum(xs) / n
        ȳ = sum(ys) / n
        b = sum((xs[i] - x̄) * (ys[i] - ȳ) for i in 1:n) / sum((xs[i] - x̄)^2 for i in 1:n)
        a = ȳ - b * x̄
        return (a, b)
    end

    model = TFIM(; J=1.0, h=1.0)
    Ns = [50, 100, 200]
    cs = [_extract_c_obc(N, model) for N in Ns]

    # 1. Reference comparison at N = 100. The OBC alternating boundary
    #    correction at fixed ℓ/N and even ℓ does not average to zero,
    #    so the single-N residual is `~1.4e-2` and shrinks as 1/√N
    #    (measured empirically across N ∈ [50, 100, 200, 400, 800]).
    #    `atol = 0.02` is a reference comparison strictly stronger than
    #    the previous `rtol = 0.05`, with margin over the empirical
    #    residual.
    @test cs[2] ≈ c_ising atol = 0.02

    # 2. Convergence: the residual shrinks monotonically with N.
    @test all(abs(cs[i + 1] - c_ising) < abs(cs[i] - c_ising) for i in 1:(length(cs) - 1))

    # 3. All three finite-N estimates lie above c_ising — the alternating
    #    boundary contribution at even ℓ is positive on the OBC critical
    #    TFIM. A sign flip in the entropy would push them below.
    @test all(c -> c > c_ising, cs)

    # 4. FSS extrapolation: the residual scales as `c(N) - c_ising ≈
    #    α / √N`. A 1/√N linear fit on the intercept matches c_ising
    #    to `atol = 5e-3` — strictly stronger than the rtol = 0.05
    #    single-N check.
    invsqrtN = [1.0 / sqrt(N) for N in Ns]
    c_inf, α = _linear_fit(invsqrtN, cs)
    @test c_inf ≈ c_ising atol = 5e-3
    @test α > 0  # boundary correction has positive sign on OBC critical TFIM
end

@testset "TFIM VonNeumannEntropy — thermal beta argument" begin
    # At β = Inf the thermal covariance reduces to the ground-state
    # covariance; the two call paths must agree to machine precision.
    N = 20
    model = TFIM(; J=1.0, h=1.0)
    ℓ = N ÷ 2
    S_gs = QAtlas.fetch(model, VonNeumannEntropy(), OBC(N); ℓ=ℓ)
    S_inf = QAtlas.fetch(model, VonNeumannEntropy(), OBC(N); ℓ=ℓ, beta=Inf)
    @test S_gs ≈ S_inf atol = 1e-12

    # Finite β < ∞: thermal entropy + entanglement entropy is larger
    # than the T = 0 value (thermal mixing adds classical entropy).
    S_thermal = QAtlas.fetch(model, VonNeumannEntropy(), OBC(N); ℓ=ℓ, beta=1.0)
    @test S_thermal > S_gs
end

@testset "TFIM VonNeumannEntropy — legacy Symbol dispatch" begin
    N = 10
    S_new = QAtlas.fetch(TFIM(; J=1.0, h=1.0), VonNeumannEntropy(), OBC(N); ℓ=N ÷ 2)
    S_legacy = @test_logs (:info, r"symbol-dispatch") QAtlas.fetch(
        :TFIM, :entanglement_entropy, OBC(); N=N, J=1.0, h=1.0, ℓ=N ÷ 2
    )
    @test S_legacy ≈ S_new atol = 1e-12

    # Symbol aliases :ee, :EE, :S_vN, :EntanglementEntropy all canonicalise
    # to :entanglement_entropy. Each raw-Symbol value keys its own one-shot
    # `@info`, so wrap every call with `@test_logs` to keep CI quiet.
    for q in (:ee, :EE, :S_vN, :EntanglementEntropy)
        S_alias = @test_logs (:info, r"symbol-dispatch") QAtlas.fetch(
            :TFIM, q, OBC(); N=N, J=1.0, h=1.0, ℓ=N ÷ 2
        )
        @test S_alias ≈ S_new atol = 1e-12
    end
end

@testset "TFIM VonNeumannEntropy — input validation" begin
    @test_throws ArgumentError QAtlas.fetch(
        TFIM(; J=1.0, h=1.0), VonNeumannEntropy(), OBC(10); ℓ=0
    )
    @test_throws ArgumentError QAtlas.fetch(
        TFIM(; J=1.0, h=1.0), VonNeumannEntropy(), OBC(10); ℓ=10
    )
end
