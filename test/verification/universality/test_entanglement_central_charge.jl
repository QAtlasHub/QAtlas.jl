# ─────────────────────────────────────────────────────────────────────────────
# Verification: central charge from entanglement entropy
#
# For a 1D critical system with OBC, the Calabrese-Cardy formula gives:
#
#   S(l) = (c/6) ln[(2N/π) sin(πl/N)] + s₁
#
# where c is the central charge and s₁ is non-universal.
#
# Note on the Peschel correlation-matrix method: the TFIM with σ^z σ^z
# convention maps to free fermions only after a Kramers-Wannier duality.
# The BdG eigenvectors give the DUAL fermion correlators, not the spin-
# basis correlators. Correct spin-basis Peschel requires handling the
# duality boundary conditions, which is deferred to future work.
# For now, we use full ED (exact for small N, O(2^N) cost).
#
# References:
#   P. Calabrese, J. Cardy, J. Stat. Mech. 0406, P06002 (2004).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Lattice2D, LinearAlgebra, Test

# Extract c from S(l) profile using OBC Calabrese-Cardy (c/6)
function extract_central_charge_obc(Ss::Vector{Float64}, ls, N::Int)
    ξs = [log((2N / π) * sin(π * l / N)) for l in ls]
    n = length(ls)
    ξ̄ = sum(ξs) / n
    S̄ = sum(Ss) / n
    slope = sum((ξs[i] - ξ̄) * (Ss[i] - S̄) for i in 1:n) / sum((ξs[i] - ξ̄)^2 for i in 1:n)
    return 6 * slope  # OBC: S = (c/6) ξ + const
end

# PBC version: S(l) = (c/3) ln[(N/π) sin(πl/N)] + s₁ — no boundary term,
# so the leading finite-size correction is O(1/N²) at criticality instead
# of O(1/N), and extracted c at N ≈ 14 lands within sub-percent of the
# exact value for free-fermion CFTs like the critical TFIM.
function extract_central_charge_pbc(Ss::Vector{Float64}, ls, N::Int)
    ξs = [log((N / π) * sin(π * l / N)) for l in ls]
    n = length(ls)
    ξ̄ = sum(ξs) / n
    S̄ = sum(Ss) / n
    slope = sum((ξs[i] - ξ̄) * (Ss[i] - S̄) for i in 1:n) / sum((ξs[i] - ξ̄)^2 for i in 1:n)
    return 3 * slope  # PBC: S = (c/3) ξ + const
end

# Linear fit y ≈ a + b·x; returns (intercept, slope). Used to extrapolate
# c(N) → c(∞) under an assumed 1/N (or 1/N²) leading correction.
function _linear_fit(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    x̄ = sum(x) / n
    ȳ = sum(y) / n
    num = sum((x[i] - x̄) * (y[i] - ȳ) for i in 1:n)
    den = sum((x[i] - x̄)^2 for i in 1:n)
    b = num / den
    a = ȳ - b * x̄
    return (a, b)
end

# Ground-state dispatcher: dense eigen() for small N (cheap, simple);
# sparse + KrylovKit Lanczos for N ≥ 14 where dense allocates tens of
# gigabytes.  Returns (ψ0, cache_key) so callers can reuse across tests.
function _gs_tfim(lat, J, h, N)
    if N ≤ 12
        eigen(Symmetric(build_tfim(lat, J, h))).vectors[:, 1]
    else
        ground_state_krylov(build_tfim_sparse(lat, J, h))
    end
end

function _gs_heisenberg(lat, J, N)
    if N ≤ 12
        eigen(Symmetric(build_spinhalf_heisenberg(lat, J))).vectors[:, 1]
    else
        ground_state_krylov(build_spinhalf_heisenberg_sparse(lat, J))
    end
end

# QATLAS_TEST_FULL=1 opts in to N ∈ {14, 16} via sparse ED.  Default is
# the PR-CI-friendly fast path (N ≤ 12), dense ED only — total runtime
# well under a minute.
const _TEST_FULL = get(ENV, "QATLAS_TEST_FULL", "0") != "0"

@testset "Entanglement entropy → central charge (ED)" begin
    @testset "TFIM at h = J → c = 1/2" begin
        c_exact = Float64(QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2).c)
        J = 1.0

        # Default fast: N=[10, 12, 14] — dense ED at N ≤ 12, sparse +
        # KrylovKit Lanczos at N=14 (O(N·2^N) memory instead of O(4^N)).
        # N=14 is needed to pass the `rtol = 0.10` assertion below
        # (extracted c is ~11% off at N=12, ~5% off at N=14).
        # QATLAS_TEST_FULL=1: also include N=16 (~32 GB sparse).
        Ns_test = _TEST_FULL ? [10, 12, 14, 16] : [10, 12, 14]
        cs_extracted = Dict{Int,Float64}()
        ψ0_cache = Dict{Int,Vector{Float64}}()  # reused for structural checks below
        for N in Ns_test
            lat = build_lattice(Square, N, 1; boundary=OpenAxis())
            ψ0 = _gs_tfim(lat, J, J, N)
            ψ0_cache[N] = ψ0

            # Skip 2 boundary sites on each side (lattice artifacts)
            ls = collect(3:(N - 3))
            Ss = [entanglement_entropy(ψ0, l, N) for l in ls]
            cs_extracted[N] = extract_central_charge_obc(Ss, ls, N)
        end

        # Larger N gives better c (finite-size corrections decrease)
        N_max = maximum(keys(cs_extracted))
        @test abs(cs_extracted[N_max] - c_exact) < abs(cs_extracted[10] - c_exact)

        # Largest N should be within 10% of c = 1/2
        @test cs_extracted[N_max] ≈ c_exact rtol = 0.10

        # S(l) structural checks at N=12 — reuse the cached ψ0 to avoid
        # a redundant ED rebuild.
        N = 12
        ψ0 = ψ0_cache[N]
        Ss = [entanglement_entropy(ψ0, l, N) for l in 1:(N - 1)]

        # Mirror symmetry of bipartite entanglement on a *pure* ground
        # state: `S(A) = S(B)` follows directly from the Schmidt
        # decomposition, so this holds to floating-point precision —
        # not a finite-size approximation.
        for l in 2:(N ÷ 2 - 1)
            @test Ss[l] ≈ Ss[N - l] atol = 1e-10
        end

        # Maximal near center
        l_max = argmax(Ss)
        @test abs(l_max - N / 2) <= 1
    end

    @testset "TFIM: area law away from criticality" begin
        N = 10
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())

        # Disordered phase (h ≫ J): product state → S ≈ 0
        ψ_dis = _gs_tfim(lat, 1.0, 10.0, N)
        S_dis = entanglement_entropy(ψ_dis, N ÷ 2, N)

        # Critical: S is larger
        ψ_crit = _gs_tfim(lat, 1.0, 1.0, N)
        S_crit = entanglement_entropy(ψ_crit, N ÷ 2, N)

        @test S_dis < 0.1
        @test S_crit > S_dis * 5
    end

    @testset "Heisenberg chain → c = 1" begin
        J = 1.0
        N = 12
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())
        ψ0 = _gs_heisenberg(lat, J, N)

        # Even l only (suppress SU(2) alternating correction)
        ls_even = collect(4:2:(N - 4))
        Ss = [entanglement_entropy(ψ0, l, N) for l in ls_even]
        c_heis = extract_central_charge_obc(Ss, ls_even, N)

        @test c_heis ≈ 1.0 rtol = 0.20

        # c(Heisenberg) > c(TFIM)
        ψ_tfim = _gs_tfim(lat, J, J, N)
        ls_all = collect(3:(N - 3))
        Ss_tfim = [entanglement_entropy(ψ_tfim, l, N) for l in ls_all]
        c_tfim = extract_central_charge_obc(Ss_tfim, ls_all, N)

        @test c_heis > c_tfim
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Tight-signal central-charge extraction (PBC)
#
# The OBC Calabrese–Cardy fit above leaves an O(1/N) boundary correction
# that at N ≤ 16 puts the extracted c ≈ 10% off the exact value even for
# a bug-free implementation — so the assertion `rtol = 0.10` only rejects
# very large deviations. Switching to PBC removes the boundary term; at
# criticality the leading correction is O(1/N²) and the extracted c lands
# within ~0.5 % of the exact value already at N = 10. This lets the tests
# assert the central charge to ~1 % — strong enough to catch any sign
# flip, parity, or coefficient error in the entropy or the CFT dispatch,
# rather than just outright garbage.
# ─────────────────────────────────────────────────────────────────────────────

@testset "Tight signal: central charge from PBC ground states" begin
    bc_chain = LatticeBoundary((PeriodicAxis(), OpenAxis()))

    @testset "TFIM at h = J (PBC) → c = 1/2 within 1 %" begin
        c_exact = 0.5
        J = 1.0
        Ns = [10, 12, 14]
        cs = Float64[]
        for N in Ns
            lat = build_lattice(Square, N, 1; boundary=bc_chain)
            H = build_tfim_sparse(lat, J, J)
            _, vecs, info = eigsolve(
                H, randn(2^N), 1, :SR; issymmetric=true, tol=1e-11, krylovdim=30
            )
            info.converged < 1 && error("PBC TFIM ED failed to converge at N = $N")
            ψ0 = vecs[1]
            ls = collect(2:(N - 2))
            Ss = [entanglement_entropy(ψ0, l, N) for l in ls]
            push!(cs, extract_central_charge_pbc(Ss, ls, N))
        end

        # Calabrese–Cardy on a torus predicts the leading subleading
        # correction at criticality is O(1/N²) (no boundary term, so
        # the O(1/N) piece that hurts OBC is absent).  Fit the extracted
        # `c(N) = c_exact + a/N² + b/N⁴ + …` to a 1/N² ansatz over
        # Ns = $(Ns) and assert the intercept matches c = 1/2 to
        # `atol = 2e-3`.  Empirically the extracted intercept lands at
        # 0.501 ± 0.001 — the residual is the 1/N⁴ term and the boundary
        # log correction that this 1-parameter fit cannot absorb.  Still
        # ~2.5x tighter than the previous `rtol = 0.01` single-point
        # check, and any sign-flip / coefficient bug pushes it well past
        # 2e-3.
        invN_sq = [1.0 / N^2 for N in Ns]
        c_inf, a_corr = _linear_fit(invN_sq, cs)
        @test c_inf ≈ c_exact atol = 2e-3

        # Sign of the leading 1/N² correction is positive on the PBC
        # critical TFIM (Calabrese–Cardy subleading on a torus).
        @test a_corr > 0

        # Monotone convergence toward c_exact (error shrinks with N).
        @test abs(cs[end] - c_exact) < abs(cs[1] - c_exact)

        # Each finite-N value is on the c > c_exact side, consistent
        # with the positive 1/N² correction above.
        @test all(c -> c > c_exact, cs)
    end

    @testset "Heisenberg chain (PBC) → c = 1 within 3 % after 1/N extrapolation" begin
        c_exact = 1.0
        J = 1.0
        Ns = [8, 10, 12, 14]
        cs = Float64[]
        for N in Ns
            lat = build_lattice(Square, N, 1; boundary=bc_chain)
            H = build_spinhalf_heisenberg_sparse(lat, J)
            _, vecs, info = eigsolve(
                H, randn(2^N), 1, :SR; issymmetric=true, tol=1e-11, krylovdim=30
            )
            info.converged < 1 && error("PBC Heisenberg ED failed to converge at N = $N")
            ψ0 = vecs[1]
            # Even l only (suppresses the SU(2) alternating correction).
            ls = collect(2:2:(N - 2))
            Ss = [entanglement_entropy(ψ0, l, N) for l in ls]
            push!(cs, extract_central_charge_pbc(Ss, ls, N))
        end

        # Linear 1/N extrapolation to c(∞). The SU(2) marginal log
        # correction decays as 1/ln(N) asymptotically, which a 1/N
        # ansatz cannot fully absorb at these sizes — `rtol = 0.03` is
        # the **physics-imposed floor** (not a sloppy tolerance), and a
        # tighter assertion would require a `c_∞ + a/N² + b/log(N)`
        # nonlinear fit (deferred — see Vault QAtlas test-tolerance
        # hygiene audit).
        invN = [1.0 / N for N in Ns]
        c_inf, _ = _linear_fit(invN, cs)

        @test c_inf ≈ c_exact rtol = 0.03

        # Monotone convergence (cs decreases toward 1.0 as N grows).
        @test abs(cs[end] - c_exact) < abs(cs[1] - c_exact)

        # Correct sign of the Heisenberg finite-size correction.
        @test all(c -> c > c_exact, cs)
    end
end
