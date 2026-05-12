# ─────────────────────────────────────────────────────────────────────────────
# Verification: XXZ1D Luttinger parameter from exact diagonalisation
#
# Source A: QAtlas XXZ1D LuttingerParameter — closed-form Bethe-ansatz
#           result K(Δ) = π / [2 (π − γ)] with γ = arccos(Δ).
#
# Source B: Bipartite spin fluctuations of the PBC ground state of
#           H = J Σ_⟨i,j⟩ [S^x S^x + S^y S^y + Δ S^z S^z], measured via
#           sparse exact diagonalisation.
#
# Bosonisation gives (Rachel, Lehur, 2012; Song, Rachel, Le Hur, 2010)
#
#     ⟨(S^z_A)²⟩ − ⟨S^z_A⟩² = (K / π²) · ln[d(A)] + const + O(1/L²)
#
# where A = {1, …, L} and d(A) = (N/π) sin(πL/N) is the PBC chord
# length, so the slope of the variance against ln[d(A)] gives K/π².
# Sub-leading finite-size corrections are well approximated by A/N in
# the regime |Δ| ≤ 1/2; a simple 1/N linear fit of K(N) then lands
# within a few percent of the Bethe-ansatz value.
#
# Not covered: Δ = +1 (Heisenberg), where the marginally-irrelevant
# SU(2)-restoring operator produces a log correction
# K(N) = K_∞ + c / ln(N) that a linear 1/N fit undershoots by ~40 %.
# An independent Δ = +1 check is already in place via the
# `Heisenberg1D GroundStateEnergyDensity` cross-reference
# (test_XXZ1D.jl), which tests the Δ = 1 dispatch through an entirely
# different physical observable.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Lattice2D, LinearAlgebra, SparseArrays, KrylovKit, Random, Test

"""
    _sz_sub_variance(ψ, L, N) -> Float64

Return `⟨(S^z_A)²⟩ − ⟨S^z_A⟩²` for subsystem `A = {1, ..., L}` in a
spin-1/2 chain of `N` sites.

Basis convention matches `test/util/spinhalf_ed.jl`: Julia's `kron`
ordering puts site 1 in the most-significant bit, and bit 0 / bit 1 of
the basis index encode `|↑⟩` / `|↓⟩` with `S^z = +1/2` / `-1/2`. Since
`S^z_A` is diagonal in the computational basis, the variance reduces to
a weighted sum over the probability distribution `|ψ_σ|²`.
"""
function _sz_sub_variance(ψ::AbstractVector{<:Real}, L::Int, N::Int)
    dim = 2^N
    m1 = 0.0
    m2 = 0.0
    @inbounds for σ in 0:(dim - 1)
        SzA = 0.0
        for k in 1:L
            bit = (σ >> (N - k)) & 1
            SzA += (bit == 0 ? 0.5 : -0.5)
        end
        p = abs2(ψ[σ + 1])
        m1 += p * SzA
        m2 += p * SzA^2
    end
    return m2 - m1^2
end

# Linear fit: returns (intercept, slope) of y ≈ a + b·x.
function _linfit(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    x̄ = sum(x) / n
    ȳ = sum(y) / n
    num = sum((x[i] - x̄) * (y[i] - ȳ) for i in 1:n)
    den = sum((x[i] - x̄)^2 for i in 1:n)
    b = num / den
    return (ȳ - b * x̄, b)
end

# Extract K from the ground state of XXZ(Δ) at size N using the
# bipartite fluctuations formula at even L ∈ 2:2:(N÷2).
function _extract_K(N::Int, Δ::Real; J::Real=1.0)
    bc = LatticeBoundary((PeriodicAxis(), OpenAxis()))
    lat = build_lattice(Square, N, 1; boundary=bc)
    H = build_xxz_sparse(lat, J, Δ)
    _, vecs, info = eigsolve(
        H,
        randn(MersenneTwister(42 + N), 2^N),
        1,
        :SR;
        issymmetric=true,
        tol=1e-11,
        krylovdim=30,
    )
    info.converged < 1 && error("XXZ sparse ED failed to converge at (N=$N, Δ=$Δ)")
    ψ0 = vecs[1]

    Ls = collect(2:2:(N ÷ 2))
    Vars = [_sz_sub_variance(ψ0, L, N) for L in Ls]
    ξs = [log((N / π) * sin(π * L / N)) for L in Ls]
    _, slope = _linfit(ξs, Vars)
    return π^2 * slope
end

const _XXZ_Ns = [8, 10, 12, 14]

@testset "XXZ1D Luttinger K — independent ED verification" begin
    @testset "Per-Δ: extrapolated K matches the closed-form K(Δ)" begin
        # Each Δ carries its own tolerance: free-fermion and FM side
        # converge cleanly; the AF side (Δ > 0, away from Heisenberg) has
        # slower convergence because sub-leading corrections pick up a
        # non-linear-in-1/N contribution that a simple linear fit does
        # not fully remove at N ≤ 14.
        cases = [
            (-0.5, 0.03),  # K_exact = 3/2
            (0.0, 0.03),  # K_exact = 1   (XX, free fermion)
            (0.5, 0.12),  # K_exact = 3/4 (subleading AF corrections)
        ]
        for (Δ, rtol) in cases
            K_exact = QAtlas.fetch(XXZ1D(; J=1.0, Δ=Δ), LuttingerParameter(), Infinite())
            Ks = [_extract_K(N, Δ) for N in _XXZ_Ns]
            K_inf, _ = _linfit([1.0 / N for N in _XXZ_Ns], Ks)

            @test K_inf ≈ K_exact rtol = rtol

            # Error at the largest N is smaller than at the smallest:
            # sanity check on the finite-size scaling direction.
            @test abs(Ks[end] - K_exact) < abs(Ks[1] - K_exact)
        end
    end

    @testset "K is monotone decreasing in Δ at each finite N" begin
        # K = π / [2(π − γ)] with γ = arccos(Δ) is a strictly decreasing
        # function of Δ on the critical regime Δ ∈ (−1, 1). The ED
        # extraction must reproduce this ordering pointwise in N, even
        # where the individual K(N) values are offset from the infinite-N
        # closed form.
        for N in _XXZ_Ns
            K_minus = _extract_K(N, -0.5)
            K_zero = _extract_K(N, 0.0)
            K_plus = _extract_K(N, 0.5)
            @test K_minus > K_zero > K_plus
        end
    end
end
