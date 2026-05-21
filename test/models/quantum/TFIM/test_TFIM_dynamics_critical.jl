# =============================================================================
# TFIM dynamics — criticality (CFT exponent scaling + envelope decay)
#
# Split out of the legacy test_TFIM_dynamics.jl (PR vs next: refactor for
# CI shard balance — original file was ~17 min on s14 because Layers 2b+2c
# ran heavy BdG Pfaffian sweeps at N=200-240 inside the same file).
#
# Helpers _build_tfim_dense, _op_site, _SZ, _SX come from test/util/
# tfim_dense_ed.jl via runtests.jl ambient include.
# =============================================================================

@testset "TFIM dynamics — critical CFT exponent scaling N → ∞" begin
    # Ising CFT predicts ⟨σ_0 σ_r⟩ ~ r^{-1/4} in the TL (Δ_σ = 1/8).  At
    # finite N with OBC the *effective* slope estimated from a doubling
    # ratio (r, 2r) at fixed bulk centre converges to -1/4 as N grows.  We
    # require monotonic decrease of the deviation and that it falls below
    # 0.10 by the largest N tested.
    J = 1.0
    h = 1.0
    Ns = (80, 160, 240)
    slope_exact = -1 / 4
    errs = Float64[]
    for N in Ns
        i0 = N ÷ 4
        v1 = real(
            QAtlas.fetch(
                TFIM(; J=J, h=h),
                ZZCorrelation{:dynamic}(),
                OBC(; N=N);
                i=i0,
                j=i0 + 10,
                t=0.0,
            ),
        )
        v2 = real(
            QAtlas.fetch(
                TFIM(; J=J, h=h),
                ZZCorrelation{:dynamic}(),
                OBC(; N=N);
                i=i0,
                j=i0 + 20,
                t=0.0,
            ),
        )
        slope = log2(v2 / v1)         # since 20/10 = 2
        push!(errs, abs(slope - slope_exact))
    end
    @test issorted(errs; rev=true)    # monotone convergence
    @test errs[end] < 0.10            # within 0.10 of the CFT value at N=240
end

@testset "TFIM dynamics — criticality temporal envelope decay" begin
    # At fixed bulk site, |⟨σᶻ(t) σᶻ(0)⟩|_GS decays slowly at criticality.
    # Inside the boundary-bounce time it should be monotone (after the
    # first sign of the lattice short-time wiggle has passed) and stay
    # well above any exponential-gap bound (the gap closes).
    N, J, h = 60, 1.0, 1.0
    i0 = N ÷ 2
    ts = [2.0, 4.0, 8.0]
    envelope = Float64[
        abs(
            QAtlas.fetch(
                TFIM(; J=J, h=h), ZZCorrelation{:dynamic}(), OBC(; N=N); i=i0, j=i0, t=t
            ),
        ) for t in ts
    ]
    @test issorted(envelope; rev=true)
    @test envelope[end] > 0.05        # not exponential
end
