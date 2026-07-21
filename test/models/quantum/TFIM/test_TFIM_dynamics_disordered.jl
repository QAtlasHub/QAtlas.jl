# =============================================================================
# TFIM dynamics — disordered phase (exact ξ⁻¹ r-scaling)
#
# Split out of the legacy test_TFIM_dynamics.jl (PR vs next: refactor for
# CI shard balance — original file was ~17 min on s14 because Layers 2b+2c
# ran heavy BdG Pfaffian sweeps at N=200-240 inside the same file).
#
# Helpers _build_tfim_dense, _op_site, _SZ, _SX come from test/util/
# tfim_dense_ed.jl via runtests.jl ambient include.
# =============================================================================

@testset "TFIM dynamics — disordered phase r-scaling toward exact ξ⁻¹" begin
    # In the disordered phase (h > J) the *exact* thermodynamic-limit inverse
    # correlation length follows from the BdG dispersion
    #   ω(k) = 2 √(J² + h² - 2 J h cos k)
    # via continuation k → i/ξ:
    #   ξ⁻¹ = arccosh((h² + J²) / (2 h J)).
    # The local log-derivative of ⟨σᶻ_0 σᶻ_r⟩ approaches -ξ⁻¹ from below as
    # r grows; we verify monotone convergence on a sequence of (r, r+Δ)
    # pairs.  Floating-point precision sets an upper limit on r (the
    # correlator decays exponentially), so we stop at r = 24 where
    # |C| ~ 10⁻⁸ for our test parameters.
    J = 1.0
    for h in (1.5, 2.0)
        slope_exact = -acosh((h^2 + J^2) / (2 * h * J))
        N = 200
        i0 = 50
        pairs = [(8, 10), (12, 14), (16, 20), (20, 24)]
        errs = Float64[]
        for (r1, r2) in pairs
            v1 = real(
                QAtlas.fetch(
                    TFIM(; J=J, h=h),
                    DynamicalCorrelation(:z, :z),
                    OBC(; N=N);
                    i=i0,
                    j=i0 + r1,
                    t=0.0,
                ),
            )
            v2 = real(
                QAtlas.fetch(
                    TFIM(; J=J, h=h),
                    DynamicalCorrelation(:z, :z),
                    OBC(; N=N);
                    i=i0,
                    j=i0 + r2,
                    t=0.0,
                ),
            )
            slope = (log(abs(v2)) - log(abs(v1))) / (r2 - r1)
            push!(errs, abs(slope - slope_exact))
        end
        @test issorted(errs; rev=true)   # monotone convergence
        @test errs[end] < 0.025          # exact within 2.5 % at r ≈ 22
    end
end
