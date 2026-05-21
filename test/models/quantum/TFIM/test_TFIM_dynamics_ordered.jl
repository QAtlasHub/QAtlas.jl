# =============================================================================
# TFIM dynamics — ordered phase (Pfeuty M² long-range plateau)
#
# Split out of the legacy test_TFIM_dynamics.jl (PR vs next: refactor for
# CI shard balance — original file was ~17 min on s14 because Layers 2b+2c
# ran heavy BdG Pfaffian sweeps at N=200-240 inside the same file).
#
# Helpers _build_tfim_dense, _op_site, _SZ, _SX come from test/util/
# tfim_dense_ed.jl via runtests.jl ambient include.
# =============================================================================

# Average of ⟨σᶻ_i σᶻ_{i+r}⟩ over a band of bulk r values, used as a
# numerical estimate of the long-range plateau (in the ordered phase).
function mean_far(N::Int, J::Float64, h::Float64, i0::Int)
    rs = filter(r -> i0 + r ≤ N - 5, [10, 12, 14, 16])
    isempty(rs) && error("no bulk r values fit; choose larger N or smaller i0")
    s = 0.0
    for r in rs
        s += real(
            QAtlas.fetch(
                TFIM(; J=J, h=h),
                ZZCorrelation{:dynamic}(),
                OBC(; N=N);
                i=i0,
                j=i0 + r,
                t=0.0,
            ),
        )
    end
    return s / length(rs)
end

@testset "TFIM dynamics — ordered phase Pfeuty M² convergence" begin
    # Pfeuty (1970) gives the exact thermodynamic-limit magnetisation
    #   M(h) = (1 - (h/J)^2)^{1/8}     (h ≤ J)
    # so the long-range plateau of ⟨σᶻ_i σᶻ_j⟩ in the ordered phase equals
    #   M² = (1 - (h/J)²)^{1/4}.
    # On a finite OBC chain the GS is Z₂-symmetric (⟨σᶻ⟩ = 0), so the
    # connected and disconnected correlators coincide and the bulk-bulk
    # plateau converges to M² *very* fast (corrections are exponential in
    # N times the inverse-gap scale).  We test the convergence on a
    # sequence of N values.
    J = 1.0
    for h in (0.6, 0.7, 0.8)
        M2_exact = (1 - (h / J)^2)^(1 / 4)
        errs = Float64[]
        for N in (30, 50, 80)
            i0 = N ÷ 4
            far = mean_far(N, J, h, i0)
            push!(errs, abs(far - M2_exact))
        end
        # Convergence: largest-N error must shrink by at least 10× over
        # the N range, and the largest-N error must be < 1e-4.  Both
        # are very loose for the actual exponential convergence.
        @test errs[end] < errs[1] / 10
        @test errs[end] < 1e-4
    end
end
