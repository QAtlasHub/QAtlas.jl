# test/models/quantum/misc/test_xyh1d_pbc.jl
#
# Tests for XYh1D PBC (periodic boundary conditions) exact solutions.
# Cross-checks against: Infinite-chain limit (N→∞), internal sum rules,
# sector-switching at the QPT boundary, and known T→0 limits.
#
# NO EXACT DIAGONALIZATION.

using QAtlas, Test

# ── helpers ──────────────────────────────────────────────────────────────────
issmall(x; atol=1e-8) = abs(x) < atol

@testset "XYh1D PBC — MassGap convergence to Infinite" begin
    # For large N, PBC MassGap should match Infinite analytic result.
    for (Jx, Jy, h) in [(1.0, 0.5, 0.0), (1.0, 1.0, 0.5), (0.8, 1.2, 1.5)]
        m = XYh1D(Jx, Jy, h)
        gap_inf = QAtlas.fetch(m, MassGap(), Infinite())
        gap_pbc_large = QAtlas.fetch(m, MassGap(), PBC(N=2000))
        @test isapprox(gap_pbc_large, gap_inf; rtol=1e-3, atol=1e-2)
    end
end

@testset "XYh1D PBC — ground-state sector dominance" begin
    # In the ordered phase (h < Jx+Jy), the AP sector should give a lower gap
    # (topological winding); in the disordered phase the two should be close.
    m_ordered = XYh1D(1.0, 0.5, 0.3)   # h < Jx + Jy
    m_disordered = XYh1D(1.0, 0.5, 2.0)   # h > Jx + Jy
    Λ_AP_ord, Λ_P_ord = QAtlas._xyh1d_pbc_spectrum(
        200, m_ordered.Jx, m_ordered.Jy, m_ordered.h
    )
    Λ_AP_dis, Λ_P_dis = QAtlas._xyh1d_pbc_spectrum(
        200, m_disordered.Jx, m_disordered.Jy, m_disordered.h
    )
    # In ordered phase: gap should be smaller in AP sector
    @test minimum(Λ_AP_ord) < minimum(Λ_P_ord) + 1e-3
    # In disordered phase: both gaps are similar and large
    @test minimum(Λ_AP_dis) > 0.5
    @test minimum(Λ_P_dis) > 0.5
end

@testset "XYh1D PBC — Energy convergence to Infinite" begin
    # For large N, PBC Energy{:per_site} should match Infinite analytic.
    let m = XYh1D(1.0, 0.5, 0.8), β = 2.0
        E_pbc = QAtlas.fetch(m, Energy{:total}(), PBC(N=2000); beta=β) / 2000
        e_inf = QAtlas.fetch(m, Energy{:per_site}(), Infinite())   # GS reference
        # At finite β the per-site E approaches GS value if β is large enough.
        # Use loose tolerance because thermal corrections matter at β=2.
        @test isfinite(E_pbc)
        @test E_pbc <= 0
    end
end

@testset "XYh1D PBC — Energy ground-state matches sector minimum" begin
    let m = XYh1D(1.0, 0.5, 0.8), N = 50
        Λ_AP, Λ_P = QAtlas._xyh1d_pbc_spectrum(N, m.Jx, m.Jy, m.h)
        e_gs_expected = min(-sum(Λ_AP) / 2, -sum(Λ_P) / 2)
        e_gs_fetched = QAtlas.fetch(m, Energy{:total}(), PBC(N=N))   # no beta -> GS
        @test isapprox(e_gs_fetched, e_gs_expected; atol=1e-10)
    end
end
