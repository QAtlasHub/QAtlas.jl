# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_Heisenberg1D_correlations_batch.jl
#
# Exact correlation function closed-form verification cards for the
# spin-½ Heisenberg chain. Pure verify(); self-contained.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Heisenberg1D — Exact Infinite Correlation Functions cards" begin
    # Exact analytical values from Sato, Shiroishi, Takahashi (2005)
    # doi:10.1016/j.nuclphysb.2005.08.045
    exact_vals = [
        -0.14771572685331508,
        0.06067976995643609,
        -0.05024862725723622,
        0.03465277698273894,
        -0.030890366644598544
    ]

    for r in 1:5
        verify(
            Heisenberg1D(),
            ZZCorrelation{:static}(),
            Infinite();
            route=:second_closed_form,
            independent=exact_vals[r],
            agree_within=1e-12,
            fetch_kw=(; i=1, j=1+r, beta=Inf),
            refs=["Sato, Shiroishi, Takahashi (2005) — Exact Bethe ansatz multiple integral correlation function for r=$r"],
        )
        
        verify(
            Heisenberg1D(),
            XXCorrelation{:connected}(),
            Infinite();
            route=:second_closed_form,
            independent=exact_vals[r],
            agree_within=1e-12,
            fetch_kw=(; i=1, j=1+r, beta=Inf),
            refs=["Sato, Shiroishi, Takahashi (2005) — Exact Bethe ansatz multiple integral correlation function for r=$r (XX connected = ZZ static for T=0)"],
        )
    end
end

@testset "Heisenberg exact correlations error handling" begin
    # Check NotImplemented for r >= 6
    @test_throws ErrorException fetch(Heisenberg1D(), ZZCorrelation{:static}(), Infinite(); i=1, j=7, beta=Inf)
    
    # Check S=1 throws properly
    @test_throws ErrorException fetch(S1Heisenberg1D(), ZZCorrelation{:static}(), Infinite(); i=1, j=2, beta=Inf)
end
