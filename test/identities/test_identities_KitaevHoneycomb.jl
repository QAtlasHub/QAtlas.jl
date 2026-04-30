using Test
using QAtlas
using QAtlas:
    KitaevHoneycomb, Energy, FreeEnergy, ThermalEntropy, SpecificHeat, OBC, Infinite
using ForwardDiff

# Self-validation harness (`test/util/thermodynamic_identities.jl`)
# applied to KitaevHoneycomb at the gapless isotropic point.  Two
# identities ship in DEFAULT_IDENTITIES:
#
#   - Gibbs:                 ε = f + T·s
#   - SpecificHeat-from-ε:   c_v = -β² ∂ε/∂β   (ForwardDiff)
#
# Notes:
#   - Both rely on `Energy(:per_site)`, which Kitaev exposes natively at
#     every BC (no `:total` conversion required).
#   - The harness runs cleanly on `Infinite()` because no boundary kwargs
#     are needed.  For `OBC` the harness invokes `fetch(model, ...,
#     OBC(N))` without forwarding `Lx, Ly` kwargs, which Kitaev's OBC
#     dispatch insists on — so we run the OBC identities by hand.

@testset "KitaevHoneycomb ε = f + T·s and c_v = -β² ∂ε/∂β  — Infinite()" begin
    model = KitaevHoneycomb(; Kx=1.0, Ky=1.0, Kz=1.0)
    βs = [0.5, 1.0, 2.0]
    results = verify_thermodynamic_identities(model, Infinite(); βs=βs)

    # 2 identities × 3 βs
    @test length(results) == 6
    @test all(r.status === :pass for r in results)

    # Spot-check both identities ran (didn't all skip).
    @test any(occursin("Gibbs", r.identity) for r in results)
    @test any(occursin("c_v", r.identity) for r in results)

    # Tolerance: matter-sector quantities go through QuadGK; the harness
    # default is `(rtol=1e-8, atol=1e-10)`.  Numerical residuals should
    # comfortably sit below 1e-7 here.
    for r in results
        @test r.abs_err < 1e-7
    end
end

@testset "KitaevHoneycomb identities — OBC manual (harness needs N kwarg)" begin
    model = KitaevHoneycomb(; Kx=1.0, Ky=1.0, Kz=1.0)
    Lx, Ly = 4, 4
    for β in (0.5, 1.0, 2.0)
        ε = QAtlas.fetch(model, Energy(:per_site), OBC(0); Lx=Lx, Ly=Ly, beta=β)
        f = QAtlas.fetch(model, FreeEnergy(), OBC(0); Lx=Lx, Ly=Ly, beta=β)
        s = QAtlas.fetch(model, ThermalEntropy(), OBC(0); Lx=Lx, Ly=Ly, beta=β)
        @test ε ≈ f + s / β atol = 1e-9

        dε = ForwardDiff.derivative(
            b -> QAtlas.fetch(model, Energy(:per_site), OBC(0); Lx=Lx, Ly=Ly, beta=b), β
        )
        c = QAtlas.fetch(model, SpecificHeat(), OBC(0); Lx=Lx, Ly=Ly, beta=β)
        @test c ≈ -β^2 * dε rtol = 1e-3
    end
end
