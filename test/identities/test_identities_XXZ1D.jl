using Test
using QAtlas
using QAtlas: XXZ1D, OBC, Energy, FreeEnergy, ThermalEntropy, SpecificHeat

# Self-validation: every (Energy, FreeEnergy, ThermalEntropy,
# SpecificHeat) at OBC routed through the dense-ED path of XXZ1D
# obeys the same Gibbs / -β² ∂ε/∂β identities as TFIM.  The dense-ED
# residuals are at machine precision (no quadrature error) so the
# AutoDiff-derived c_v should match analytic c_v exactly up to the
# step-size error of the ForwardDiff dual.

@testset "XXZ1D ε = f + T·s and c_v = -β² ∂ε/∂β — OBC(N=6)" begin
    βs = [0.5, 1.0, 2.0]
    for Δ in (-0.5, 0.0, 0.7, 1.0)
        model = XXZ1D(; J=1.0, Δ=Δ)
        results = verify_thermodynamic_identities(model, OBC(6); βs=βs)

        @test length(results) == 6
        @test all(r.status === :pass for r in results)
        @test any(occursin("Gibbs", r.identity) for r in results)
        @test any(occursin("c_v", r.identity) for r in results)
        for r in results
            @test r.abs_err < 1e-7
        end
    end
end

@testset "XXZ1D Δ = -1 (FM) — identities hold at small N" begin
    # Δ = -1 is the gapless ferromagnetic point.  Dense-ED still gives
    # exact results at finite N; check that the harness picks them up.
    model = XXZ1D(; J=1.0, Δ=-1.0)
    βs = [0.5, 2.0]
    results = verify_thermodynamic_identities(model, OBC(5); βs=βs)
    @test all(r.status === :pass for r in results)
end
