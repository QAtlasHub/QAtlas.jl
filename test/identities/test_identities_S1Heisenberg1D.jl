using Test
using QAtlas
using QAtlas: S1Heisenberg1D, OBC

# Self-validation harness applied to spin-1 Heisenberg (Haldane chain)
# at OBC.  Dense-ED gives finite-N exact thermal data so residuals
# should be at the eigendecomposition noise floor (≪ 1e-8).

@testset "S1Heisenberg1D ε = f + T·s and c_v = -β² ∂ε/∂β  — OBC(N=4)" begin
    model = S1Heisenberg1D(; J=1.0)
    βs = [0.5, 1.0, 2.0]
    results = verify_thermodynamic_identities(model, OBC(4); βs=βs)

    @test length(results) == 6
    @test all(r.status === :pass for r in results)
    @test any(occursin("Gibbs", r.identity) for r in results)
    @test any(occursin("c_v", r.identity) for r in results)

    for r in results
        @test r.abs_err < 1e-7
    end
end

@testset "S1Heisenberg1D — identities at non-unit J (J = 0.7)" begin
    model = S1Heisenberg1D(; J=0.7)
    βs = [0.5, 2.0]
    results = verify_thermodynamic_identities(model, OBC(4); βs=βs)
    @test all(r.status === :pass for r in results)
end
