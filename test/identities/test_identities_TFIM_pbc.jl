using Test
using QAtlas
using QAtlas: TFIM, PBC, Energy, FreeEnergy, ThermalEntropy, SpecificHeat

# Self-validation harness applied to the TFIM at PBC.  PBC adds the
# parity-projected fermion sectors on top of the BdG free-fermion
# picture, so identity violations would surface here as :fail rather
# than silent wrong numbers.

@testset "TFIM ε = f + T·s and c_v = -β² ∂ε/∂β  — PBC(N=8)" begin
    model = TFIM(; J=1.0, h=0.5)
    βs = [0.5, 1.0, 2.0]
    results = verify_thermodynamic_identities(model, PBC(8); βs=βs)

    @test length(results) == 6
    @test all(r.status === :pass for r in results)

    @test any(occursin("Gibbs", r.identity) for r in results)
    @test any(occursin("c_v", r.identity) for r in results)

    # Closed-form parity-projected free-fermion thermo: residuals should
    # be at the noise floor, well below the harness's default thresholds.
    for r in results
        @test r.abs_err < 1e-8
    end
end

@testset "TFIM PBC at the disordered phase h > J — identities still hold" begin
    model = TFIM(; J=1.0, h=1.5)
    βs = [0.7, 2.0]
    results = verify_thermodynamic_identities(model, PBC(8); βs=βs)
    @test all(r.status === :pass for r in results)
end

@testset "TFIM PBC at the critical point h = J — identities still hold" begin
    # Critical: NS gap closes as N→∞, R has a soft k=0 mode at small β.
    # Stay at moderate β to keep both sectors numerically tame.
    model = TFIM(; J=1.0, h=1.0)
    βs = [0.7, 2.0]
    results = verify_thermodynamic_identities(model, PBC(8); βs=βs)
    @test all(r.status === :pass for r in results)
end
