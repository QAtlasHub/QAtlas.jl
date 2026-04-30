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

    # 4 default identities × 3 βs (Kubo χ_xx is opt-in).
    @test length(results) == 12
    @test all(r.status === :pass for r in results)

    @test any(occursin("Gibbs", r.identity) for r in results)
    @test any(occursin("c_v", r.identity) for r in results)
    @test any(occursin("m_x", r.identity) for r in results)

    # Closed-form parity-projected free-fermion thermo + central-diff
    # `O(δ²) ~ 1e-10` truncation: stay under 1e-7 to be safe.
    for r in results
        @test r.abs_err < 1e-7
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
    # The `m_x = -∂f/∂h` central-diff identity has truncation error
    # `O(δ² · f''')`; near criticality `f'''` is large on small lattices,
    # so loosen the harness atol from 1e-10 → 1e-5 for this slice.
    # The closed-form thermo identities (Gibbs, c_v from ε / from s)
    # remain at machine precision since they do not perturb h.
    model = TFIM(; J=1.0, h=1.0)
    βs = [0.7, 2.0]
    results = verify_thermodynamic_identities(model, PBC(8); βs=βs, atol=1e-5)
    @test all(r.status === :pass for r in results)
end
