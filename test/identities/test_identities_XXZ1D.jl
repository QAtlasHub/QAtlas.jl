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

        # 4 default identities × 3 βs.  XXZ1D has no `h` field, so the
        # m_x identity is recorded as `:skipped` rather than evaluated.
        @test length(results) == 12
        @test all(r.status !== :fail for r in results)
        @test any(occursin("Gibbs", r.identity) for r in results)
        @test any(occursin("c_v", r.identity) for r in results)
        # m_x identity is skipped on XXZ1D (no h field).
        @test count(r -> r.status === :skipped, results) == 3
        # Numerical tightness on the rows that did run.
        for r in filter(r -> r.status === :pass, results)
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
    @test all(r.status !== :fail for r in results)
end

@testset "XXZ1D Δ = 1 (Heisenberg point) — SU(2) symmetry identities" begin
    # At Δ = 1 the Hamiltonian is fully SU(2) invariant, so per-site
    # susceptibilities along all three axes coincide and every uniform
    # magnetisation vanishes in the canonical ensemble.  The
    # `SYMMETRY_IDENTITIES` set covers the χ_xx = χ_yy = χ_zz pair plus
    # m_α = 0 for α ∈ {x, y, z}.
    model = XXZ1D(; J=1.0, Δ=1.0)
    βs = [0.5, 1.0, 2.0]
    results = verify_thermodynamic_identities(
        model, OBC(6); βs=βs, identities=SYMMETRY_IDENTITIES
    )

    # 5 symmetry identities × 3 βs = 15 results; all should pass at the
    # Heisenberg point (full SU(2) invariance).
    @test length(results) == 15
    @test all(r.status === :pass for r in results)
    for r in results
        @test r.abs_err < 1e-12
    end
end

@testset "XXZ1D Δ ≠ 1 — SU(2) identities skip but m_y = 0 holds" begin
    # Off the Heisenberg point the model is only U(1) × Z₂, so the SU(2)
    # axis-equality identities are :skipped.  The real-Hamiltonian
    # `m_y = 0` identity still applies and passes everywhere.
    for Δ in (-0.5, 0.0, 0.5)
        model = XXZ1D(; J=1.0, Δ=Δ)
        results = verify_thermodynamic_identities(
            model, OBC(6); βs=[1.0], identities=SYMMETRY_IDENTITIES
        )
        # 4 SU(2) identities skip + 1 m_y identity passes
        @test count(r -> r.status === :skipped, results) == 4
        @test count(r -> r.status === :pass, results) == 1
        @test all(r.status !== :fail for r in results)
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "XXZ1D identities — verification cards" begin
    # FM point Δ = -1: exact saturated ground state, e0 = -J/4.
    verify(
        XXZ1D(; J=1.0, Δ=-1.0),
        GroundStateEnergyDensity(),
        Infinite();
        route=:second_closed_form,
        independent=-0.25,
        agree_within=1e-12,
        refs=["XXZ FM point Δ=-1: aligned state exact, e0 = -J/4"],
    )

    # Δ = 1 ≡ Heisenberg1D (independent code path).
    verify(
        XXZ1D(; J=1.0, Δ=1.0),
        GroundStateEnergyDensity(),
        Infinite();
        route=:delegation_invariant,
        independent=QAtlas.fetch(Heisenberg1D(), GroundStateEnergyDensity(), Infinite()),
        agree_within=1e-12,
        refs=["XXZ1D(Δ=1) ≡ Heisenberg1D: independent code paths must agree"],
    )
end
