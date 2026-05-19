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

    # 4 default identities × 3 βs; m_x identity skipped (no h field on
    # S1Heisenberg1D).
    @test length(results) == 12
    @test all(r.status !== :fail for r in results)
    @test any(occursin("Gibbs", r.identity) for r in results)
    @test any(occursin("c_v", r.identity) for r in results)
    @test count(r -> r.status === :skipped, results) == 3

    for r in filter(r -> r.status === :pass, results)
        @test r.abs_err < 1e-7
    end
end

@testset "S1Heisenberg1D — identities at non-unit J (J = 0.7)" begin
    model = S1Heisenberg1D(; J=0.7)
    βs = [0.5, 2.0]
    results = verify_thermodynamic_identities(model, OBC(4); βs=βs)
    @test all(r.status !== :fail for r in results)
end

@testset "S1Heisenberg1D — SU(2) symmetry identities" begin
    # Spin-1 Heisenberg is SU(2) symmetric: χ_xx = χ_yy = χ_zz, m_α = 0.
    model = S1Heisenberg1D(; J=1.0)
    βs = [0.5, 1.0, 2.0]
    results = verify_thermodynamic_identities(
        model, OBC(4); βs=βs, identities=SYMMETRY_IDENTITIES
    )
    @test length(results) == 15
    @test all(r.status === :pass for r in results)
    for r in results
        @test r.abs_err < 1e-12
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "S1Heisenberg1D identities — verification cards" begin
    # White-Huse 1993 DMRG Haldane-chain energy density (literature).
    verify(
        S1Heisenberg1D(; J=1.0),
        Energy(:per_site),
        Infinite();
        route=:literature_value,
        independent=-1.401484038971,
        agree_within=1e-6,
        refs=["White-Huse 1993 DMRG: e ≈ -1.401484 J (spin-1 Haldane chain)"],
    )

    # Haldane gap literature value.
    verify(
        S1Heisenberg1D(; J=1.0),
        MassGap(),
        Infinite();
        route=:literature_value,
        independent=0.41048,
        agree_within=1e-4,
        refs=["White-Huse 1993 DMRG: Haldane gap Δ ≈ 0.41048 J"],
    )
end
