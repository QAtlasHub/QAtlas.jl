using Test
using QAtlas

# Self-validation of Heisenberg1D OBC delegators (which forward to
# XXZ1D(Δ=1) under the hood).  Identities must hold both as a thermo
# self-check and as the SU(2)-symmetric isotropic point.

@testset "Heisenberg1D OBC — DEFAULT thermo identities" begin
    βs = [0.5, 1.0, 2.0]
    results = verify_thermodynamic_identities(Heisenberg1D(), OBC(6); βs=βs)
    # Heisenberg1D has no `h` field on its struct (J is a kwarg of fetch),
    # so `MAGNETIZATION_X_FROM_FREE_ENERGY` is :skipped.  The other 3
    # default identities (Gibbs, c_v from ε, c_v from s) all run.
    @test length(results) == 12
    @test all(r.status !== :fail for r in results)
    @test count(r -> r.status === :skipped, results) == 3
    for r in filter(r -> r.status === :pass, results)
        @test r.abs_err < 1e-7
    end
end

@testset "Heisenberg1D OBC — SU(2) symmetry identities" begin
    # Spin-1/2 Heisenberg is the canonical SU(2)-invariant chain;
    # `is_su2_symmetric(::Heisenberg1D) = true` is hard-coded in
    # test/util/thermodynamic_identities.jl.  Every χ-axis equality and
    # every m_α = 0 identity must pass.
    βs = [0.5, 1.0, 2.0]
    results = verify_thermodynamic_identities(
        Heisenberg1D(), OBC(6); βs=βs, identities=SYMMETRY_IDENTITIES
    )
    @test length(results) == 15
    @test all(r.status === :pass for r in results)
    for r in results
        @test r.abs_err < 1e-12
    end
end
