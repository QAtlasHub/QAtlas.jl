# Model-specific relations re-homed from AbstractQAtlas (#730).
#
# The load-bearing property here is REGISTRY VISIBILITY, not arithmetic.
# `@relation` registers by mutating AbstractQAtlas's state at load time, which
# does not survive QAtlas's precompilation — so the failure mode this file
# guards is silent and session-dependent: the relations work perfectly when
# called directly, while `all_relations()` / `relations_constraining(...)` never
# see them.  A test that only checked `residual(...)` would pass in exactly the
# broken state, which is why the first testset below is the important one.

using Test
using QAtlas
using QAtlas:
    EdwardsAndersonOrderParameter,
    NishimoriEnergy,
    NishimoriMagnetizationOverlap,
    AlmeidaThoulessStability,
    DrudeMobility,
    SingleBandHall,
    MODEL_SPECIFIC_RELATIONS,
    EdwardsAndersonParameter
# Reach AbstractQAtlas THROUGH QAtlas (which `import`s the module for its
# `__init__`), matching the house pattern -- `QAtlas.AbstractQAtlasModel` and
# friends elsewhere in test/.  A bare `using AbstractQAtlas` would need it as a
# direct dep of the TEST environment, which it is not: it is a dep of the
# package, and `Pkg.test`'s sandbox does not promote those.
const ABQ = QAtlas.AbstractQAtlas

@testset "#730 relations are visible in AbstractQAtlas's registry" begin
    registered = ABQ.all_relations()
    for r in MODEL_SPECIFIC_RELATIONS
        @test any(x -> x isa typeof(r), registered)
    end

    # The declared set and the re-registered list must not drift apart: a
    # relation added to model_specific.jl but forgotten in
    # MODEL_SPECIFIC_RELATIONS would be invisible in every fresh session.
    declared = Set{Type}()
    for T in (
        EdwardsAndersonOrderParameter,
        NishimoriEnergy,
        NishimoriMagnetizationOverlap,
        AlmeidaThoulessStability,
        DrudeMobility,
        SingleBandHall,
    )
        push!(declared, T)
    end
    @test Set(typeof(r) for r in MODEL_SPECIFIC_RELATIONS) == declared

    # Domain tags survive, so `all_relations(domain=...)` can still separate them.
    @test ABQ.domain(EdwardsAndersonOrderParameter()) === :spinglass
    @test ABQ.domain(DrudeMobility()) === :transport
end

@testset "#730 type-keyed slots make the relations queryable" begin
    # The point of typing q_EA: asking the quantity what constrains it now
    # reaches the re-homed relations.
    constraining = ABQ.relations_constraining(EdwardsAndersonParameter)
    @test any(r -> r isa EdwardsAndersonOrderParameter, constraining)
    @test any(r -> r isa NishimoriMagnetizationOverlap, constraining)
end

@testset "#730 physics is preserved verbatim" begin
    # Nishimori line: U = −J tanh(βJ), exact for any lattice/dimension.
    J, β = 1.3, 0.7
    @test ABQ.residual(NishimoriEnergy(); U=(-J * tanh(β * J)), J=J, β=β) ≈ 0 atol = 1e-14
    @test ABQ.solve(NishimoriEnergy(), Val(:U); J=J, β=β) ≈ -J * tanh(β * J)
    # β-or-T convention comes along for free.
    @test ABQ.residual(NishimoriEnergy(); U=(-J * tanh(β * J)), J=J, T=1 / β) ≈ 0 atol =
        1e-14

    # q = m on the Nishimori line; q_EA = overlap by definition.
    @test ABQ.residual(NishimoriMagnetizationOverlap(); q=0.42, m=0.42) == 0
    @test ABQ.residual(EdwardsAndersonOrderParameter(); q_EA=0.31, overlap=0.31) == 0

    # Drude / single-band Hall.
    @test ABQ.residual(DrudeMobility(); μ=2.0 * 3.0 / 4.0, e=2.0, τ=3.0, m=4.0) ≈ 0 atol =
        1e-14
    @test ABQ.residual(SingleBandHall(); R_H=1 / (5.0 * 2.0), n=5.0, e=2.0) ≈ 0 atol = 1e-14

    # de Almeida–Thouless is an INEQUALITY: slack is the replicon eigenvalue,
    # ≥ 0 in the replica-symmetric phase, < 0 once RSB sets in.
    @test ABQ.slack(AlmeidaThoulessStability(); βJ=0.5, sech4_avg=1.0) > 0   # RS stable
    @test ABQ.slack(AlmeidaThoulessStability(); βJ=1.0, sech4_avg=1.0) == 0  # on the AT line
    @test ABQ.slack(AlmeidaThoulessStability(); βJ=2.0, sech4_avg=1.0) < 0   # RSB
end
