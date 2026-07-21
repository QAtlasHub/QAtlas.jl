# Generated BOUND checks (#734 Phase B) — the inequality slice of the
# generated-check protocol, materialized from the @bound edges of
# src/bound_registry.jl.
#
# These are one-sided: `slack(inequality) ≥ -atol`, with the criterion and the
# arithmetic owned by AbstractQAtlas.  Unlike an identity, a bound catches a
# failure that two equally-wrong quantities could hide from each other — a sign
# error, a bad analytic continuation, a mis-normalized thermal state.

include("util_run_checks.jl")
using QAtlas: generated_checks, BOUNDS

@testset "generated bound checks" begin
    checks = generated_checks(; kinds=(:bound,))
    @test !isempty(checks)

    ids = [c.id for c in checks]
    # Both shipped bounds materialize, and on more than one model each — the
    # "a new model gets constraint coverage for free" property, same as the
    # identity slice.
    @test any(startswith("bound/specific_heat_positivity/"), ids)
    @test any(startswith("bound/susceptibility_positivity/"), ids)
    @test any(startswith("bound/specific_heat_positivity/TFIM/"), ids)

    # The susceptibility bound's slot is the parametric FAMILY `Susceptibility`,
    # so the generator must have expanded it to concrete axis pairs rather than
    # emitting one dead check against the UnionAll.
    chi = filter(startswith("bound/susceptibility_positivity/"), ids)
    @test !isempty(chi)
    @test all(id -> occursin("Susceptibility{", id), chi)

    # Ids are the sharding/reporting key: they must be unique and stable.
    @test length(unique(ids)) == length(ids)

    run_generated_suite(checks; label="generated bound checks")
end

@testset "bound edges derive their participants from the inequality" begin
    @test !isempty(BOUNDS)
    for e in BOUNDS
        # Nothing restates the participant list: it comes from the relation's
        # typed slots, so it cannot drift from what `slack` actually consumes.
        @test !isempty(e.quantities)
        @test Set(keys(e.quantities)) ⊆ Set(QAtlas.variable_slots(e.inequality) .|> first)
        @test e.atol ≥ 0
    end
end

@testset "bound! rejects inequalities it could never materialize" begin
    # An untyped slot is a derived input (`var_E`, `S_AB`, …) the generator
    # cannot obtain from `fetch`.  Declaring such a bound would produce a check
    # that silently never runs, so it must fail loudly at declaration instead.
    @test_throws ArgumentError QAtlas.bound!(
        :_probe_untyped_slot; inequality=QAtlas.AbstractQAtlas.EntropyNonNegativity
    )
    # Duplicate names would collide the generated ids.
    @test_throws ArgumentError QAtlas.bound!(
        :specific_heat_positivity; inequality=QAtlas.SpecificHeatPositivity
    )
    # A negative tolerance is meaningless for a one-sided test.
    @test_throws ArgumentError QAtlas.bound!(
        :_probe_bad_atol; inequality=QAtlas.SpecificHeatPositivity, atol=-1.0
    )
end
