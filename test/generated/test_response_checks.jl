# Generated RESPONSE checks (#734 Phase B) — the derivative-supplied slice.
#
# These are the AbstractQAtlas relations that were reachable on QAtlas hubs but
# unusable because one slot is a derivative rather than a fetched value.
#
# Note the backend: `runtests.jl` loads ForwardDiff, so the extension is active
# here and these run on AD, at AD tolerance.  Without it the suite still works,
# on finite differences at `default_rtol(FiniteDifference())` — that fallback is
# covered in test/core/test_derivative.jl.

include("util_run_checks.jl")
using QAtlas: generated_checks, RESPONSES, preferred_backend, ForwardDiffBackend

@testset "generated response checks" begin
    # The test environment loads ForwardDiff, so the AD extension must be the
    # one that ran — otherwise these would silently be FD results reported at
    # an AD tolerance.
    @test preferred_backend() isa ForwardDiffBackend

    checks = generated_checks(; kinds=(:response,))
    @test !isempty(checks)
    ids = [c.id for c in checks]
    @test any(startswith("response/entropy_response/"), ids)
    @test any(startswith("response/specific_heat_from_entropy/"), ids)
    @test length(unique(ids)) == length(ids)

    run_generated_suite(checks; label="generated response checks")
end

@testset "response edges resolve their slots from the relation" begin
    @test !isempty(RESPONSES)
    for e in RESPONSES
        # The subject is derived, never declared; it must be a real slot name.
        @test e.subject in first.(QAtlas.variable_slots(e.relation))
        # Every untyped slot is supplied — an unsupplied one could not run.
        untyped = [n for (n, T) in QAtlas.variable_slots(e.relation) if T === nothing]
        @test Set(untyped) == Set(keys(e.derived))
    end
end

@testset "response! refuses what it cannot materialize" begin
    # A field derivative needs a model-parameter mechanism that does not exist.
    @test_throws ArgumentError QAtlas.∂(QAtlas.FreeEnergy, :h)
    # An untyped slot left unsupplied would be a check that never runs.
    @test_throws ArgumentError QAtlas.response!(
        :_probe_unsupplied; relation=QAtlas.EntropyResponse, derived=NamedTuple()
    )
    @test_throws ArgumentError QAtlas.response!(
        :entropy_response;
        relation=QAtlas.EntropyResponse,
        derived=(dF_dT=QAtlas.∂(QAtlas.FreeEnergy, :T),),
    )
end
