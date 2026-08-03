using Test
using QAtlas
using QAtlas:
    REGISTRY,
    COST_VALUES,
    Infinite,
    OBC,
    PBC,
    FreeEnergy,
    Energy,
    SpecificHeat,
    ThermalEntropy,
    MassGap,
    CriticalTemperature,
    TFIM,
    IsingSquare,
    SixVertex,
    register!,
    check_cost_coherence,
    coherence_errors

# Tests for the `cost` axis of the implementation registry (`COST_VALUES` in
# `src/core/registry.jl`, C14 in `src/core/coherence.jl`).
#
# The axis exists because a route can be perfectly `:exact` and still be
# unaffordable, and a caller cannot tell those apart from `status` alone.  It
# only earns that keep if a label is a claim the implementation has to honour,
# so the properties asserted here are the ones that make a wrong label FAIL
# rather than merely read oddly:
#
#   1. the vocabulary is closed, and `:exponential` cannot be registered
#      without the `max_size` that says where it stops;
#   2. no `Infinite` row claims `:exponential` — there is no system size there
#      to be exponential in;
#   3. specific routes that were MEASURED to run a quadrature are not labelled
#      `:closed_form`.  This is the regression that matters: the tempting rule
#      "`Infinite` has no `N`, so it is closed form" is false for a large part
#      of the atlas, and a future tidy-up sweep that applies it would publish
#      "free" on routes that run a solver on every call.

@testset "COST_VALUES is a closed vocabulary" begin
    @test COST_VALUES == (:unknown, :closed_form, :polynomial, :exponential)
    @test all(e.cost in COST_VALUES for e in REGISTRY)
end

@testset "register! rejects a cost outside the vocabulary" begin
    @test_throws ArgumentError register!(TFIM, MassGap, Infinite; cost=:cheap)
end

@testset ":exponential must say where it stops" begin
    # An exponential route that will not name a cap is hiding a limit rather
    # than documenting one — the whole reason the axis was split off `status`.
    @test_throws ArgumentError register!(TFIM, MassGap, OBC; cost=:exponential)
    @test all(e.max_size !== nothing for e in REGISTRY if e.cost === :exponential)
end

@testset "no Infinite row is exponential" begin
    # `Infinite` carries no system size, so nothing can scale exponentially in
    # one.  Such a row would have either the wrong bc or the wrong cost.
    @test isempty([e for e in REGISTRY if e.cost === :exponential && e.bc <: Infinite])
end

@testset "an exponential route never shadows a cheaper one at the same hub" begin
    # The docs say the ED rows do not duplicate the closed forms — that where a
    # (model, quantity) carries both, the ED row is a finite OBC chain and the
    # cheap row is the thermodynamic limit, a DIFFERENT question rather than a
    # second route to one number.  That claim is only worth making if it stays
    # true, so it is asserted at the level it could break: the same
    # (model, quantity, bc), where two rows would be two routes to one answer.
    #
    # If this ever fires it is not necessarily a bug — someone may have added a
    # closed form beside an ED row, which is GOOD news.  The right response is to
    # make the cheap route canonical and demote or delete the ED one, not to
    # relax the test.
    bycost = Dict{Tuple{Type,Type,Type},Vector{Symbol}}()
    for e in REGISTRY
        push!(get!(bycost, (e.model, e.quantity, e.bc), Symbol[]), e.cost)
    end
    @test !isempty(bycost)                      # not vacuous
    shadowed = sort!([
        string(nameof(m), "/", nameof(q), "/", nameof(bc), " ", cs) for
        ((m, q, bc), cs) in bycost if
        (:exponential in cs) && any(c -> c in (:closed_form, :polynomial), cs)
    ],)
    isempty(shadowed) || @info "a hub offers BOTH an exponential and a cheaper \
                                route for one quantity — make the cheap one \
                                canonical rather than relaxing this test" shadowed
    @test isempty(shadowed)
end

@testset "C14 reports no cost-coherence errors" begin
    findings = check_cost_coherence()
    errs = coherence_errors(findings)
    @test isempty(errs)
    # Unclassified rows are reported as GAPS, not errors: they are a visible
    # backlog, not a violated invariant.  This asserts the reporting channel
    # works, not that the backlog is empty.
    @test all(f.check === :cost_coherence for f in findings)
end

@testset "quadrature routes at Infinite are not labelled :closed_form" begin
    # Each of these was measured under the sampling profiler to reach
    # `QuadGK.do_quadgk` on a live `fetch`; the `method` tag says otherwise on
    # two of the three (`SixVertex/FreeEnergy` is tagged `:analytic`,
    # `IsingSquare/FreeEnergy` `:onsager`), which is exactly why the tag is not
    # the discriminator.
    quad_at_infinite = [
        (TFIM, FreeEnergy),          # _tfim_quad, over the BdG dispersion
        (IsingSquare, FreeEnergy),   # the Onsager k-integral, evaluated by quadgk
        (SixVertex, FreeEnergy),     # quadgk, despite the :analytic method tag
    ]
    for (M, Q) in quad_at_infinite
        rows = [
            e for e in REGISTRY if e.model === M && e.quantity === Q && e.bc <: Infinite
        ]
        @test !isempty(rows)
        @test all(e.cost !== :closed_form for e in rows)
    end
end
