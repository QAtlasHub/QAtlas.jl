# test/lint/test_derivation_reach.jl: the reach of the AbstractQAtlas relation
# network over this atlas, pinned.
#
# `derivation_reach()` answers "which hubs could the network cross-check", and the
# answer moves with things outside this file: a quantity added to a hub, a relation
# retyped upstream, a family erased differently.  Left unpinned it can only shrink
# quietly: nothing fails when a hub stops being reachable, because a cross-check
# that is not generated is not a cross-check that fails.
#
# MEASURED on this tree: 27 of 112 hubs reach at least one relation, over 9
# distinct relations, and 22 of the 27 reach `FreeEnergyLegendre` alone.  What the
# five remaining hubs carry is the whole marginal value of this plane over the
# `@identity_edge :gibbs` edge that already exists, so they are named.
#
# The scaling plane, which is where the coverage actually is today, is checked in
# test/generated/test_derivation_checks.jl, not here.

using QAtlas, Test
using QAtlas: derivation_reach, REGISTRY

@testset "the network's reach over this atlas has not shrunk" begin
    reach = derivation_reach()
    @test length(reach) ≥ 27

    # Exact, not a superset.  A relation newly closing here is a change in what the
    # network can check and is worth being told about; a superset test would let the
    # set drift in the direction this file exists to watch.
    rels = sort!(unique!(reduce(vcat, (r.relations for r in reach); init=Symbol[])))
    @test rels == [
        :ActivatedSpecificHeat,
        :CasimirCentralCharge,
        :FreeEnergyFromZ,
        :FreeEnergyLegendre,
        :GriffithsAutocorrelation,
        :GriffithsSpecificHeat,
        :GriffithsSusceptibility,
        :LoschmidtRate,
        :OffCriticalEntanglementSaturation,
    ]
end

@testset "the five hubs that reach more than the Gibbs relation" begin
    reach = derivation_reach()
    # Keyed by NAME: `Disordered{…}` is parametric, so the reach row holds the
    # concrete decoration and `=== Disordered` would never match it.
    byhub = Dict((nameof(r.model), nameof(r.bc)) => r.relations for r in reach)

    for ((M, BC), rels) in (
        ((:IsingSquare, :PBC), [:FreeEnergyFromZ]),
        ((:TFIM, :Infinite), [:OffCriticalEntanglementSaturation]),
        ((:TFIM, :OBC), [:LoschmidtRate]),
        ((:XXZ1D, :Infinite), [:CasimirCentralCharge]),
        (
            (:Disordered, :Infinite),
            [
                :ActivatedSpecificHeat,
                :GriffithsAutocorrelation,
                :GriffithsSpecificHeat,
                :GriffithsSusceptibility,
            ],
        ),
    )
        @test haskey(byhub, (M, BC))
        for r in rels
            @test r in byhub[(M, BC)]
        end
    end

    # The count matches the names, so a sixth hub gaining a relation fails HERE
    # rather than passing under the assertions above.
    @test count(r -> r.relations != [:FreeEnergyLegendre], reach) == 5
end

# Positive control.  Every assertion above is "X is present", so a
# `derivation_reach()` that returned everything would pass them all and say
# nothing.  The bound in the other direction is that the reach is SMALL: most hubs
# reach no relation at all, which is the measurement that decided the staging (the
# scaling plane first, the hub plane second).
@testset "the reach is a small fraction of the atlas, not everything" begin
    reach = derivation_reach()
    hubs = length(
        unique((e.model, e.bc) for e in REGISTRY if e.method !== :not_implemented)
    )
    @test length(reach) < hubs ÷ 2

    # ...and it is dominated by one relation, so the marginal value over the
    # existing :gibbs edge is the five hubs above rather than the count.
    @test count(r -> r.relations == [:FreeEnergyLegendre], reach) ≥ length(reach) - 5
end

# The reach is STRUCTURAL and over-approximates on purpose, which is worth pinning
# because the number reads like a to-do list otherwise.  Four of the nine relations
# above sit in test_abq_conformance.jl's `MATERIALIZABLE_BUT_UNWIRED` allow-list:
# their typed quantity slots close on the hub while the SUPPLIED slots they also
# need (a disorder-averaged χ(T), a c_V(T) curve) do not exist here.  The two
# tests ask different questions and must not be read as disagreeing.
@testset "a reachable relation is not therefore a wireable one" begin
    reach = derivation_reach()
    rels = Set(reduce(vcat, (r.relations for r in reach); init=Symbol[]))
    for r in (
        :GriffithsSusceptibility,
        :GriffithsSpecificHeat,
        :GriffithsAutocorrelation,
        :ActivatedSpecificHeat,
    )
        @test r in rels
    end
end
