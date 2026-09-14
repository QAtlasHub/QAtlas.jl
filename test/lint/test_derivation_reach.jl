# test/lint/test_derivation_reach.jl — the reach of the AbstractQAtlas relation
# network over this atlas, pinned.
#
# `derivation_reach()` answers "which hubs could the network cross-check", and
# the answer moves with things outside this file: a quantity added to a hub, a
# relation retyped upstream, a family erased differently.  Left unpinned it can
# only shrink quietly — nothing fails when a hub stops being reachable, because
# a cross-check that is not generated is not a cross-check that fails.
#
# What is pinned is a FLOOR plus the named rows that carry something beyond the
# Gibbs relation.  Those four are the whole marginal value of the hub plane over
# the `@identity_edge :gibbs` edge that already exists, so losing one of them
# without noticing is the regression this file exists for.
#
# MEASURED (AbstractQAtlas 0.7, this atlas): 26 of 111 hubs reach at least one
# relation, over 5 distinct relations, and 24 of the 26 reach FreeEnergyLegendre
# alone.  The scaling plane, which is where the coverage actually is, is checked
# in test/generated/test_derivation_checks.jl, not here.

using QAtlas, Test
using QAtlas: derivation_reach, REGISTRY

@testset "the network's reach over this atlas has not shrunk" begin
    reach = derivation_reach()
    @test length(reach) ≥ 26

    # Exact, not a superset.  A relation newly closing here is a change in what the
    # network can check and is worth being told about; a superset test would let
    # the set drift in the direction this file exists to watch.
    rels = sort!(unique!(reduce(vcat, (r.relations for r in reach); init=Symbol[])))
    @test rels == [
        :CasimirCentralCharge,
        :FreeEnergyFromZ,
        :FreeEnergyLegendre,
        :LoschmidtRate,
        :OffCriticalEntanglementSaturation,
    ]
end

@testset "the four hubs that reach more than the Gibbs relation" begin
    reach = derivation_reach()
    byhub = Dict((r.model, r.bc) => r.relations for r in reach)

    for ((M, BC), rel) in (
        ((QAtlas.IsingSquare, PBC), :FreeEnergyFromZ),
        ((QAtlas.TFIM, Infinite), :OffCriticalEntanglementSaturation),
        ((QAtlas.TFIM, OBC), :LoschmidtRate),
        ((QAtlas.XXZ1D, Infinite), :CasimirCentralCharge),
    )
        @test haskey(byhub, (M, BC))
        @test rel in byhub[(M, BC)]
    end
end

# Positive control.  Every assertion above is "X is present", so a
# `derivation_reach()` that returned everything would pass them all and say
# nothing.  The bound in the other direction is that the reach is SMALL: most
# hubs reach no relation at all, which is the measurement that decided the
# staging (the scaling plane first, the hub plane second).
@testset "the reach is a small fraction of the atlas, not everything" begin
    reach = derivation_reach()
    hubs = length(
        unique((e.model, e.bc) for e in REGISTRY if e.method !== :not_implemented)
    )
    @test length(reach) < hubs ÷ 2

    # ...and it is dominated by one relation, so the hub plane's marginal value
    # over the existing :gibbs edge is the handful above rather than the count.
    gibbs_only = count(r -> r.relations == [:FreeEnergyLegendre], reach)
    @test gibbs_only ≥ length(reach) - 6
end
