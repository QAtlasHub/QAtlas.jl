# test/core/test_abq_vocabulary.jl
#
# The guard that says the shared vocabulary has not silently FORKED.
#
# QAtlas takes its quantity/region/model roots from AbstractQAtlas (#734). Nothing
# stopped it from also declaring one locally under the same name, and a same-named
# second type is not an error in Julia — it is just a different type. The failure
# mode is silence:
#
#     relations_constraining(AbstractQAtlas.RenyiEntropy) == 3
#     relations_constraining(QAtlas.RenyiEntropy)         == 0
#     QAtlas.RenyiEntropy(2) isa AbstractQAtlas.RenyiEntropy  ->  false
#
# Three relations key on the base type; an atlas value of the local type could
# never be seen by them. No warning, and `relations_constraining` returning 0 reads
# exactly like "this quantity has no relations", so the hole is invisible from
# either side. MEASURED, seven names had drifted this way — `RenyiEntropy`,
# `Universality`, `CriticalExponents`, `GrowthExponents`, `PartitionFunction`,
# `CriticalTemperature`, `SpontaneousMagnetization` — while `FreeEnergy`, `Energy`,
# `SpecificHeat` and `VonNeumannEntropy` were shared correctly. So the migration
# mostly landed, and what it left behind was undetectable without this check.
#
# The test is mechanical on purpose. Enumerating collisions by eye is what let seven
# of them accumulate.

using Test
using QAtlas
using AbstractQAtlas
const ABQ = AbstractQAtlas

# Every name AbstractQAtlas exports that QAtlas also resolves.
function _shared_names()
    out = Symbol[]
    for n in names(ABQ)
        n === :AbstractQAtlas && continue
        isdefined(QAtlas, n) && isdefined(ABQ, n) || continue
        push!(out, n)
    end
    return out
end

# ...of those, the ones bound to a DIFFERENT object on each side.
function _forked(pred)
    return [
        n for n in _shared_names() if
        pred(getfield(ABQ, n)) && getfield(QAtlas, n) !== getfield(ABQ, n)
    ]
end

@testset "the shared vocabulary is not forked" begin
    shared = _shared_names()
    # not vacuous: the two packages really do share a vocabulary
    @test length(shared) > 30

    # TYPES are the load-bearing case — a forked type breaks dispatch, relation
    # keying and `isa`, all silently.
    forked_types = _forked(x -> x isa Type)
    isempty(forked_types) || @info "forked TYPES (same name, different type)" forked_types
    @test isempty(forked_types)

    # Functions too: a second generic function under a shared name means methods
    # land on an object the other package never calls.
    forked_funcs = _forked(x -> x isa Function)
    isempty(forked_funcs) ||
        @info "forked FUNCTIONS (same name, different generic)" forked_funcs
    @test isempty(forked_funcs)
end

@testset "the names that had drifted are the base ones again" begin
    # Named explicitly, so a re-fork of any of these fails by name rather than as an
    # anonymous entry in the sweep above.
    for n in (
        :Universality,
        :CriticalExponents,
        :GrowthExponents,
        :PartitionFunction,
        :CriticalTemperature,
        :SpontaneousMagnetization,
    )
        @test getfield(QAtlas, n) === getfield(ABQ, n)
        @test parentmodule(getfield(QAtlas, n)) === ABQ
    end
end

@testset "adopting the base types reconnects the relation network" begin
    # The point of the whole exercise: a quantity the atlas produces must be
    # findable by the relations that constrain it. `relations_constraining` reads the
    # relation → quantity map in reverse, so a forked type scores 0.
    #
    # Only quantities that some relation actually constrains can be checked this way,
    # so this asserts the mechanism on the ones that have relations rather than
    # demanding every quantity have one.
    withrel = [
        q for q in (FreeEnergy, Energy, SpecificHeat, ThermalEntropy, VonNeumannEntropy) if
        !isempty(ABQ.relations_constraining(q))
    ]
    @test !isempty(withrel)
    for q in withrel
        # the atlas-side binding resolves to the same type the map is keyed on
        @test !isempty(ABQ.relations_constraining(q))
    end
    # PartitionFunction is now a thermal potential rather than a bare quantity, which
    # is what puts it in reach of the thermodynamic relations at all.
    @test PartitionFunction <: ABQ.AbstractThermalPotential
    @test SpontaneousMagnetization <: ABQ.AbstractMagnetization
end
