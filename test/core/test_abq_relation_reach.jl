# test/core/test_abq_relation_reach.jl — the atlas's own fetched values must reach
# AbstractQAtlas's relation network.
#
# This is a WIRING test, not a physics one: it pins that a velocity fetched from a
# QAtlas hub lands in the `v` slot of the CFT relations and that the family slot
# resolves to the right component. The Casimir formula itself is AbstractQAtlas's to
# verify, and there is no independent Casimir energy anywhere in this atlas to check
# it against — no hub registers a velocity and a Casimir energy together.
#
# Why it can exist at all: before AbstractQAtlas 0.5.0, `LuttingerVelocity` was its own
# struct, `_as_key` keys on `typeof`, and the relations declare `v::Velocity` — so a
# fetched Luttinger velocity was a bag key no relation could match. It is now
# `Velocity{:luttinger}`, a component of the family the slot names.

using Test
using QAtlas
using AbstractQAtlas:
    CasimirCentralCharge, Velocity, bag, solve, applicable_relations, relation_report

@testset "a fetched Luttinger velocity reaches the CFT relations" begin
    m = XXZ1D(; J=1.0, Δ=0.0)                       # XX point: critical, free-fermion
    u = fetch(m, LuttingerVelocity(), Infinite())
    c = fetch(m, CentralCharge(), Infinite())
    @test u ≈ 1.0                                   # u = J at Δ = 0 (= the free-fermion v_F)
    @test c ≈ 1.0                                   # Luttinger liquid

    L, dE = 32.0, -π * 1.0 * 1.0 / (6 * 32.0^2)
    b = bag(CentralCharge => c, LuttingerVelocity() => u)

    # the relation is now reachable from atlas values...
    @test any(r -> r isa CasimirCentralCharge, applicable_relations(b; dE=dE, L=L))
    # ...and the family slot resolves to the LUTTINGER component, not to some
    # anonymous velocity — the row names which velocity it used.
    rows = [r for r in relation_report(b; dE=dE, L=L) if r.relation isa CasimirCentralCharge]
    @test length(rows) == 1
    @test rows[1].subject.type === Velocity{:luttinger}
    # solving back through the family slot needs the component named (ABQ §8a)
    @test solve(
        CasimirCentralCharge(), CentralCharge, b; subject=LuttingerVelocity, dE=dE, L=L
    ) ≈ c
end
