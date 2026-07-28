using QAtlas, Test

# The :region generator (#780 step 3).
#
# What it must get right is not arithmetic -- `region_report` does that upstream
# -- but BOOKKEEPING: that the instances it generates are the ones upstream would
# discover, that every hub which can answer does answer, and that every hub which
# cannot says so instead of quietly contributing nothing.

const _REGION_CHECKS = generated_checks(; kinds=(:region,))

_regions_sig(rs) = [join(sort(collect(Int, r.sites)), "+") for r in rs]

@testset "region generator — it generates, and every instance reaches a verdict" begin
    @test !isempty(_REGION_CHECKS)
    outcomes = [QAtlas.run_generated_check(c) for c in _REGION_CHECKS]
    tally = Dict{Symbol,Int}()
    for o in outcomes
        tally[o.status] = get(tally, o.status, 0) + 1
    end
    @test get(tally, :fail, 0) == 0
    # A generator that only ever skips is not coverage -- the failure mode the
    # #734 conformance gate exists to catch.  The majority must be verdicts.
    @test get(tally, :pass, 0) > get(tally, :skip, 0)

    # every skip carries a REASON, never a bare skip
    for (c, o) in zip(_REGION_CHECKS, outcomes)
        o.status === :skip || continue
        @test !isempty(o.detail)
    end
end

@testset "region generator — the discovered relations are upstream's, not a restatement" begin
    # The generator enumerates instances from a PLACEHOLDER bag, on the premise
    # that which rows exist depends on the REGION SET alone and not on the
    # values.  Assert that premise directly: a bag of real-looking values must
    # discover exactly the same (relation, regions) pairs as a bag of zeros.
    e = only(QAtlas.REGION_EDGES)
    family = QAtlas._region_family(e.blocks)
    ABQ = QAtlas.AbstractQAtlas
    sig(rows) = sort([
        string(nameof(typeof(r.relation)), "|", join(_regions_sig(r.regions), ";")) for
        r in rows
    ])
    zeros_bag = ABQ.bag((ABQ.entanglement_entropy(Region(s...)) => 0.0 for s in family)...)
    varied_bag = ABQ.bag(
        (ABQ.entanglement_entropy(Region(s...)) => 0.1 * length(s) for s in family)...
    )
    @test sig(ABQ.region_report(zeros_bag)) == sig(ABQ.region_report(varied_bag))
end

@testset "region generator — all four discovered relations are exercised" begin
    rels = Set{String}()
    for c in _REGION_CHECKS
        parts = split(c.id, "/")
        length(parts) >= 5 && push!(rels, parts[5])
    end
    for r in ("Subadditivity", "ArakiLieb", "StrongSubadditivity", "WeakMonotonicity")
        @test r in rels
    end
end

@testset "region generator — the Jordan-Wigner guard is what limits the free-fermion hub" begin
    # TFIM/OBC is a free-fermion hub, so it refuses the ONE non-contiguous member
    # of the family, A ∪ C = {1,2,5,6} (#783).  Every instance it skips must be
    # one that needs that region -- if it ever skips an instance built only from
    # contiguous regions, something else is broken and this says so.
    e = only(QAtlas.REGION_EDGES)
    family = QAtlas._region_family(e.blocks)
    m, bc = TFIM(), OBC(e.finite_N)
    b = QAtlas._region_entropy_bag(m, bc, family, (beta=Inf,))
    @test length(b) == length(family) - 1        # exactly one region refused

    dense = QAtlas._region_entropy_bag(XXZ1D(), OBC(e.finite_N), family, (beta=Inf,))
    @test length(dense) == length(family)        # dense ED refuses nothing

    # The id lists the instance's REGIONS, not the unions the relation forms from
    # them, so a skip cannot be recognised by string-matching: the
    # strong-subadditivity instance over (A={1,2}, B={5,6}, C={3,4}) needs
    # A ∪ B = {1,2,5,6} without that ever appearing in its id.
    #
    # Nor is "names {1,2} and {5,6}" the right test.  NONE of these relations
    # forms A ∪ C: the bipartite pair needs A ∪ B, and both triples are stated
    # over A ∪ B and B ∪ C (strong subadditivity adds A ∪ B ∪ C and B; weak
    # monotonicity adds A and C).  So the triple (A={1,2}, B={3,4}, C={5,6})
    # needs only contiguous unions and must NOT skip, even though it names both
    # {1,2} and {5,6}.  Build the unions the relation actually forms.
    #
    # And the REGIONS THEMSELVES count, not only the unions: a relation needs
    # S(A) and S(B) as well as S(A∪B), so the bipartite instance over
    # (A={1,2,5,6}, B={3,4}) is refused for A alone, even though A ∪ B = {1..6}
    # is contiguous.  Measured, those two instances are exactly the ones a
    # unions-only predicate got wrong.
    parse_regions(id) = [parse.(Int, split(t, "+")) for t in split(split(id, "/")[6], "_")]
    contiguous(v) = (maximum(v) - minimum(v) + 1) == length(unique(v))
    function needs_refused_region(id)
        rs = parse_regions(id)
        unions = if length(rs) == 2
            [vcat(rs[1], rs[2])]
        else
            [vcat(rs[1], rs[2]), vcat(rs[2], rs[3]), vcat(rs...)]
        end
        return any(u -> !contiguous(sort(u)), vcat(rs, unions))
    end

    tfim_checks = filter(c -> occursin("/TFIM/", c.id), _REGION_CHECKS)
    @test !isempty(tfim_checks)
    n_skipped = 0
    for c in tfim_checks
        skipped = QAtlas.run_generated_check(c).status === :skip
        # both directions: every skip needs the refused union, and every
        # instance needing it is skipped.  One direction alone would let the hub
        # go silently blind to contiguous instances too.
        @test skipped == needs_refused_region(c.id)
        n_skipped += skipped
    end
    @test n_skipped > 0
end

@testset "region generator — Infinite hubs are excluded structurally" begin
    # A region NAMES SITES, so it needs a positioned chain; an Infinite hub is a
    # closed form in the block LENGTH and has no site 1.  Measured, such a hub
    # produces an EMPTY bag, so leaving it in would generate instances that every
    # one comes back :skip -- coverage in the listing, nothing in the report.
    @test !any(c -> occursin("/Infinite/", c.id), _REGION_CHECKS)
end

@testset "region generator — the declared exclusion is visible, not silent" begin
    excluded = filter(c -> occursin("/S1Heisenberg1D/", c.id), _REGION_CHECKS)
    @test length(excluded) == 1                   # one visible skip, not zero checks
    o = QAtlas.run_generated_check(only(excluded))
    @test o.status === :skip
    @test occursin("MEASURED", o.detail)          # the reason carries its measurement
end
