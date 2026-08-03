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

@testset "region generator — the free-fermion hub no longer refuses anything" begin
    # This testset used to assert the OPPOSITE: TFIM/OBC refused the one
    # non-contiguous family member, A ∪ C = {1,2,5,6}, because the Jordan-Wigner
    # string does not factorise there (#783).  #832 reinstates the string
    # explicitly, so the hub answers every member and the skip is gone.
    #
    # Asserting "nothing is skipped" alone would be a weak test — an edge that
    # generated no checks at all would satisfy it.  So the counts are pinned on
    # both sides, and the previously-refused region is checked to return a number
    # that actually DIFFERS from the fermionic one, which is what says the string
    # was reinstated rather than dropped.
    e = only(QAtlas.REGION_EDGES)
    family = QAtlas._region_family(e.blocks)
    m, bc = TFIM(), OBC(e.finite_N)
    b = QAtlas._region_entropy_bag(m, bc, family, (beta=Inf,))
    dense = QAtlas._region_entropy_bag(XXZ1D(), OBC(e.finite_N), family, (beta=Inf,))
    @test length(b) == length(family)            # nothing refused any more
    @test length(b) == length(dense)             # and it matches the dense-ED hub

    tfim_checks = filter(c -> occursin("/TFIM/", c.id), _REGION_CHECKS)
    @test !isempty(tfim_checks)
    @test all(c -> QAtlas.run_generated_check(c).status !== :skip, tfim_checks)
    @test all(c -> QAtlas.run_generated_check(c).status === :pass, tfim_checks)

    # the region that used to be refused now answers, and answers the SPIN
    # question — a route that silently returned the fermionic number would agree
    # with `FermionicEntanglementEntropy` here, and it must not
    r = Region(1, 2, 5, 6)
    S_spin = QAtlas.fetch(m, VonNeumannEntropy(), bc; region=r)
    S_ferm = QAtlas.fetch(m, FermionicEntanglementEntropy(), bc; region=r)
    @test isfinite(S_spin) && S_spin > 0
    @test abs(S_spin - S_ferm) > 0.05
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
