# test/core/test_abq_conformance.jl
#
# #734's acceptance criterion, as a mechanism instead of a promise.
#
# The issue states it as: QAtlas is a conforming atlas iff, for the quantities it
# claims, AbstractQAtlas's universal relations hold on QAtlas's values.  Nothing
# computed the left-hand side of that "iff".  Which relations get materialized is
# a hand-maintained list spread across identity_registry.jl / bound_registry.jl /
# response_registry.jl, and nothing related it back to `ABQ.all_relations()`.
#
# This test closes that hole.  It checks no physics — the generated checks do
# that.  It checks that the set of universal relations we materialize is the set
# we COULD materialize, minus an explicit allow-list that names what is missing.
#
# ── What "could" means, and why it is not `quantities(rel)` ──────────────────
#
# `quantities(rel)` is a MERGED set: it folds in `also_constrains` and the
# auto-derived family members, so a relation can report a quantity while having
# no slot bound to a type.  `EntropyNonNegativity` is the worked example —
# `quantities` says `(VonNeumannEntropy,)`, but `variable_slots` says
# `((:S, nothing),)`, and `bound!` rejects it outright:
#
#     EntropyNonNegativity has no type-keyed slot, so there is nothing to fetch
#
# So the metric here is the generators' own discriminator, lifted verbatim from
# `core/response.jl`:
#
#     _quantity_slot(T) = (T isa Type && T <: AbstractQuantity) ? T : nothing
#
# A slot typed with a quantity must be FETCHED, so a hub has to implement it.  A
# slot typed with something else (`InverseTemperature`) comes from the sweep.  An
# untyped slot needs `derived =`, which only `@response` supplies — `@bound`
# rejects it.  So: materializable on a hub iff it has at least one QUANTITY slot
# and every quantity slot is implemented there.
#
# A relation with no quantity slot at all cannot be hung on a hub by any
# generator, and no amount of QAtlas-side work changes that — the block is
# upstream, in how AbstractQAtlas declares the relation.  That is where the whole
# entanglement family sits: `Subadditivity`'s S_A / S_B / S_AB are all untyped,
# because there is no region-parameterized entropy quantity to type them WITH.
#
# The metric is self-checking: it reproduces exactly the seven relations that are
# actually wired and emitting checks.  Two looser variants tried first did not —
# one over-counted by using `quantities()`, the other under-counted by treating
# `InverseTemperature` as something a hub must implement.
using Test
using QAtlas

const ABQ = QAtlas.AbstractQAtlas

# ── hubs ─────────────────────────────────────────────────────────────────────

"Map every `(model, bc)` hub to the set of quantity types the registry implements there."
function _hub_quantities()
    d = Dict{Tuple{Type,Type},Set{Type}}()
    for row in QAtlas.REGISTRY
        push!(get!(d, (row.model, row.bc), Set{Type}()), row.quantity)
    end
    return d
end

# A hub satisfies a slot if it implements the slot type, a subtype of it, or a
# supertype (family slots such as `Susceptibility` cover `Susceptibility{(:x,:x)}`).
_satisfies(have::Set{Type}, q::Type) = any(h -> h <: q || q <: h, have)

# Lifted from core/response.jl so the two cannot drift apart.
_is_quantity_slot(T) = (T isa Type && T <: QAtlas.AbstractQuantity)

function _quantity_slots(rel)
    return [(n, T) for (n, T) in ABQ.variable_slots(rel) if _is_quantity_slot(T)]
end
_untyped_slots(rel) = [n for (n, T) in ABQ.variable_slots(rel) if T === nothing]

"How many hubs implement every QUANTITY slot of `rel`.  0 when it has none."
function _materializable_hub_count(rel, hubs)
    qs = _quantity_slots(rel)
    isempty(qs) && return 0
    return count(have -> all(((_, T),) -> _satisfies(have, T), qs), values(hubs))
end

# ── what is wired today ──────────────────────────────────────────────────────

"""
Relations QAtlas materializes, read back off the edge stores rather than from a
second hand-maintained list — a name in a comment must not be able to pass for a
wiring.
"""
function _wired_relations()
    wired = Set{Symbol}()
    for spec in QAtlas.EDGE_STORES, edge in spec.store
        for field in (:relation, :inequality)
            hasproperty(edge, field) || continue
            push!(wired, nameof(typeof(getproperty(edge, field))))
        end
    end
    return wired
end

# One edge kind delegates to an ABQ relation without storing the relation object,
# so it cannot be read back.  Recorded explicitly rather than papered over; if
# `TupleIdentityEdge` ever gains a `relation` field this entry becomes stale and
# the staleness assertion below will say so.
const WIRED_WITHOUT_A_STORED_RELATION = Dict{Symbol,String}(
    :FreeEnergyLegendre =>
        "the `:gibbs` @identity calls solve(…, Val(:F)) inside its " *
        "check closure; TupleIdentityEdge has no `relation` field",
)

# ── the allow-list: materializable, not wired, with the reason ───────────────

"""
Relations that COULD be hung on a hub today — every QUANTITY slot is implemented
somewhere — but are not.  The value names what is missing, so the list reads as a
work queue rather than a suppression list.  Remove an entry in the same PR that
wires the relation.

This list is the whole of QAtlas's remaining #734 surface, and it is three
entries long.  The large body of unmaterialized ABQ relations is not here because
it is not QAtlas's to fix: 68 of the 129 have no quantity slot at all (blocked
upstream), and 51 more name quantities this atlas does not compute.
"""
const MATERIALIZABLE_BUT_UNWIRED = Dict{Symbol,String}(
    :SusceptibilityFDT => "needs var_M, the magnetization variance — the sibling of the \
                           var_E supplier that SpecificHeatFDT already uses"
)

@testset "AbstractQAtlas relation conformance" begin
    hubs = _hub_quantities()
    @test !isempty(hubs)

    rels = ABQ.all_relations()
    @test !isempty(rels)

    wired = union(_wired_relations(), keys(WIRED_WITHOUT_A_STORED_RELATION))

    materializable = Dict{Symbol,Int}()   # name => hub count
    no_quantity_slot = Symbol[]           # blocked upstream in ABQ
    unimplemented_here = Symbol[]         # QAtlas does not compute those quantities
    for rel in rels
        name = nameof(typeof(rel))
        if isempty(_quantity_slots(rel))
            push!(no_quantity_slot, name)
            continue
        end
        n = _materializable_hub_count(rel, hubs)
        n > 0 ? (materializable[name] = n) : push!(unimplemented_here, name)
    end

    # 1. Nothing materializable may be silently unwired.
    unaccounted = sort(
        collect(setdiff(keys(materializable), wired, keys(MATERIALIZABLE_BUT_UNWIRED)))
    )
    if !isempty(unaccounted)
        @error "AbstractQAtlas relations can be materialized on QAtlas hubs but are \
                neither wired nor allow-listed — wire them, or add an entry to \
                MATERIALIZABLE_BUT_UNWIRED saying what is missing" unaccounted = [
            (r, get(materializable, r, 0)) for r in unaccounted
        ]
    end
    @test isempty(unaccounted)

    # 2. The allow-list must not rot.  An entry that is now wired, or that is no
    #    longer materializable, is stale — it would otherwise keep excusing a
    #    relation nobody is thinking about any more.
    stale_wired = sort(collect(intersect(keys(MATERIALIZABLE_BUT_UNWIRED), wired)))
    if !isempty(stale_wired)
        @error "allow-listed relations are now wired — delete their \
                MATERIALIZABLE_BUT_UNWIRED entries" stale_wired
    end
    @test isempty(stale_wired)

    gone = sort([r for r in keys(MATERIALIZABLE_BUT_UNWIRED) if !haskey(materializable, r)])
    if !isempty(gone)
        @error "allow-listed relations are no longer materializable on any hub — the \
                entry describes a situation that no longer exists" gone
    end
    @test isempty(gone)

    # 3. A wired relation must actually emit checks.  Declaring an edge whose
    #    generator produces nothing is coverage on paper only — and it is easy to
    #    do by accident: `MagnetizationResponse` was declared exactly that way,
    #    because its subject `Magnetization{:z}` was implemented on no hub whose
    #    `h` is longitudinal.
    #
    #    A check that RUNS but always SKIPS is the same failure wearing a
    #    disguise, and it is the one that actually got through: wiring
    #    `SusceptibilityResponse` generated six checks that every one came back
    #    `skip — no fetch method for Magnetization{:z}`.  Generated-but-never-run
    #    is not coverage either, so `emitted` counts only checks that reach a
    #    verdict.
    # Map each edge's relation to the checks it generated, then ask only whether
    # ONE of them reaches a verdict — short-circuiting, because running all ~500
    # generated checks here would turn a bookkeeping guard into the slowest test
    # in the suite.
    checks_by_relation = Dict{Symbol,Vector{QAtlas.GeneratedCheck}}()
    all_checks = generated_checks()
    for spec in QAtlas.EDGE_STORES, edge in spec.store
        hasproperty(edge, :name) || continue
        rels_here = Symbol[]
        for field in (:relation, :inequality)
            hasproperty(edge, field) || continue
            push!(rels_here, nameof(typeof(getproperty(edge, field))))
        end
        isempty(rels_here) && continue
        mine = filter(c -> occursin("/$(edge.name)/", c.id), all_checks)
        for r in rels_here
            append!(get!(checks_by_relation, r, QAtlas.GeneratedCheck[]), mine)
        end
    end

    "Does any check of this relation reach a verdict?  A skip is not one."
    function _reaches_a_verdict(rel::Symbol)
        for c in get(checks_by_relation, rel, QAtlas.GeneratedCheck[])
            run_generated_check(c).status === :skip || return true
        end
        return false
    end

    silent = sort([
        r for
        r in intersect(keys(materializable), _wired_relations()) if !_reaches_a_verdict(r)
    ])
    if !isempty(silent)
        @error "relations are wired but no check of theirs reaches a verdict — every \
                one is absent or skipped, which is not coverage.  Fix the edge (usually \
                a missing fetch method for a derived input) or drop it and record the \
                reason in MATERIALIZABLE_BUT_UNWIRED" silent
    end
    @test isempty(silent)

    # 4. Same for the stored-relation exemption.
    for (name, _) in WIRED_WITHOUT_A_STORED_RELATION
        if name in _wired_relations()
            @error "relation is now readable off an edge store — delete its \
                    WIRED_WITHOUT_A_STORED_RELATION entry" name
        end
        @test !(name in _wired_relations())
    end

    # 5. Report the standing position, so a CI log carries the numbers #734 is
    #    ultimately measured against — including the ones that are NOT ours to
    #    close, so the split stays visible instead of folding into one
    #    discouraging total.
    @info "ABQ relation conformance" relations = length(rels) materializable = length(
        materializable
    ) wired = length(intersect(keys(materializable), wired)) allow_listed = length(
        MATERIALIZABLE_BUT_UNWIRED
    ) blocked_upstream_no_quantity_slot = length(no_quantity_slot) quantities_not_implemented_here = length(
        unimplemented_here
    )
end
