# test/core/test_abq_conformance.jl
#
# #734's acceptance criterion, as a mechanism instead of a promise.
#
# The issue states it as: QAtlas is a conforming atlas iff, for the quantities it
# claims, AbstractQAtlas's universal relations hold on QAtlas's values.  Nothing
# computed the left-hand side of that "iff".  Which relations get checked is a
# hand-maintained list spread across identity_registry.jl / bound_registry.jl /
# response_registry.jl, and nothing related it back to `ABQ.all_relations()` — so
# a relation that became reachable (because a quantity landed on a hub, or
# because ABQ declared a new one) could sit unchecked indefinitely, with no
# symptom.  It did: the reachable count was 22 on 2026-07-21 and still 22 on
# 2026-07-27, unnoticed either way.
#
# This test closes that hole.  It does NOT check any physics — the generated
# checks do that.  It checks that the SET of universal relations we materialize
# is the set we could materialize, minus an explicit, reasoned allow-list.
#
# Reachability here means: some (model, bc) hub in the registry implements every
# QUANTITY slot of the relation.  That is necessary but not sufficient for
# checkability — a relation can also need a non-quantity slot (`var_M`, `log_d`,
# a region label) that no supplier provides.  Those are exactly what the
# allow-list records, one entry per blocking reason.
using Test
using QAtlas

const ABQ = QAtlas.AbstractQAtlas

# ── reachability ─────────────────────────────────────────────────────────────

"Map every `(model, bc)` hub to the set of quantity types the registry implements there."
function _hub_quantities()
    d = Dict{Tuple{Type,Type},Set{Type}}()
    for row in QAtlas.REGISTRY
        push!(get!(d, (row.model, row.bc), Set{Type}()), row.quantity)
    end
    return d
end

# A hub satisfies a slot if it implements the slot type, a subtype of it, or a
# supertype (the taxonomy quantities are abstract; `Susceptibility` covers
# `Susceptibility{(:x,:x)}` and vice versa at the family level).
_satisfies(have::Set{Type}, q::Type) = any(h -> h <: q || q <: h, have)

"How many hubs implement every quantity slot of `rel`.  0 for symbol-only relations."
function _reachable_hub_count(rel, hubs)
    qs = ABQ.quantities(rel)
    isempty(qs) && return 0
    return count(have -> all(q -> _satisfies(have, q), qs), values(hubs))
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

# ── the allow-list: reachable, deliberately not wired, with the reason ───────

"""
Every entry is a relation whose quantity slots ARE implemented somewhere, but
which cannot be materialized yet.  The value names what is missing, so the list
reads as a work queue rather than a suppression list.  Remove an entry in the
same PR that wires the relation.
"""
const REACHABLE_BUT_UNWIRED = Dict{Symbol,String}(
    # ── needs nothing; next to be wired ──
    :EntropyNonNegativity => "nothing blocks it — one slot, an existing quantity; \
                              same shape as the shipped SpecificHeatPositivity bound",

    # ── needs a derived input along a MODEL-parameter axis (open PR #760) ──
    :MagnetizationResponse => "needs dF/dh; core/response.jl supplies the temperature \
                               axis only, and `h` is a model parameter, not a fetch kwarg",
    :SusceptibilityResponse => "needs dM/dh; same model-parameter axis as \
                                MagnetizationResponse",

    # ── needs a derived input that has no supplier yet ──
    :MicrocanonicalTemperature => "needs dS/dE — a derivative of one quantity with \
                                   respect to another; expressible as (dS/dβ)/(dE/dβ), \
                                   both halves of which core/derivative.jl already supplies",
    :SusceptibilityFDT => "needs var_M, the magnetization variance — the sibling of the \
                           var_E supplier",
    :MaxEntropyBound => "needs log_d, the Hilbert-space dimension — model metadata \
                         (N·log d), not a fetched quantity",

    # ── needs REGION as a fetch argument (all VonNeumannEntropy on different regions) ──
    :Subadditivity => "needs region-resolved entropy: S_A, S_B, S_AB",
    :ArakiLieb => "needs region-resolved entropy: S_AB, S_A, S_B",
    :StrongSubadditivity => "needs region-resolved entropy: S_AB, S_BC, S_ABC, S_B",
    :WeakMonotonicity => "needs region-resolved entropy: S_AB, S_BC, S_A, S_C",
    :MutualInformationDefinition => "needs region-resolved entropy: S_A, S_B, S_AB",
    :EntropyMixingConcavity => "needs an ensemble-mixing vocabulary: S_mix vs S_avg",
    :HolevoMixingBound => "needs an ensemble-mixing vocabulary: S_avg, H_weights, S_mix",
    :KitaevPreskillTEE => "needs region-resolved entropy over a tripartition: \
                           S_A…S_ABC (eight slots)",
    :CFTEntanglementSlope => "needs dS/dlogℓ — the entropy's slope in SUBSYSTEM SIZE, \
                              which presupposes region as a fetch argument",
)

@testset "AbstractQAtlas relation conformance" begin
    hubs = _hub_quantities()
    @test !isempty(hubs)

    rels = ABQ.all_relations()
    @test !isempty(rels)

    wired = union(_wired_relations(), keys(WIRED_WITHOUT_A_STORED_RELATION))
    reach = Dict{Symbol,Int}()
    for rel in rels
        n = _reachable_hub_count(rel, hubs)
        n > 0 && (reach[nameof(typeof(rel))] = n)
    end

    # 1. Nothing reachable may be silently unchecked.
    unaccounted = sort(collect(setdiff(keys(reach), wired, keys(REACHABLE_BUT_UNWIRED))))
    if !isempty(unaccounted)
        @error "AbstractQAtlas relations are reachable on QAtlas hubs but neither wired \
                nor allow-listed — wire them, or add an entry to REACHABLE_BUT_UNWIRED \
                saying what is missing" unaccounted = [
            (r, get(reach, r, 0)) for r in unaccounted
        ]
    end
    @test isempty(unaccounted)

    # 2. The allow-list must not rot.  An entry that is now wired, or that is no
    #    longer reachable at all, is stale — it would otherwise keep excusing a
    #    relation nobody is thinking about any more.
    stale_wired = sort(collect(intersect(keys(REACHABLE_BUT_UNWIRED), wired)))
    if !isempty(stale_wired)
        @error "allow-listed relations are now wired — delete their \
                REACHABLE_BUT_UNWIRED entries" stale_wired
    end
    @test isempty(stale_wired)

    unreachable_entries = sort([
        r for r in keys(REACHABLE_BUT_UNWIRED) if !haskey(reach, r)
    ])
    if !isempty(unreachable_entries)
        @error "allow-listed relations are no longer reachable on any hub — the entry \
                is describing a situation that no longer exists" unreachable_entries
    end
    @test isempty(unreachable_entries)

    # 3. Same for the stored-relation exemption.
    for (name, _) in WIRED_WITHOUT_A_STORED_RELATION
        if name in _wired_relations()
            @error "relation is now readable off an edge store — delete its \
                    WIRED_WITHOUT_A_STORED_RELATION entry" name
        end
        @test !(name in _wired_relations())
    end

    # 4. Report the standing position, so a CI log carries the number that #734
    #    is ultimately measured against.
    @info "ABQ relation conformance" relations = length(rels) reachable = length(reach) wired = length(
        intersect(keys(reach), wired)
    ) allow_listed = length(REACHABLE_BUT_UNWIRED)
end
