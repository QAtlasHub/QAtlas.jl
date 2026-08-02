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
# untyped slot needs `derived =`, which only `@response_edge` supplies — `@bound_edge`
# rejects it.  So: materializable on a hub iff it has at least one QUANTITY slot
# and every quantity slot is implemented there.
#
# A relation with no quantity slot at all cannot be hung on a hub by a
# TYPE-DISPATCHING generator (`@identity_edge` / `@bound_edge` / `@response_edge`), because
# there is nothing to look up.  This file used to say the block was therefore
# upstream and unfixable, and named the entanglement family as the example:
# `Subadditivity`'s S_A / S_B / S_AB are untyped, since there is no
# region-parameterized entropy quantity to type them WITH.
#
# THAT WAS WRONG, and #780 is the correction.  The slots are untyped BY DESIGN:
# they are the same quantity — `VonNeumannEntropy` — on different REGIONS, and
# AbstractQAtlas keys those by VALUE (`RegionSupport`), not by type.  A generator
# that dispatches on regions instead of on types reaches them fine; `@region`
# (core/region_checks.jl) is that generator, and `_region_wired_relations()`
# below reads its relations back.  So "no quantity slot" means "not reachable by
# type dispatch", NOT "not reachable".
#
# The metric is self-checking: it reproduces exactly the relations that are
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
    return union(wired, _region_wired_relations())
end

# ── relations wired by the :region generator ────────────────────────────────
#
# A `RegionEdge` stores BLOCKS, not a relation object: which inequalities apply
# is DISCOVERED by ABQ's `region_report` from the region set, at runtime.  So the
# field scan above cannot see them, and they would read as unwired.
#
# Read them back exactly the way the generator computes them — a placeholder bag
# over the declared blocks — so the two cannot drift.  Change the blocks and this
# follows; teach ABQ a new region inequality and this picks it up with no edit.
function _region_wired_relations()
    wired = Set{Symbol}()
    for e in QAtlas.REGION_EDGES
        family = QAtlas._region_family(e.blocks)
        probe = ABQ.bag((ABQ.entanglement_entropy(Region(s...)) => 0.0 for s in family)...)
        for row in ABQ.region_report(probe)
            push!(wired, nameof(typeof(row.relation)))
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
        "the `:gibbs` @identity_edge calls solve(…, Val(:F)) inside its " *
        "check closure; TupleIdentityEdge has no `relation` field",
)

# ── the allow-list: materializable, not wired, with the reason ───────────────

"""
Relations that COULD be hung on a hub today — every QUANTITY slot is implemented
somewhere — but are not.  The value names what is missing, so the list reads as a
work queue rather than a suppression list.  Remove an entry in the same PR that
wires the relation.

This list is the whole of QAtlas's remaining #734 surface reachable by TYPE
dispatch.  The large body of unmaterialized ABQ relations is not here because it
is not reachable that way: most have no quantity slot at all, and many more name
quantities this atlas does not compute.  "No quantity slot" is not a verdict of
unreachability — the region inequalities had none and are wired now, by
dispatching on regions instead (#780) — so that bucket is a QUEUE of relations
needing a differently-shaped generator, not a closed set.
"""
const MATERIALIZABLE_BUT_UNWIRED = Dict{Symbol,String}(
    # NOT simply "needs a var_M supplier".  χ = β·var_M, and the only supplier
    # shape available here is a derivative — but Var(M) = (1/β)·∂M/∂h, so
    # supplying it that way makes the relation assert χ = ∂M/∂h, which is exactly
    # what `susceptibility_response` already checks.  It would read as a new
    # check and be the same one rearranged.
    #
    # `SpecificHeatFDT` is independent for a reason that has no analogue here: it
    # reaches C through the ENERGY while `:specific_heat_from_entropy` reaches it
    # through the ENTROPY.  There is no second route to χ.
    #
    # Closing this needs a GENUINE fluctuation — ⟨M²⟩ − ⟨M⟩² computed from the
    # state, not a rearranged derivative.  The dense-ED hubs can produce one.
    :SusceptibilityFDT => "needs a genuine var_M = ⟨M²⟩ − ⟨M⟩²; supplying it as a \
                           derivative collapses the relation onto \
                           susceptibility_response — see the note above",

    # Both entries below became materializable only when QAtlas stopped declaring its
    # own `PartitionFunction` and adopted AbstractQAtlas's (#734 leftovers). That is
    # the intended effect — a forked type is invisible to every relation keyed on the
    # base one — and it is why the count moved 10 -> 12 here.
    #
    # F = -ln(Z)/β is materializable on exactly one hub, and there it would not be a
    # check: IsingSquare registers BOTH `PartitionFunction` at PBC and `FreeEnergy` at
    # PBC through the same `:kaufman_free_fermion` routine, so the relation restates
    # the implementation. The pair that WOULD be independent is Z at PBC against
    # `FreeEnergy` at Infinite (`:onsager`) — different regions, so not instantiable as
    # it stands. DimerLattice has the same split (Z at OBC, F at Infinite).
    :FreeEnergyFromZ => "materializable only where Z and F share the routine that \
                         produced them (IsingSquare/PBC, :kaufman_free_fermion), so it \
                         restates the implementation; the independent pair (Z at PBC \
                         vs F at Infinite/:onsager) spans two regions. Closing this \
                         needs a finite-size F from a route that did not produce Z.",

    # Z = D · tpq_weight is a statement about a TPQ construction, and the atlas has
    # neither the weight nor the Hilbert dimension as quantities — the three hubs that
    # carry a PartitionFunction (IsingSquare, DimerLattice, ChernSimons3D) supply Z and
    # nothing else the relation needs.
    :CanonicalTPQ => "needs tpq_weight and the Hilbert dimension D, neither of which \
                      is a quantity this atlas computes; the hubs with a \
                      PartitionFunction supply Z alone. Closing this needs a TPQ hub.",

    # Became materializable when QAtlas stopped declaring its own `LoschmidtAmplitude`
    # / `LoschmidtRateFunction` and adopted AbstractQAtlas's (ABQ#128). Same intended
    # effect as the PartitionFunction entries above.
    #
    # Unlike them, this one is wirable IN PRINCIPLE and is queued rather than closed.
    # MEASURED: exactly one hub carries both participants — TFIM at OBC (XXZ1D has the
    # rate function only, at Infinite, so it cannot close the identity). Two things
    # stop the generator as it stands:
    #   * both fetches need `initial::TFIM`, a MODEL-valued kwarg, and the sweep grid
    #     is splatted into `fetch` across every hub implementing the participants — a
    #     hub-specific model instance in a hub-generic grid;
    #   * `λ = -log L / N` needs `N`, an untyped slot the generator cannot fetch. The
    #     check lambda does not receive the boundary condition, so `N` would have to be
    #     restated beside `finite_N` — a second copy of a number that must agree.
    # Both are generator-shape gaps, not physics gaps.
    :LoschmidtRate => "wirable in principle but not by this generator: the only hub \
                       with both participants (TFIM/OBC) needs a model-valued \
                       `initial=` kwarg in a hub-generic sweep grid, and `N` would \
                       have to be restated beside `finite_N`. See #815.",

    # Became materializable when AbstractQAtlas 0.4.2 typed the bound's `v_LR` slot on
    # `LiebRobinsonVelocity` -- the intended effect, and the reason the count moved
    # again. It still cannot be WIRED, for the same shape of reason as
    # SusceptibilityFDT: the relation is `v <= v_LR`, and `v` is an independently
    # measured information velocity, which is not a quantity this atlas computes. The
    # TFIM `verify_bound` card in test_TFIM_status_examples.jl does supply one -- the max
    # group velocity of the BdG dispersion, computed in the test -- so the CHECK exists;
    # what is missing is a quantity to hang the relation's other slot on.
    :LiebRobinsonBound => "the `v` slot is an independently measured information \
                           velocity and no quantity names it; the TFIM card supplies one \
                           ad hoc from the dispersion. Closing this needs that \
                           measurement to become a quantity.",

    # ── the eight universal bounds (AbstractQAtlas 0.6.1 / ABQ#132) ─────────────
    #
    # These are materializable the moment they exist, because QAtlas registers every
    # bounding VALUE on a `Bound{Domain}` hub. They are unwirable for a reason that
    # is not a gap but the DESIGN: each states `measured OP fetched-bound`, and the
    # measured side is deliberately an untyped slot. No quantity in either package
    # names "the CHSH correlator I measured" or "the Lyapunov exponent I extracted",
    # and inventing one would be naming a measurement protocol, not a physical
    # quantity. `bound!` refuses an untyped slot by design, so this is the guard
    # working, not failing.
    #
    # Same shape as `LiebRobinsonBound` below, which has carried this reason since
    # AbstractQAtlas 0.4.2 — these eight are its siblings, and the reason generalises
    # exactly. Each closes when a quantity for the measured side arrives, at which
    # point the entry must be deleted (the stale-entry assertions below enforce that).
    :CHSHInequality => "the `S` slot is a MEASURED CHSH correlator and no quantity \
                        names it; only the Tsirelson/classical/PR bound VALUE is \
                        fetchable. Closing this needs the measurement to become a \
                        quantity.",
    :MerminInequality => "same shape as CHSHInequality: `M` is a measured Mermin \
                          operator value, untyped by design.",
    :LyapunovChaosBound => "`λ_L` is a Lyapunov exponent extracted from an OTOC fit; \
                            the atlas fetches the MSS ceiling, not the exponent.",
    :BekensteinEntropyBound => "`S` is the entropy of a bounded system supplied by \
                                the caller; the atlas fetches 2πRE, not S.",
    :OrthogonalizationTimeBound => "`τ` is a measured orthogonalization time; the \
                                    atlas fetches the limit it must clear. (The \
                                    energy-data forms, MargolusLevitinBound and \
                                    MandelstamTammBound, are separately unwired for \
                                    the same reason.)",
    :FastScramblingBound => "`t_scr` is a measured scrambling time; the atlas \
                             fetches the Sekino-Susskind floor.",
    :SecretKeyRateBound => "`r` is an achieved key rate; the atlas fetches the \
                            Shor-Preskill rate it must exceed.",
    :CloningFidelityBound => "`F` is a reported cloner fidelity; the atlas fetches \
                              the Buzek-Hillery optimum it must not exceed.",

    # Became materializable when AbstractQAtlas 0.5.0 folded FermiVelocity /
    # LuttingerVelocity into `Velocity{K}`: the relation's `v` slot is typed on the
    # family, and XXZ1D/Infinite is the ONE hub carrying both a Velocity component
    # (`Velocity{:luttinger}`) and a CentralCharge. Same intended effect as the two
    # PartitionFunction entries above — a name that was not part of the family could
    # not fill the slot.
    #
    # What is missing is not a quantity but a BOUNDARY CONDITION. Both remaining slots
    # are supplied values: `L`, a finite system size, and `dE = e₀(L) − e_∞`, the
    # finite-size correction to the ground-state energy density. At `Infinite` there is
    # no L, and no hub registers a finite-size energy correction. Closing this needs a
    # finite-L hub that also knows its c and its velocity — XXZ1D registers both only at
    # `Infinite`.
    :CasimirCentralCharge => "both open slots (`L`, `dE = e₀(L) − e_∞`) are supplied \
                              finite-size values, and the only hub that materializes \
                              the relation (XXZ1D/Infinite) has no L to give them.",
)

@testset "AbstractQAtlas relation conformance" begin
    hubs = _hub_quantities()
    @test !isempty(hubs)

    rels = ABQ.all_relations()
    @test !isempty(rels)

    wired = union(_wired_relations(), keys(WIRED_WITHOUT_A_STORED_RELATION))

    materializable = Dict{Symbol,Int}()   # name => hub count, by TYPE dispatch
    wired_by_region = _region_wired_relations()
    no_quantity_slot = Symbol[]           # not reachable by type dispatch, not yet wired
    unimplemented_here = Symbol[]         # QAtlas does not compute those quantities
    for rel in rels
        name = nameof(typeof(rel))
        if isempty(_quantity_slots(rel))
            # No quantity slot means "no TYPE to dispatch on", which is a
            # statement about the generator shape, not about reachability.  A
            # relation already wired by region dispatch is not blocked by
            # anything and must not be counted as if it were (#780) — that
            # bucket is a queue, and leaving a solved item in it would misreport
            # the remaining work as larger than it is.
            name in wired_by_region || push!(no_quantity_slot, name)
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

    # The :region checks carry their relation in the id rather than on an edge,
    # so they are collected by the same discovered set, not by a field scan.
    for r in _region_wired_relations()
        append!(
            get!(checks_by_relation, r, QAtlas.GeneratedCheck[]),
            filter(c -> c.kind === :region && occursin("/$(r)/", c.id), all_checks),
        )
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
    ) wired_by_type = length(intersect(keys(materializable), wired)) allow_listed = length(
        MATERIALIZABLE_BUT_UNWIRED
    ) wired_by_region = length(wired_by_region) no_type_slot_not_yet_wired = length(
        no_quantity_slot
    ) quantities_not_implemented_here = length(unimplemented_here)
end
