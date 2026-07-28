# ─────────────────────────────────────────────────────────────────────────────
# core/region_checks.jl — @region: the entropy inequalities over REGIONS.
#
# Step 3 of #780.  The nine entanglement relations in AbstractQAtlas are a
# different SHAPE from everything `@identity` / `@bound` / `@response` handle:
# their slots (`S_A`, `S_B`, `S_AB`, …) are not quantity TYPES to look up on a
# hub, they are the SAME quantity — `VonNeumannEntropy` — evaluated on different
# REGIONS.  `@bound` correctly rejects them; a slot named `S_AB` has no type to
# dispatch on.
#
# AbstractQAtlas already solves the discovery half.  `region_report(bag)` takes a
# bag of `entanglement_entropy(R) => S(R)` and works out for itself which
# inequalities are instantiable over the regions present: subadditivity and
# Araki–Lieb for every disjoint pair whose union is also present, weak
# monotonicity and strong subadditivity for every pairwise-disjoint triple.
#
# WHAT THIS FILE DOES NOT DO, deliberately: it does not restate that discovery.
# Which relations apply is a function of the REGION SET ALONE — the values only
# decide pass/fail, never which rows exist.  So generation builds a bag of
# PLACEHOLDER values purely to enumerate `(relation, regions)` pairs, and the
# runner builds the real bag and reads the matching row.  Upstream's rules are
# used verbatim at both ends; if AbstractQAtlas learns a new region inequality,
# it appears here with no change on this side.
#
# WHY THE BLOCKS ARE ADJACENT.  `A = 1:2, B = 3:4, C = 5:6` makes every union the
# report needs (`A∪B`, `B∪C`, `A∪B∪C`) a single contiguous interval, which is the
# only shape the free-fermion hubs can answer for (#783 — outside one interval
# the Jordan–Wigner string does not factorise and the answer is the FERMIONIC
# entropy).  Dense-ED hubs have no such restriction (#786), but the blocks are
# shared so that every hub runs the same instances and their results are
# comparable.
# ─────────────────────────────────────────────────────────────────────────────

using AbstractQAtlas: bag, entanglement_entropy, region_report

"""
    RegionEdge

A declared family of region-entropy inequality checks: a set of `blocks` to
instantiate the regions from, the chain length to place them on, and a sweep.
"""
struct RegionEdge
    name::Symbol
    blocks::Vector{Vector{Int}}
    finite_N::Int
    sweep::NamedTuple
    exclusions::Vector{Pair{Type,String}}
    atol::Float64
    notes::String
    location::String
end

const REGION_EDGES = RegionEdge[]

"""
    @region(name; blocks, finite_N, sweep=(;), exclusions=[], atol=1e-9, notes="")

Declare a family of region entropy-inequality checks.  `blocks` are the disjoint
site groups the regions are built from; every union `region_report` can form
from them is fetched and handed to it.

```julia
@region(
    :entanglement_regions,
    blocks = [[1, 2], [3, 4], [5, 6]],
    finite_N = 8,
    sweep = (beta=[Inf],),
)
```
"""
macro region(name, kwargs...)
    opts = Dict{Symbol,Any}()
    for kw in kwargs
        kw isa Expr && kw.head === :(=) ||
            throw(ArgumentError("@region: expected keyword arguments"))
        opts[kw.args[1]] = kw.args[2]
    end
    haskey(opts, :blocks) || throw(ArgumentError("@region: `blocks` is required"))
    haskey(opts, :finite_N) || throw(ArgumentError("@region: `finite_N` is required"))
    loc = string(__source__.file, ":", __source__.line)
    return quote
        push!(
            REGION_EDGES,
            RegionEdge(
                $(esc(name)),
                [collect(Int, b) for b in $(esc(opts[:blocks]))],
                $(esc(opts[:finite_N])),
                $(esc(get(opts, :sweep, :((;))))),
                Pair{Type,String}[$(esc(get(opts, :exclusions, :([]))))...],
                $(esc(get(opts, :atol, 1e-9))),
                $(esc(get(opts, :notes, ""))),
                $loc,
            ),
        )
    end
end

# ──────────────────────────────────────────────────────────────────────
# Generator — the :region kind of generated_checks()
# ──────────────────────────────────────────────────────────────────────

# Every region the report could want: each block, and every union of blocks.
# Returned sorted so the generated ids are stable across runs.
function _region_family(blocks::Vector{Vector{Int}})
    out = Vector{Int}[]
    n = length(blocks)
    for mask in 1:((1 << n) - 1)
        sites = Int[]
        for i in 1:n
            (mask >> (i - 1)) & 1 == 1 && append!(sites, blocks[i])
        end
        push!(out, sort!(unique(sites)))
    end
    sort!(out; by=s -> (length(s), s))
    return out
end

_region_tag(sites::Vector{Int}) = join(sites, "+")
_regions_tag(rs) = join([_region_tag(sort!(collect(Int, r.sites))) for r in rs], "_")

# Fetch S(R) for every region in the family.  A hub that cannot answer for some
# region (the free-fermion hubs refuse multi-interval ones) simply contributes
# nothing for it — `region_report` then instantiates whatever the remaining
# regions support, which is the honest behaviour: fewer instances, never a
# wrong one.
const _REGION_BAG_CACHE = Dict{Any,Any}()

function _region_entropy_bag(m, bc, family::Vector{Vector{Int}}, point::NamedTuple)
    key = (typeof(m), bc, family, point)
    haskey(_REGION_BAG_CACHE, key) && return _REGION_BAG_CACHE[key]
    b = _build_region_entropy_bag(m, bc, family, point)
    _REGION_BAG_CACHE[key] = b
    return b
end

function _build_region_entropy_bag(m, bc, family::Vector{Vector{Int}}, point::NamedTuple)
    pairs = Any[]
    for sites in family
        s = try
            fetch(m, VonNeumannEntropy(), bc; region=Region(sites...), point...)
        catch
            continue
        end
        (s isa Real && isfinite(s)) || continue
        push!(pairs, entanglement_entropy(Region(sites...)) => Float64(s))
    end
    return bag(pairs...)
end

function _region_exclusion_reason(e::RegionEdge, model_T::Type)
    for (T, reason) in e.exclusions
        model_T === T && return reason
    end
    return nothing
end

function region_checks()
    out = GeneratedCheck[]
    for e in REGION_EDGES
        family = _region_family(e.blocks)
        # Discovery is a function of the REGION SET alone, so a placeholder bag
        # enumerates exactly the rows a real one will produce.  0.0 everywhere
        # is fine: `region_report` reads only the KEYS to decide which rows
        # exist.  (It also means generation costs no physics at all.)
        probe = bag((entanglement_entropy(Region(s...)) => 0.0 for s in family)...)
        template = region_report(probe)
        for hub in _implemented_hubs(Type[VonNeumannEntropy])
            model_T, bc_T = hub.model, hub.bc
            # A region NAMES SITES, so it presupposes a finite, positioned
            # chain.  `Infinite` hubs are translation-invariant closed forms in
            # the block LENGTH -- there is no site 1 for `Region(1,2)` to mean --
            # so they are skipped structurally, not because they measured badly.
            # (Measured, they produce an empty bag and every instance would have
            # come back `:skip`: 72 checks reporting nothing, the exact "a skip
            # is not a verdict" failure the #734 conformance gate exists to
            # catch.)
            bc_T <: Infinite && continue
            hub_id = string("region/", e.name, "/", _kgshort(model_T), "/", _kgshort(bc_T))
            reason = _region_exclusion_reason(e, model_T)
            if reason !== nothing
                _push_excluded_check!(out, :region, hub_id, reason)
                continue
            end
            for point in _sweep_points(e.sweep)
                for row in template
                    relname = nameof(typeof(row.relation))
                    id =
                        string(hub_id, "/", relname, "/", _regions_tag(row.regions)) *
                        _point_suffix(point)
                    want_rel, want_regions = typeof(row.relation), row.regions
                    N, atol = e.finite_N, e.atol
                    runner = function ()
                        m = model_T()
                        bc = _bc_instance(bc_T; finite_N=N)
                        b = _region_entropy_bag(m, bc, family, point)
                        rows = region_report(b)
                        i = findfirst(
                            r ->
                                typeof(r.relation) === want_rel &&
                                r.regions == want_regions,
                            rows,
                        )
                        if i === nothing
                            return _skip_outcome(
                                "hub cannot supply every entropy this instance " *
                                "needs (a free-fermion hub refuses multi-interval " *
                                "regions — #783)",
                            )
                        end
                        return _bound_outcome(rows[i].slack; atol=atol)
                    end
                    push!(
                        out,
                        GeneratedCheck(
                            :region,
                            id,
                            string(
                                "region :",
                                e.name,
                                " (",
                                relname,
                                ") on ",
                                _kgshort(model_T),
                                " at ",
                                _kgshort(bc_T),
                                " over ",
                                _regions_tag(row.regions),
                                isempty(keys(point)) ? "" : " ($(point))",
                            ),
                            runner,
                        ),
                    )
                end
            end
        end
    end
    return out
end

register_check_generator!(:region, region_checks)
