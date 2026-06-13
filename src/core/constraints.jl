# core/constraints.jl — the shared kernel of the constraint-edge layer (#697).
#
# Every meta-edge in QAtlas plays one of three roles with respect to the single
# runtime surface `fetch(model, quantity, bc; kwargs...)`:
#
#   1. describe it — `@register`   (a hub has an implementation)
#   2. route    it — `@reduces`    (delegation)
#   3. constrain it — `@identity` / `@dual` / `@limits_to` (+ `@symmetry` node
#      attributes): declared relations the implementations must satisfy.
#
# The *constrain* role is what this kernel serves.  Its central design decision
# (the integration point of #698/#699/#700/#701) is that the per-edge-type
# recipe of #697 — store + macro + query + coherence + test generator — is NOT
# implemented five times.  The recipe's shared mechanics live here ONCE:
#
#   * STORE REGISTRATION (`EDGE_STORES`): every declarative store — the legacy
#     four (REGISTRY/REALIZES/REDUCES/ABOUT) and each constraint store —
#     registers itself with accessors for its references and row locations.
#     Generic graph-wide passes (C1 reference integrity, future drift guards)
#     iterate the registration instead of hard-coding each store, so adding an
#     edge type never edits coherence.jl again.
#
#   * GENERATED-CHECK PROTOCOL (`GeneratedCheck` / `CheckOutcome`): the common
#     output type of every constraint generator.  A generator enumerates
#     (edge × implementations present in REGISTRY) and emits executable
#     cross-checks; `generated_checks()` aggregates them deterministically and
#     `test/generated/` runs them.  A new edge type buys test generation by
#     registering one generator function (`register_check_generator!`).
#
#   * HUB ENUMERATION (`_implemented_hubs`): the shared question "which
#     (model, bc) pairs implement ALL of these quantities natively?" — the
#     engine of the identity/dual/limit generators, including the
#     non-delegating filter that kills circular verification (#699/#701).
#
# Edge types stay thin, one file each: symmetry.jl (node attributes + LSM-type
# theorem checks), identity.jl (quantity↔quantity), duality.jl (model↔model
# parameter-mapped equivalence), limits.jl (model→model asymptotic limits).
# The planned `@measured` experimental anchors (#702) slot in the same way: a
# store registered here + a generator whose comparator is `|theory−exp| ≲ kσ`.
#
# Connection to instances stays one-directional by construction: edges
# describe/route/constrain `fetch`; implementations never reference the graph.

# ──────────────────────────────────────────────────────────────────────
# Store registration
# ──────────────────────────────────────────────────────────────────────

"""
    EdgeStoreSpec(name, store, references_of, location_of)

Registration record for one declarative edge store: `store` is the module-level
`const` vector, `references_of(row)` returns the bibkeys a row cites, and
`location_of(row)` renders a human-readable row locator for findings.  See
[`register_edge_store!`](@ref).
"""
struct EdgeStoreSpec
    name::Symbol
    store::Vector
    references_of::Function
    location_of::Function
end

"""
    EDGE_STORES :: Vector{EdgeStoreSpec}

Self-registration of every declarative store in the knowledge graph — the
legacy describe/route stores (`REGISTRY`, `REALIZES`, `REDUCES`, `ABOUT`) and
the constraint stores (`SYMMETRY_PROFILES`, `IDENTITIES`, `DUALITIES`,
`LIMIT_EDGES`).  Graph-wide structural passes (reference integrity C1, drift
guards) iterate THIS list, so a new edge store is covered the moment it
registers itself.
"""
const EDGE_STORES = EdgeStoreSpec[]

"""
    register_edge_store!(name, store; references_of, location_of)

Register a declarative store for generic graph-wide passes.  `references_of`
defaults to a `references`-field reader that yields `String[]` for a row type
without that field (so a future bibkey-less store contributes nothing to C1
instead of throwing inside it); `location_of` should pin the row precisely
enough that a dangling-bibkey finding is actionable.
"""
function register_edge_store!(
    name::Symbol,
    store::Vector;
    references_of::Function=row ->
        hasproperty(row, :references) ? row.references : String[],
    location_of::Function=row -> string(name, " row"),
)
    any(s -> s.name === name, EDGE_STORES) &&
        throw(ArgumentError("register_edge_store!: store :$(name) already registered"))
    push!(EDGE_STORES, EdgeStoreSpec(name, store, references_of, location_of))
    return nothing
end

# The legacy describe/route stores register here (constraints.jl loads after
# about.jl, before the constraint edge-type files).  Locations mirror the
# wording the hand-rolled C1 loops used to emit.
register_edge_store!(
    :registry, REGISTRY; location_of=e -> "$(_kgshort(e.model))/$(_kgshort(e.quantity))"
)
register_edge_store!(
    :realizes, REALIZES; location_of=r -> "realizes $(_kgshort(r.model))→:$(r.class)"
)
register_edge_store!(
    :reduces,
    REDUCES;
    location_of=r -> "reduces $(_kgshort(r.source))→$(_kgshort(r.target))",
)
register_edge_store!(:about, ABOUT; location_of=c -> "about $(_kgshort(c.model))")

# ──────────────────────────────────────────────────────────────────────
# Generated-check protocol
# ──────────────────────────────────────────────────────────────────────

"""
    CheckOutcome

The result of running one [`GeneratedCheck`](@ref).  `status` is one of:

  * `:pass`  — the check ran and the two values agree;
  * `:fail`  — the check ran and the values DISAGREE (a genuine numerical
               contradiction); `lhs`/`rhs`/`abs_err`/`rel_err` are meaningful;
  * `:skip`  — the check is declared inapplicable (an exclusion); numerics are
               `NaN`, `detail` is the reason;
  * `:error` — the runner THREW (a config/dispatch bug, NOT a numerical
               disagreement); numerics are `NaN`, `detail` carries the
               exception.  Kept distinct from `:fail` so a broken edge is not
               mis-read as a physics contradiction.

`:fail` and `:error` are both test failures, but only `:fail` means "the
physics disagrees".
"""
struct CheckOutcome
    status::Symbol
    lhs::Float64
    rhs::Float64
    abs_err::Float64
    rel_err::Float64
    detail::String
end

"""
    GeneratedCheck

One executable cross-check derived from (a constraint edge × the
implementations present in `REGISTRY`).  `kind` names the generating edge type
(`:identity` / `:dual` / `:limit` / `:symmetry`), `id` is a deterministic
identifier (stable across runs — the sharding/reporting key), and `run` is a
zero-argument callable returning a [`CheckOutcome`](@ref).
"""
struct GeneratedCheck
    kind::Symbol
    id::String
    description::String
    run::Function
end

function Base.show(io::IO, c::GeneratedCheck)
    return print(io, "GeneratedCheck(", c.kind, ", \"", c.id, "\")")
end

"""
    _outcome(lhs, rhs; rtol, atol, detail="") -> CheckOutcome

Compare two scalars with the harness's pass criterion
(`abs_err ≤ atol || rel_err ≤ rtol`).
"""
function _outcome(lhs::Real, rhs::Real; rtol::Real, atol::Real, detail::String="")
    l, r = Float64(lhs), Float64(rhs)
    abs_err = abs(l - r)
    rel_err = abs_err / max(abs(l), abs(r), eps())
    status = (abs_err ≤ atol || rel_err ≤ rtol) ? :pass : :fail
    return CheckOutcome(status, l, r, abs_err, rel_err, detail)
end

_skip_outcome(reason::String) = CheckOutcome(:skip, NaN, NaN, NaN, NaN, reason)

# Reject tolerance footguns at declaration time: an `rtol ≥ 1` passes every
# check (a common atol↔rtol slot mistake), and a negative tolerance is
# meaningless.  Shared by the identity / dual / limit factories.
function _check_tolerances(who::Symbol, rtol::Real, atol::Real)
    (0 ≤ rtol < 1) || throw(
        ArgumentError(
            "$(who)!: rtol must be in [0, 1); got $(rtol) (rtol ≥ 1 passes everything)"
        ),
    )
    atol ≥ 0 || throw(ArgumentError("$(who)!: atol must be ≥ 0; got $(atol)"))
    return nothing
end

"""
    run_generated_check(c::GeneratedCheck) -> CheckOutcome

Run `c`, converting a thrown exception into an `:error` outcome (NOT `:fail`)
whose `detail` carries the exception — a runner only ever throws on a
config/dispatch bug, so it must not be conflated with a numerical `:fail`.
Generated suites report every check rather than aborting at the first
throwing hub.
"""
function run_generated_check(c::GeneratedCheck)
    try
        return c.run()
    catch err
        return CheckOutcome(:error, NaN, NaN, NaN, NaN, string(typeof(err), ": ", err))
    end
end

"""
    CHECK_GENERATORS :: Vector{Pair{Symbol,Function}}

`kind => generator` registration, one per constraint edge type.  Each
generator is a zero-argument function returning `Vector{GeneratedCheck}`,
enumerated lazily (at call time, NOT load time) so it sees the fully-populated
`REGISTRY` and edge stores.
"""
const CHECK_GENERATORS = Pair{Symbol,Function}[]

"""
    register_check_generator!(kind::Symbol, gen::Function)

Register the test generator of a constraint edge type for
[`generated_checks`](@ref).
"""
function register_check_generator!(kind::Symbol, gen::Function)
    any(p -> first(p) === kind, CHECK_GENERATORS) &&
        throw(ArgumentError("register_check_generator!: :$(kind) already registered"))
    push!(CHECK_GENERATORS, kind => gen)
    return nothing
end

"""
    generated_checks(; kinds=nothing) -> Vector{GeneratedCheck}

Every executable cross-check the constraint layer derives from the current
registry state — the union over all registered generators, deterministically
sorted by `id`.  Pass `kinds` (e.g. `(:identity,)`) to select a subset; the
per-kind test files in `test/generated/` are exactly such selections, so the
union of the generated test suite equals this list (the universe.jl
philosophy applied to generated tests).
"""
function generated_checks(; kinds=nothing)
    out = GeneratedCheck[]
    for (kind, gen) in CHECK_GENERATORS
        (kinds === nothing || kind in kinds) || continue
        # A throwing generator is a bug, not a runtime check failure; surface
        # it loudly AND name the culprit kind (a bare rethrow would leave the
        # caller guessing which of the registered generators died).
        local emitted
        try
            emitted = gen()::Vector{GeneratedCheck}
        catch err
            error("generated_checks: the :$(kind) generator threw — $(err)")
        end
        append!(out, emitted)
    end
    sort!(out; by=c -> c.id)
    if !allunique(c.id for c in out)
        dups = [id for id in (c.id for c in out) if count(==(id), (c.id for c in out)) > 1]
        error(
            "generated_checks: duplicate check ids $(sort!(unique(dups))) — a " *
            "generator is not deterministic (ids must be unique across kinds)",
        )
    end
    return out
end

# ──────────────────────────────────────────────────────────────────────
# Shared instance / hub machinery
# ──────────────────────────────────────────────────────────────────────

"""
    _quantity_instance(Q::Type{<:AbstractQuantity}) -> AbstractQuantity

The canonical instance of a quantity type for generated fetches.  Works for
every field-less leaf (including parametric ones like `Energy{:per_site}`);
quantity types that REQUIRE constructor arguments (`RenyiEntropy(α)`, …) are
not instantiable from a bare type and raise an informative error — a
constraint edge over such a quantity must carry the instance itself.
"""
function _quantity_instance(::Type{Q}) where {Q<:AbstractQuantity}
    fieldcount(Q) == 0 || throw(
        ArgumentError(
            "_quantity_instance: $(Q) has fields; a constraint edge over it " *
            "must declare the quantity instance, not the bare type",
        ),
    )
    return Q()
end

"""
    _bc_instance(bc_T::Type{<:BoundaryCondition}; finite_N::Int) -> BoundaryCondition

Materialize a boundary condition from its registry type: `Infinite()` as-is,
`OBC`/`PBC` at the declared `finite_N` (constraint edges carry the size their
generated finite-N checks run at).
"""
function _bc_instance(::Type{BC}; finite_N::Int) where {BC<:BoundaryCondition}
    BC === Infinite && return Infinite()
    BC === OBC && return OBC(finite_N)
    BC === PBC && return PBC(; N=finite_N)
    throw(ArgumentError("_bc_instance: unsupported boundary-condition type $(BC)"))
end

"""
    _with_param(model, field::Symbol, val) -> model′

Reconstruct `model` with the named `field` replaced by `val` via the
positional constructor (the same mechanism as the identity harness's
`_perturb_field`).  Errors if the field is absent — a constraint edge naming a
non-existent parameter is a declaration bug, not a skip.
"""
function _with_param(model::AbstractQAtlasModel, field::Symbol, val)
    field in propertynames(model) || throw(
        ArgumentError(
            "_with_param: $(typeof(model)) has no field :$(field) " *
            "(fields: $(propertynames(model)))",
        ),
    )
    args = map(f -> f === field ? val : getfield(model, f), propertynames(model))
    return typeof(model).name.wrapper(args...)
end

# A registry row is *independent* (not delegation-backed) — the filter that
# keeps generated cross-checks non-circular: comparing a delegating row
# against its own delegation target verifies nothing.
_is_independent_row(e::Implementation) = !_is_delegation(e.method)

"""
    _canonical_row(model_T, quantity_T, bc_T) -> Union{Implementation,Nothing}

The canonical `REGISTRY` row of an exact `(model, quantity, bc)` hub, or
`nothing` — the shared lookup behind the duality/limit coherence checks and
generators (one definition instead of four copies of the scan loop).
"""
function _canonical_row(model_T::Type, quantity_T::Type, bc_T::Type)
    for e in REGISTRY
        if e.model === model_T && e.quantity === quantity_T && e.bc === bc_T && e.canonical
            return e
        end
    end
    return nothing
end

"""
    _both_endpoints_independent(source_T, target_T, quantity_T, bc_T) -> Bool

Both model endpoints have a canonical, independent (non-delegating) registry
row for `(quantity_T, bc_T)` — the precondition shared by the duality (#699)
and limit (#701) generators for emitting a genuine two-implementation
cross-check.
"""
function _both_endpoints_independent(source_T, target_T, quantity_T, bc_T)
    for m_T in (source_T, target_T)
        row = _canonical_row(m_T, quantity_T, bc_T)
        (row !== nothing && _is_independent_row(row)) || return false
    end
    return true
end

"""
    _check_endpoint_rows!(out, source_T, target_T, quantity_T, bc_T, tag, label)

Shared C12/C13 endpoint-row coherence: append a `:gap` finding (tagged `tag`,
prefixed `label`) for each of the two model endpoints whose `(quantity, bc)`
row is missing or delegation-backed — the cross-check the edge promises
cannot be generated, or would be circular.
"""
function _check_endpoint_rows!(
    out, source_T, target_T, quantity_T, bc_T, tag::Symbol, label::AbstractString
)
    for (side, m_T) in ((:source, source_T), (:target, target_T))
        row = _canonical_row(m_T, quantity_T, bc_T)
        if row === nothing
            push!(
                out,
                CoherenceFinding(
                    tag,
                    :gap,
                    "$(label) lists $(_kgshort(quantity_T)) at $(_kgshort(bc_T)) but " *
                    "the $(side) side ($(_kgshort(m_T))) has no canonical row — the " *
                    "promised cross-check cannot be generated",
                ),
            )
        elseif !_is_independent_row(row)
            push!(
                out,
                CoherenceFinding(
                    tag,
                    :gap,
                    "$(label): the $(side) row for $(_kgshort(quantity_T)) at " *
                    "$(_kgshort(bc_T)) is delegation-backed — the cross-check would " *
                    "be circular and is not generated",
                ),
            )
        end
    end
    return out
end

"""
    _equivalent_rows(rowq::Type, Q::Type) -> Bool

Extension point of [`_row_covers`](@ref): declare that a registered quantity
type covers a *different* requested type because `fetch` routes between them
automatically.  New equivalence axes add a method right below this fallback
(co-located with the quantity that owns the routing would be cleaner but the
kernel is the single import point all generators see), NOT an edit to the
hub-enumeration loop.
"""
_equivalent_rows(::Type, ::Type) = false
# Energy granularity: per-site/total conversion is automatic by design and
# conversion fallbacks are deliberately NOT registered (see the
# core/registry.jl header), so an `Energy{:total}` row covers an
# `Energy{:per_site}` request and vice versa.
_equivalent_rows(::Type{<:Energy}, ::Type{<:Energy}) = true

"""
    _row_covers(rowq::Type, Q::Type) -> Bool

Whether a registry row carrying quantity `rowq` covers a request for `Q`:
exact match, or any declared [`_equivalent_rows`](@ref) routing equivalence.
"""
_row_covers(rowq::Type, Q::Type) = rowq === Q || _equivalent_rows(rowq, Q)

"""
    _implemented_hubs(quantities; require_independent=false) -> Vector{NamedTuple}

The `(model, bc)` pairs whose canonical registry rows cover ALL of
`quantities` (exact quantity types, modulo the Energy-granularity
equivalence of `_row_covers`) — the hubs a constraint generator can
materialize checks on.  With `require_independent=true`, hubs where ANY of the
participating rows is delegation-backed are dropped (the #699/#701
circularity rule).  Universality / Bound namespaces are excluded: constraint
checks compare concrete-model implementations.

Deterministic: sorted by `(model name, bc name)`, one entry per hub.
"""
function _implemented_hubs(quantities; require_independent::Bool=false)
    qs = collect(Type, quantities)
    isempty(qs) && return NamedTuple[]
    hubs = NamedTuple[]
    seen = Set{Tuple{Type,Type}}()
    for seed in REGISTRY
        (seed.canonical && _row_covers(seed.quantity, first(qs))) || continue
        (_is_universality(seed.model) || _is_bound(seed.model)) && continue
        m_T, bc_T = seed.model, seed.bc
        (m_T, bc_T) in seen && continue
        push!(seen, (m_T, bc_T))
        covered = true
        independent = true
        for q in qs
            row = nothing
            for e in REGISTRY
                if e.model === m_T &&
                    e.bc === bc_T &&
                    e.canonical &&
                    _row_covers(e.quantity, q)
                    row = e
                    break
                end
            end
            row === nothing && (covered=false; break)
            _is_independent_row(row) || (independent = false)
        end
        covered || continue
        (require_independent && !independent) && continue
        push!(hubs, (model=m_T, bc=bc_T))
    end
    sort!(hubs; by=h -> (_kgshort(h.model), _kgshort(h.bc)))
    return hubs
end

# Cartesian sweep over a NamedTuple of value vectors:
# `(beta=[0.5, 1.0],)` → [(beta=0.5,), (beta=1.0,)].  Deterministic order.
function _sweep_points(sweep::NamedTuple)
    isempty(keys(sweep)) && return NamedTuple[NamedTuple()]
    pts = NamedTuple[NamedTuple()]
    for k in keys(sweep)
        vals = getfield(sweep, k)
        pts = [merge(p, NamedTuple{(k,)}((v,))) for p in pts for v in vals]
    end
    return pts
end

# Compact deterministic id fragment for a sweep point: "beta=0.5,q=3.14".
_point_id(p::NamedTuple) = join([string(k, "=", getfield(p, k)) for k in keys(p)], ",")

# The "/<point>" id suffix for a sweep point — empty for the no-sweep point.
# Collapses the `isempty(keys(p)) ? "" : "/" * _point_id(p)` idiom repeated by
# every generator.
_point_suffix(p::NamedTuple) = isempty(keys(p)) ? "" : "/" * _point_id(p)

# Emit a visible :skip GeneratedCheck for an excluded hub (declared exclusions
# never silently vanish).  Shared by the identity generators.
function _push_excluded_check!(
    out, kind::Symbol, id::AbstractString, reason::AbstractString
)
    push!(out, GeneratedCheck(kind, id, "EXCLUDED: $(reason)", () -> _skip_outcome(reason)))
    return out
end
