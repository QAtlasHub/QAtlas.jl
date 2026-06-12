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
defaults to `row -> row.references` (every current store carries a
`references` field); `location_of` should pin the row precisely enough that a
dangling-bibkey finding is actionable.
"""
function register_edge_store!(
    name::Symbol,
    store::Vector;
    references_of::Function=row -> row.references,
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

The result of running one [`GeneratedCheck`](@ref): `status` is `:pass` /
`:fail` / `:skip`, `lhs`/`rhs` are the two compared values (`NaN` when
skipped), and `detail` carries the skip reason or failure context.
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

"""
    run_generated_check(c::GeneratedCheck) -> CheckOutcome

Run `c`, converting a thrown exception into a `:fail` outcome whose `detail`
carries the error — generated suites report every check rather than aborting
at the first throwing hub.
"""
function run_generated_check(c::GeneratedCheck)
    try
        return c.run()
    catch err
        return CheckOutcome(:fail, NaN, NaN, NaN, NaN, string(typeof(err), ": ", err))
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
        append!(out, gen()::Vector{GeneratedCheck})
    end
    sort!(out; by=c -> c.id)
    allunique(c.id for c in out) ||
        error("generated_checks: duplicate check ids — a generator is not deterministic")
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
    _equivalent_rows(rowq::Type, Q::Type) -> Bool

Extension point of [`_row_covers`](@ref): declare that a registered quantity
type covers a *different* requested type because `fetch` routes between them
automatically.  New equivalence axes add a method HERE (next to the quantity
that owns the routing), not an edit to the kernel's hub enumeration.
"""
_equivalent_rows(::Type, ::Type) = false
# Energy granularity: per-site/total conversion is automatic by design and
# conversion fallbacks are deliberately NOT registered (see the
# core/registry.jl header), so an `Energy{:total}` row covers an
# `Energy{:per_site}` request and vice versa.
_equivalent_rows(::Type{<:Energy}, ::Type{<:Energy}) = true

# Does a registry row's quantity cover the requested one?  Exact match, plus
# any declared `_equivalent_rows` routing equivalence.
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
