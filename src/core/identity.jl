# core/identity.jl — quantity ↔ quantity identity edges (#698, over the #697
# kernel).
#
# An identity edge declares that N quantities are NOT independent: they are
# bound by an exact relation (Gibbs f = ε − T·s, sum rules, …) that every
# implementation must satisfy.  The edge is the single source of truth; the
# kernel's generator protocol turns it into executable checks for every
# (model, bc) hub that implements all participants — so a NEW model gets
# identity coverage with zero hand-written tests, and a change to any
# participating quantity has a declared blast radius (selective CI).
#
# Two edge types under `AbstractIdentityEdge`, one store (queries and
# generators dispatch on the type — no `mode` tag, no nothing-filled halves):
#
#   * `TupleIdentityEdge` — an explicit relation over named quantities:
#         @identity(:gibbs,
#             quantities = (f=FreeEnergy, e=Energy{:per_site}, s=ThermalEntropy),
#             check = (v, p) -> (v.f, v.e - v.s / p.beta),
#             sweep = (beta = [0.5, 1.0, 2.0],))
#     `check` receives the fetched values and the sweep point as NamedTuples
#     and returns the `(lhs, rhs)` pair that must agree — returning both sides
#     (rather than a Bool) keeps failures diagnosable.
#
#   * `IsotropyIdentityEdge` — the #690 × #700 integration: a quantity FAMILY
#     (an abstract supertype from the taxonomy layer) whose components must
#     coincide for models declaring a symmetry:
#         @identity(:su2_susceptibility_isotropy,
#             family = AbstractSusceptibility,
#             requires_internal = :SU2,
#             sweep = (beta = [0.5, 1.0],))
#     generates pairwise `component` equalities (χ_xx = χ_yy = χ_zz) for every
#     model whose `@symmetry` profile declares `internal=:SU2` — the registry
#     replaces the harness's hand-coded `is_su2_symmetric` model filter.
#
# Identities are *internal consistency* relations, so delegation-backed rows
# are legitimate participants (the relation must hold whatever the route);
# the independence filter of duality/limit edges does not apply here.

"""
    AbstractIdentityEdge

A quantity↔quantity identity — one of the two concrete modes below.  Splitting
the two modes into distinct types (rather than one struct with a `mode` tag and
half its fields left `nothing`) makes the inactive-half states unrepresentable
and lets the queries/generators dispatch instead of branch.  Shared fields,
present on both: `name`, `sweep` (fetch-kwargs grid), `finite_N` (OBC/PBC hub
size), `rtol`/`atol`, `exclusions` (`Model`/`(Model,BC) => reason` pairs
emitted as visible `:skip` checks), `notes`, `references`.
"""
abstract type AbstractIdentityEdge end

"""
    TupleIdentityEdge <: AbstractIdentityEdge

An explicit relation over named quantities: `check(vals, point) -> (lhs, rhs)`
must agree for every hub implementing all of `quantities` (a NamedTuple
`name => quantity Type`).  See [`@identity`](@ref).
"""
struct TupleIdentityEdge <: AbstractIdentityEdge
    name::Symbol
    quantities::NamedTuple
    check::Function
    sweep::NamedTuple
    finite_N::Int
    rtol::Float64
    atol::Float64
    exclusions::Vector{Pair{Any,String}}
    notes::String
    references::Vector{String}
end

"""
    IsotropyIdentityEdge <: AbstractIdentityEdge

A component-isotropy relation over a quantity `family` (an abstract taxonomy
supertype): the family's components must coincide, optionally gated on a
`requires_internal` symmetry.  See [`@identity`](@ref).
"""
struct IsotropyIdentityEdge <: AbstractIdentityEdge
    name::Symbol
    family::Type
    requires_internal::Union{Symbol,Nothing}
    sweep::NamedTuple
    finite_N::Int
    rtol::Float64
    atol::Float64
    exclusions::Vector{Pair{Any,String}}
    notes::String
    references::Vector{String}
end

"""
    IDENTITIES :: Vector{AbstractIdentityEdge}

The quantity↔quantity identity store, populated at include-time by
[`identity!`](@ref) / [`@identity`](@ref).  Query with
[`identities_for`](@ref) / [`participants`](@ref); the generated checks are
the `:identity` kind of [`generated_checks`](@ref).
"""
const IDENTITIES = AbstractIdentityEdge[]

"""
    identity!(name; quantities=nothing, check=nothing, family=nothing,
              requires_internal=nothing, sweep=(;), finite_N=8,
              rtol=1e-8, atol=1e-10, exclusions=[], notes="", references=String[])

Record an identity edge (a [`TupleIdentityEdge`](@ref) or
[`IsotropyIdentityEdge`](@ref)).  Exactly one of the two mode signatures must
be given: (`quantities` + `check`) for a tuple identity, or `family` (an
abstract quantity supertype, optionally `requires_internal` symmetry-gated)
for component isotropy.  Tuple participants must be field-less quantity types
(instantiable from the bare type).

`finite_N` is the OBC/PBC hub size the generated checks run at; the default 8
suits spin-1/2 hubs but is a 3ᴺ dense-ED trap for spin-1 models — prefer 6 for
families that reach S=1 chains (see `src/identity_registry.jl`).  `rtol` must
be in `[0, 1)` (a value ≥ 1 would pass every check) and `atol ≥ 0`.
"""
function identity!(
    name::Symbol;
    quantities::Union{NamedTuple,Nothing}=nothing,
    check::Union{Function,Nothing}=nothing,
    family::Union{Type,Nothing}=nothing,
    requires_internal::Union{Symbol,Nothing}=nothing,
    sweep::NamedTuple=NamedTuple(),
    finite_N::Int=8,
    rtol::Real=1e-8,
    atol::Real=1e-10,
    exclusions::AbstractVector=Pair{Any,String}[],
    notes::AbstractString="",
    references::AbstractVector{<:AbstractString}=String[],
)
    any(e -> e.name === name, IDENTITIES) &&
        throw(ArgumentError("identity!: :$(name) already declared"))
    _check_tolerances(:identity, rtol, atol)
    tuple_mode = quantities !== nothing
    family_mode = family !== nothing
    tuple_mode ⊻ family_mode || throw(
        ArgumentError(
            "identity!: declare EITHER quantities+check (tuple identity) OR " *
            "family (component isotropy), not " *
            (tuple_mode ? "both" : "neither"),
        ),
    )
    if tuple_mode
        check isa Function ||
            throw(ArgumentError("identity!: a tuple identity requires a check function"))
        length(quantities) ≥ 2 ||
            throw(ArgumentError("identity!: an identity relates ≥ 2 quantities"))
        for (k, Q) in pairs(quantities)
            Q isa Type && Q <: AbstractQuantity || throw(
                ArgumentError("identity!: participant $(k)=$(Q) is not a quantity type")
            )
            _quantity_instance(Q)   # throws informatively for field-carrying types
        end
        requires_internal === nothing || throw(
            ArgumentError(
                "identity!: requires_internal is the symmetry gate of the " *
                "component-isotropy mode; a tuple identity holds unconditionally",
            ),
        )
    else
        check === nothing ||
            throw(ArgumentError("identity!: family mode derives its checks; drop `check`"))
        (family isa Type && isabstracttype(family) && family <: AbstractQuantity) || throw(
            ArgumentError(
                "identity!: family must be an abstract quantity supertype " *
                "(taxonomy layer); got $(family)",
            ),
        )
        length(_family_members(family)) ≥ 2 || throw(
            ArgumentError(
                "identity!: family $(family) has < 2 concrete members with a " *
                "declared component — nothing to equate",
            ),
        )
    end
    excl = Pair{Any,String}[]
    for p in exclusions
        k = first(p)
        k isa Type ||
            (k isa Tuple && length(k) == 2 && k[1] isa Type && k[2] isa Type) ||
            throw(
                ArgumentError(
                    "identity!: exclusion keys are a Model type (all bcs) or a " *
                    "(Model, BC) tuple; got $(repr(k))",
                ),
            )
        push!(excl, Pair{Any,String}(k, String(last(p))))
    end
    edge = if tuple_mode
        TupleIdentityEdge(
            name,
            quantities,
            check,
            sweep,
            finite_N,
            Float64(rtol),
            Float64(atol),
            excl,
            String(notes),
            String[r for r in references],
        )
    else
        IsotropyIdentityEdge(
            name,
            family,
            requires_internal,
            sweep,
            finite_N,
            Float64(rtol),
            Float64(atol),
            excl,
            String(notes),
            String[r for r in references],
        )
    end
    push!(IDENTITIES, edge)
    return nothing
end

# The declared skip reason for a hub, if any: an exclusion keyed by the model
# type covers every bc; a (model, bc) tuple covers that hub only.
# Untyped in `e`: this reads only `.exclusions`, so it serves the bound edges of
# core/bound.jl as well — the exclusion vocabulary is shared, not identity-specific.
function _exclusion_reason(e, model_T::Type, bc_T::Type)
    for (k, reason) in e.exclusions
        if k === model_T || (k isa Tuple && k[1] === model_T && k[2] === bc_T)
            return reason
        end
    end
    return nothing
end

"""
    @identity :name key=value …

Macro sugar around [`identity!`](@ref): the positional `:name` is the edge's
identifier; the remaining `key=value` pairs are forwarded as keyword
arguments.  See `src/identity_registry.jl` for the declared catalog.
"""
macro identity(name, kwargs...)
    return _forward_kw_macro(identity!, :identity, (name,), kwargs)
end

# Concrete family members that participate in component identities: the
# REGISTERED quantity types under the family supertype that carry a
# `component` (the `…Local` variants deliberately do not — their extra site
# argument is a different fetch shape).  Registry-driven rather than
# reflection-driven (`subtypes` would need InteractiveUtils and would see
# leaves no model implements anyway).  Deterministic order.
function _family_members(family::Type)
    members = Type[]
    for e in REGISTRY
        e.quantity <: family || continue
        component(e.quantity) === nothing && continue
        e.quantity in members || push!(members, e.quantity)
    end
    sort!(members; by=string)
    return members
end

_participates(e::TupleIdentityEdge, q_T::Type) = any(Q -> Q === q_T, values(e.quantities))
function _participates(e::IsotropyIdentityEdge, q_T::Type)
    return q_T <: e.family && component(q_T) !== nothing
end

"""
    identities_for(quantity) -> Vector{AbstractIdentityEdge}

The identity edges `quantity` (instance or type) participates in — tuple
identities naming it, and family identities whose family it belongs to.
"""
function identities_for(quantity)
    q_T = _as_type(quantity)
    return AbstractIdentityEdge[e for e in IDENTITIES if _participates(e, q_T)]
end

"""
    participants(edge::AbstractIdentityEdge) -> Vector{Type}

The quantity types `edge` relates: the declared tuple, or the family's
component-carrying concrete members.
"""
participants(e::TupleIdentityEdge) = Type[Q for Q in values(e.quantities)]
participants(e::IsotropyIdentityEdge) = Type[T for T in _family_members(e.family)]

# ──────────────────────────────────────────────────────────────────────
# C11 — static identity coherence (no fetch execution)
# ──────────────────────────────────────────────────────────────────────

"""
    check_identity_coverage() -> Vector{CoherenceFinding}

C11: every identity edge should be *exercised* — an edge that generates zero
checks constrains nothing, a self-reported `:gap`.  Dispatches per edge type:
a tuple identity no hub implements; a gated isotropy identity whose
`requires_internal` matches no [`@symmetry`](@ref) profile (gate closed); or
ANY isotropy identity — gated or not — with no hub implementing ≥ 2 distinct
family components (the ungated case the modal version silently skipped).
"""
function check_identity_coverage()
    out = CoherenceFinding[]
    for e in IDENTITIES
        f = _coverage_finding(e)
        f === nothing || push!(out, f)
    end
    return out
end

function _coverage_finding(e::TupleIdentityEdge)
    isempty(_implemented_hubs(values(e.quantities))) || return nothing
    return CoherenceFinding(
        :identity_coverage,
        :gap,
        "identity :$(e.name) has no (model, bc) hub implementing all of " *
        "$(join(_kgshort.(collect(values(e.quantities))), ", ")) — it generates no checks",
    )
end

function _coverage_finding(e::IsotropyIdentityEdge)
    if e.requires_internal !== nothing && isempty(models_with_symmetry(e.requires_internal))
        return CoherenceFinding(
            :identity_coverage,
            :gap,
            "identity :$(e.name) requires internal :$(e.requires_internal) but no " *
            "@symmetry profile declares it — the gate is closed",
        )
    end
    # Generates a check only where some (gated) hub implements ≥ 2 distinct
    # components; if none does, the edge is inert — covers the ungated case too.
    isempty(_isotropy_hubs(e)) || return nothing
    gate = if e.requires_internal === nothing
        ""
    else
        " (gated on internal :$(e.requires_internal))"
    end
    return CoherenceFinding(
        :identity_coverage,
        :gap,
        "identity :$(e.name) over $(_kgshort(e.family))$(gate) has no hub implementing " *
        "≥ 2 distinct components — it generates no checks",
    )
end

# ──────────────────────────────────────────────────────────────────────
# Generator — the :identity kind of generated_checks()
# ──────────────────────────────────────────────────────────────────────

function _identity_tuple_checks(e::TupleIdentityEdge)
    out = GeneratedCheck[]
    for hub in _implemented_hubs(values(e.quantities))
        hub_id = string(
            "identity/", e.name, "/", _kgshort(hub.model), "/", _kgshort(hub.bc)
        )
        reason = _exclusion_reason(e, hub.model, hub.bc)
        if reason !== nothing
            _push_excluded_check!(out, :identity, hub_id, reason)
            continue
        end
        for point in _sweep_points(e.sweep)
            id = hub_id * _point_suffix(point)
            names = keys(e.quantities)
            qtypes = values(e.quantities)
            model_T, bc_T = hub.model, hub.bc
            runner = function ()
                m = model_T()
                bc = _bc_instance(bc_T; finite_N=e.finite_N)
                vals = NamedTuple{names}(
                    map(Q -> fetch(m, _quantity_instance(Q), bc; point...), qtypes)
                )
                lhs, rhs = e.check(vals, point)
                return _outcome(lhs, rhs; rtol=e.rtol, atol=e.atol)
            end
            push!(
                out,
                GeneratedCheck(
                    :identity,
                    id,
                    string(
                        "identity :",
                        e.name,
                        " on ",
                        _kgshort(model_T),
                        " at ",
                        _kgshort(bc_T),
                        isempty(keys(point)) ? "" : " ($(point))",
                    ),
                    runner,
                ),
            )
        end
    end
    return out
end

# The hubs an isotropy edge generates checks on: `(model, bc, reps)` where
# `reps` is one representative quantity per distinct component (≥ 2 of them),
# restricted to the symmetry gate when present.  Shared by the generator and
# the C11 coverage check so the two cannot disagree on "does this generate?".
# One representative per component because isotropy equates DIFFERENT
# components; same-component family variants (e.g. two modes of an :xx
# correlator) are not an isotropy statement and would collide the
# component-keyed ids.  Deterministic order.
function _isotropy_hubs(e::IsotropyIdentityEdge)
    members = _family_members(e.family)
    gated = if e.requires_internal === nothing
        nothing
    else
        Set(models_with_symmetry(e.requires_internal))
    end
    impl = Dict{Tuple{Type,Type},Vector{Type}}()
    for Q in members, hub in _implemented_hubs((Q,))
        push!(get!(impl, (hub.model, hub.bc), Type[]), Q)
    end
    out = Tuple{Type,Type,Vector{Type}}[]
    for key in sort!(collect(keys(impl)); by=k -> (_kgshort(k[1]), _kgshort(k[2])))
        model_T, bc_T = key
        gated === nothing || model_T in gated || continue
        comps = Set{Symbol}()
        reps = Type[]
        for Q in sort!(impl[key]; by=string)
            c = component(Q)::Symbol
            c in comps || (push!(reps, Q); push!(comps, c))
        end
        length(reps) ≥ 2 && push!(out, (model_T, bc_T, reps))
    end
    return out
end

function _identity_isotropy_checks(e::IsotropyIdentityEdge)
    out = GeneratedCheck[]
    for (model_T, bc_T, reps) in _isotropy_hubs(e)
        ref_Q = reps[1]
        for other_Q in reps[2:end], point in _sweep_points(e.sweep)
            pair = string(component(ref_Q), "=", component(other_Q))
            id = string(
                "identity/",
                e.name,
                "/",
                _kgshort(model_T),
                "/",
                _kgshort(bc_T),
                "/",
                pair,
                _point_suffix(point),
            )
            reason = _exclusion_reason(e, model_T, bc_T)
            if reason !== nothing
                _push_excluded_check!(out, :identity, id, reason)
                continue
            end
            runner = function ()
                m = model_T()
                bc = _bc_instance(bc_T; finite_N=e.finite_N)
                a = fetch(m, _quantity_instance(ref_Q), bc; point...)
                b = fetch(m, _quantity_instance(other_Q), bc; point...)
                return _outcome(a, b; rtol=e.rtol, atol=e.atol)
            end
            push!(
                out,
                GeneratedCheck(
                    :identity,
                    id,
                    string(
                        "isotropy :",
                        e.name,
                        " ",
                        pair,
                        " on ",
                        _kgshort(model_T),
                        " at ",
                        _kgshort(bc_T),
                    ),
                    runner,
                ),
            )
        end
    end
    return out
end

function identity_checks()
    out = GeneratedCheck[]
    for e in IDENTITIES
        append!(out, _edge_checks(e))
    end
    return out
end

_edge_checks(e::TupleIdentityEdge) = _identity_tuple_checks(e)
_edge_checks(e::IsotropyIdentityEdge) = _identity_isotropy_checks(e)

register_check_generator!(:identity, identity_checks)
register_edge_store!(:identity, IDENTITIES; location_of=e -> "identity :$(e.name)")
