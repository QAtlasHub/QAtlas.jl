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
# Two modes, one store:
#
#   * `:tuple` — an explicit relation over named quantities:
#         @identity(:gibbs,
#             quantities = (f=FreeEnergy, e=Energy{:per_site}, s=ThermalEntropy),
#             check = (v, p) -> (v.f, v.e - v.s / p.beta),
#             sweep = (beta = [0.5, 1.0, 2.0],))
#     `check` receives the fetched values and the sweep point as NamedTuples
#     and returns the `(lhs, rhs)` pair that must agree — returning both sides
#     (rather than a Bool) keeps failures diagnosable.
#
#   * `:component_isotropy` — the #690 × #700 integration: a quantity FAMILY
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
    IdentityEdge

One quantity↔quantity identity — see [`@identity`](@ref).  `mode` is `:tuple`
(explicit `quantities` + `check`) or `:component_isotropy` (`family` +
`requires_internal`).  `sweep` is the fetch-kwargs grid the generated checks
run on; `finite_N` the size used for `OBC`/`PBC` hubs; `exclusions` lists
`Model => reason` pairs that are emitted as visible `:skip` checks rather
than silently dropped.
"""
struct IdentityEdge
    name::Symbol
    mode::Symbol
    quantities::NamedTuple                 # :tuple mode — name => quantity Type
    check::Union{Function,Nothing}         # :tuple mode — (vals, point) -> (lhs, rhs)
    family::Union{Type,Nothing}            # :component_isotropy mode
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
    IDENTITIES :: Vector{IdentityEdge}

The quantity↔quantity identity store, populated at include-time by
[`identity!`](@ref) / [`@identity`](@ref).  Query with
[`identities_for`](@ref) / [`participants`](@ref); the generated checks are
the `:identity` kind of [`generated_checks`](@ref).
"""
const IDENTITIES = IdentityEdge[]

"""
    identity!(name; quantities=nothing, check=nothing, family=nothing,
              requires_internal=nothing, sweep=(;), finite_N=8,
              rtol=1e-8, atol=1e-10, exclusions=[], notes="", references=String[])

Record an identity edge.  Exactly one of the two mode signatures must be
given: (`quantities` + `check`) for a `:tuple` identity, or `family` (an
abstract quantity supertype, optionally `requires_internal` symmetry-gated)
for `:component_isotropy`.  Tuple participants must be field-less quantity
types (instantiable from the bare type).
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
    push!(
        IDENTITIES,
        IdentityEdge(
            name,
            tuple_mode ? :tuple : :component_isotropy,
            tuple_mode ? quantities : NamedTuple(),
            check,
            family,
            requires_internal,
            sweep,
            finite_N,
            Float64(rtol),
            Float64(atol),
            excl,
            String(notes),
            String[r for r in references],
        ),
    )
    return nothing
end

# The declared skip reason for a hub, if any: an exclusion keyed by the model
# type covers every bc; a (model, bc) tuple covers that hub only.
function _exclusion_reason(e::IdentityEdge, model_T::Type, bc_T::Type)
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

"""
    identities_for(quantity) -> Vector{IdentityEdge}

The identity edges `quantity` (instance or type) participates in — tuple
identities naming it, and family identities whose family it belongs to.
"""
function identities_for(quantity)
    q_T = _as_type(quantity)
    out = IdentityEdge[]
    for e in IDENTITIES
        if e.mode === :tuple
            any(Q -> Q === q_T, values(e.quantities)) && push!(out, e)
        else
            q_T <: e.family && component(q_T) !== nothing && push!(out, e)
        end
    end
    return out
end

"""
    participants(edge::IdentityEdge) -> Vector{Type}

The quantity types `edge` relates: the declared tuple, or the family's
component-carrying concrete members.
"""
function participants(edge::IdentityEdge)
    edge.mode === :tuple && return Type[Q for Q in values(edge.quantities)]
    return Type[T for T in _family_members(edge.family)]
end

# ──────────────────────────────────────────────────────────────────────
# C11 — static identity coherence (no fetch execution)
# ──────────────────────────────────────────────────────────────────────

"""
    check_identity_coverage() -> Vector{CoherenceFinding}

C11: every identity edge should be *exercised* — an edge none of whose
participant sets is implemented by any hub generates zero checks, i.e. the
declared relation constrains nothing.  Self-reported as a `:gap` (a
missing-but-expected hub, not an invariant violation).  A symmetry-gated
family identity whose `requires_internal` matches no [`@symmetry`](@ref)
profile is likewise a `:gap` (the gate is closed for want of profiles).
"""
function check_identity_coverage()
    out = CoherenceFinding[]
    for e in IDENTITIES
        if e.mode === :tuple
            isempty(_implemented_hubs(values(e.quantities))) && push!(
                out,
                CoherenceFinding(
                    :identity_coverage,
                    :gap,
                    "identity :$(e.name) has no (model, bc) hub implementing all of " *
                    "$(join(_kgshort.(collect(values(e.quantities))), ", ")) — it " *
                    "generates no checks",
                ),
            )
        elseif e.requires_internal !== nothing
            isempty(models_with_symmetry(e.requires_internal)) && push!(
                out,
                CoherenceFinding(
                    :identity_coverage,
                    :gap,
                    "identity :$(e.name) requires internal :$(e.requires_internal) " *
                    "but no @symmetry profile declares it — the gate is closed",
                ),
            )
        end
    end
    return out
end

# ──────────────────────────────────────────────────────────────────────
# Generator — the :identity kind of generated_checks()
# ──────────────────────────────────────────────────────────────────────

function _identity_tuple_checks(e::IdentityEdge)
    out = GeneratedCheck[]
    for hub in _implemented_hubs(values(e.quantities))
        hub_id = string(
            "identity/", e.name, "/", _kgshort(hub.model), "/", _kgshort(hub.bc)
        )
        reason = _exclusion_reason(e, hub.model, hub.bc)
        if reason !== nothing
            push!(
                out,
                GeneratedCheck(
                    :identity, hub_id, "EXCLUDED: $(reason)", () -> _skip_outcome(reason)
                ),
            )
            continue
        end
        for point in _sweep_points(e.sweep)
            id = isempty(keys(point)) ? hub_id : string(hub_id, "/", _point_id(point))
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

function _identity_isotropy_checks(e::IdentityEdge)
    out = GeneratedCheck[]
    members = _family_members(e.family)
    gated_models = if e.requires_internal === nothing
        nothing
    else
        Set(models_with_symmetry(e.requires_internal))
    end
    # one pass: (model, bc) hub => family members it implements
    impl = Dict{Tuple{Type,Type},Vector{Type}}()
    for Q in members, hub in _implemented_hubs((Q,))
        push!(get!(impl, (hub.model, hub.bc), Type[]), Q)
    end
    for key in sort!(collect(keys(impl)); by=k -> (_kgshort(k[1]), _kgshort(k[2])))
        model_T, bc_T = key
        gated_models === nothing || model_T in gated_models || continue
        present = sort!(impl[key]; by=string)
        # One representative per component: isotropy equates DIFFERENT
        # components, so same-component family variants (e.g. two modes of an
        # :xx correlator) are not an isotropy statement — and including them
        # would collide the component-keyed check ids.
        comps = Set{Symbol}()
        reps = Type[]
        for Q in present
            c = component(Q)::Symbol
            c in comps || (push!(reps, Q); push!(comps, c))
        end
        length(reps) ≥ 2 || continue
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
                isempty(keys(point)) ? "" : "/" * _point_id(point),
            )
            reason = _exclusion_reason(e, model_T, bc_T)
            if reason !== nothing
                push!(
                    out,
                    GeneratedCheck(
                        :identity, id, "EXCLUDED: $(reason)", () -> _skip_outcome(reason)
                    ),
                )
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
        append!(
            out,
            e.mode === :tuple ? _identity_tuple_checks(e) : _identity_isotropy_checks(e),
        )
    end
    return out
end

register_check_generator!(:identity, identity_checks)
register_edge_store!(:identity, IDENTITIES; location_of=e -> "identity :$(e.name)")
