# core/bound.jl — quantity BOUND edges: the inequality sibling of core/identity.jl.
#
# An identity edge says N quantities must satisfy an exact equality.  A bound
# edge says a quantity (or a tuple of them) must satisfy an INEQUALITY that no
# implementation may violate: `C_v ≥ 0`, `χ_T ≥ 0`, entropy non-negativity, …
# Physically these are stability and information-theoretic bounds, and they are
# exactly the checks that catch a sign error, a bad analytic continuation, or a
# thermal-state normalization mistake — failures an equality identity between
# two equally-wrong quantities can miss.
#
# WHERE THE MATH LIVES.  Nowhere here.  AbstractQAtlas already states these as
# `@inequality` relations with a `slack` verb whose sign IS the criterion, so a
# bound edge names the AbstractQAtlas inequality and QAtlas supplies only the
# values and the harness — the same split #734 Phase B established for `:gibbs`,
# where the arithmetic moved to `FreeEnergyLegendre` and the sweep / finite_N /
# exclusions stayed here.  The participants are DERIVED from the inequality's
# typed slots (`variable_slots`), so a bound edge cannot drift from the relation
# it claims to check: there is no second copy of the participant list to update.
#
# FAMILY SLOTS.  AbstractQAtlas types some slots as a parametric FAMILY rather
# than a concrete leaf (`SusceptibilityPositivity(χT::Susceptibility)` covers
# every axis pair).  Its own bag matching cannot bind those (its C4 note calls
# such a relation permanently dead there), but a generator that fetches by
# concrete type can: a family slot expands to one check per concrete member the
# hub actually implements, the same way `IsotropyIdentityEdge` walks a family.
# That is why these bounds are usable from here even though they are inert
# upstream.

"""
    BoundEdge

A one-sided constraint on a hub's fetched values: `slack(inequality; …) ≥ 0`.

Fields mirror [`TupleIdentityEdge`](@ref) — `name`, `sweep` (fetch-kwargs grid),
`finite_N` (OBC/PBC hub size), `atol`, `exclusions` (`Model` / `(Model, BC) =>
reason`, emitted as visible `:skip` checks), `notes`, `references` — plus:

- `inequality` — the AbstractQAtlas `AbstractInequality` instance that owns the
  statement and computes the slack.
- `quantities` — slot name => quantity type, **derived** from the inequality's
  typed slots, not restated by the caller.

Only `atol` is meaningful (no `rtol`): the criterion is `slack ≥ -atol`, a
one-sided test against zero, where a relative tolerance has no reference scale.
"""
struct BoundEdge
    name::Symbol
    inequality::AbstractInequality
    quantities::NamedTuple
    sweep::NamedTuple
    finite_N::Int
    atol::Float64
    exclusions::Vector{Pair{Any,String}}
    notes::String
    references::Vector{String}
end

"""
    BOUNDS :: Vector{BoundEdge}

Every declared bound edge, in declaration order.
"""
const BOUNDS = BoundEdge[]

"""
    bound!(name; inequality, sweep, finite_N, atol, exclusions, notes, references)

Declare a bound edge.  `inequality` is an AbstractQAtlas `AbstractInequality`
type or instance; its typed slots become the participants.

Rejected at declaration time, loudly, because each would otherwise degrade into
a silently-useless check:

- a duplicate `name`;
- a negative `atol`;
- an inequality with **no** typed slots (nothing to fetch — it could never be
  materialized on a hub);
- an inequality with any **untyped** slot (a derived input such as `var_E` or
  `S_AB` that the generator cannot obtain from `fetch`).  Those relations are
  real and worth covering, but they need a supplier, not a silent skip.
"""
function bound!(
    name::Symbol;
    inequality,
    sweep::NamedTuple=NamedTuple(),
    finite_N::Int=8,
    atol::Real=1e-10,
    exclusions::AbstractVector=Pair{Any,String}[],
    notes::AbstractString="",
    references::AbstractVector{<:AbstractString}=String[],
)
    any(e -> e.name === name, BOUNDS) &&
        throw(ArgumentError("bound!: :$(name) already declared"))
    atol ≥ 0 || throw(ArgumentError("bound!: atol must be ≥ 0; got $(atol)"))
    ineq = inequality isa Type ? inequality() : inequality
    ineq isa AbstractInequality || throw(
        ArgumentError(
            "bound!: :$(name) — `inequality` must be an AbstractQAtlas " *
            "AbstractInequality; got $(typeof(ineq))",
        ),
    )
    slots = variable_slots(ineq)
    # `variable_slots` yields (name, type-or-nothing) pairs directly — do NOT wrap
    # it in `pairs()`, which would re-key them by integer index.
    typed = [(n, T) for (n, T) in slots if T !== nothing]
    untyped = [n for (n, T) in slots if T === nothing]
    isempty(typed) && throw(
        ArgumentError(
            "bound!: :$(name) — $(nameof(typeof(ineq))) has no type-keyed slot, so " *
            "there is nothing to fetch; it cannot be materialized on a hub",
        ),
    )
    isempty(untyped) || throw(
        ArgumentError(
            "bound!: :$(name) — $(nameof(typeof(ineq))) needs the untyped slot(s) " *
            "$(untyped), which the generator cannot obtain from `fetch`. Supply " *
            "them via a dedicated edge rather than declaring a bound that would " *
            "never run.",
        ),
    )
    qs = NamedTuple{Tuple(first.(typed))}(Tuple(last.(typed)))
    push!(
        BOUNDS,
        BoundEdge(
            name,
            ineq,
            qs,
            sweep,
            finite_N,
            Float64(atol),
            Pair{Any,String}[p for p in exclusions],
            String(notes),
            String[r for r in references],
        ),
    )
    return nothing
end

"""
    @bound(:name, inequality = SomeInequality, sweep = (...,), ...)

Keyword-macro front end for [`bound!`](@ref), matching [`@identity`](@ref)'s
call shape so the two read alike at the declaration site.
"""
macro bound(name, kwargs...)
    kws = [esc(k) for k in kwargs]
    return :(bound!($(esc(name)); $(kws...)))
end

# ──────────────────────────────────────────────────────────────────────
# Generator — the :bound kind of generated_checks()
# ──────────────────────────────────────────────────────────────────────

"""
    _bound_outcome(slack; atol, detail="") -> CheckOutcome

One-sided pass criterion: `slack ≥ -atol`.  Reported as `lhs = slack`,
`rhs = 0`, so the recorded `abs_err` is the VIOLATION magnitude (zero for any
satisfied bound, however slack), not a distance from a target — a bound that
passes with room to spare is not "less accurate".
"""
function _bound_outcome(slack::Real; atol::Real, detail::String="")
    s = Float64(slack)
    status = s ≥ -atol ? :pass : :fail
    violation = status === :pass ? 0.0 : abs(s)
    return CheckOutcome(status, s, 0.0, violation, violation, detail)
end

# A family slot (`Susceptibility`) stands for every concrete member; a concrete
# slot stands for itself.  Returns the per-slot candidate lists, in a
# deterministic order, so the generated ids are stable across runs.
function _bound_slot_candidates(e::BoundEdge)
    return map(values(e.quantities)) do Q
        return isconcretetype(Q) ? Type[Q] : sort!(_family_members(Q); by=string)
    end
end

function bound_checks()
    out = GeneratedCheck[]
    for e in BOUNDS
        names = keys(e.quantities)
        for combo in Iterators.product(_bound_slot_candidates(e)...)
            qtypes = collect(combo)
            for hub in _implemented_hubs(qtypes)
                tag =
                    length(qtypes) == 1 ? _kgshort(qtypes[1]) : join(_kgshort.(qtypes), "-")
                hub_id = string(
                    "bound/",
                    e.name,
                    "/",
                    _kgshort(hub.model),
                    "/",
                    _kgshort(hub.bc),
                    "/",
                    tag,
                )
                reason = _exclusion_reason(e, hub.model, hub.bc)
                if reason !== nothing
                    _push_excluded_check!(out, :bound, hub_id, reason)
                    continue
                end
                for point in _sweep_points(e.sweep)
                    id = hub_id * _point_suffix(point)
                    model_T, bc_T, ineq = hub.model, hub.bc, e.inequality
                    runner = function ()
                        m = model_T()
                        bc = _bc_instance(bc_T; finite_N=e.finite_N)
                        vals = NamedTuple{names}(
                            map(Q -> fetch(m, _quantity_instance(Q), bc; point...), qtypes),
                        )
                        return _bound_outcome(slack(ineq; vals...); atol=e.atol)
                    end
                    push!(
                        out,
                        GeneratedCheck(
                            :bound,
                            id,
                            string(
                                "bound :",
                                e.name,
                                " (",
                                nameof(typeof(ineq)),
                                ") on ",
                                _kgshort(model_T),
                                " at ",
                                _kgshort(bc_T),
                                " for ",
                                tag,
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

register_check_generator!(:bound, bound_checks)
register_edge_store!(:bound, BOUNDS; location_of=e -> "bound :$(e.name)")
