# core/response.jl — constraint edges whose relation needs a DERIVED input.
#
# The third edge kind, after @identity (equalities between fetched values) and
# @bound (inequalities on them).  A response edge covers the AbstractQAtlas
# relations stated with a supplied derivative — `EntropyResponse(S, dF_dT)`,
# `SpecificHeatFromEntropy(C, dS_dT, T)` — which the other two cannot express
# because one of their slots is not fetchable.
#
# The split of labour is the same as everywhere else in this layer: the relation
# and its algebra are AbstractQAtlas's, the values and the harness are QAtlas's,
# and the derivative comes from core/derivative.jl (AD when an extension is
# loaded, finite differences otherwise).
#
# A relation has three kinds of slot, and this file is mostly about telling them
# apart:
#
#   * QUANTITY slots (`S::ThermalEntropy`) — fetched from the hub.
#   * FIELD slots (`T::Temperature`, `β::InverseTemperature`) — these are
#     AbstractField, NOT AbstractQuantity: they are the state point, so they come
#     from the sweep, not from `fetch`.
#   * UNTYPED slots (`dF_dT`) — the derived input this edge supplies.
#
# SUBJECT-VS-DERIVED, not residual-vs-zero.  Like `:gibbs` (#734 Phase B), the
# check `solve`s the relation for its quantity slot and compares that against
# the hub's own fetched value, so both sides of a failure are physical numbers.
# A bare `residual ≈ 0` would be shorter and much harder to read when it breaks.

"""
    DerivedInput(quantity, wrt)

A slot supplied by differentiating `fetch(model, quantity(), bc)` with respect
to `wrt`, which is `:T` (temperature) or `:β` (inverse temperature).

Two optional transforms cover the slots that are not simply `d⟨Q⟩/dx`:

- `of(value, x)` — what to differentiate, as a function of the fetched value and
  the axis variable.  `GibbsHelmholtz` wants `d(βF)/dβ`, not `dF/dβ`, so it
  passes `of = (F, β) -> β * F`.
- `then(derivative)` — post-processing.  `SpecificHeatFDT`'s energy variance is
  `var_E = -∂U/∂β`, so it passes `then = d -> -d`.

Both default to the identity, and both are deliberately plain functions rather
than a symbolic mini-language: the whole point of this layer is that the algebra
lives in AbstractQAtlas, so anything expressible here should stay small enough to
read at the declaration site.

Only the temperature axis is supported here.  A field derivative (`dF/dh`, for
`MagnetizationResponse`) is deliberately absent: `h` is a MODEL parameter, not a
fetch kwarg, so differentiating along it means reconstructing the model at each
step and every model spells its field differently.  That needs a model-parameter
mechanism, and inventing one silently for two relations would be worse than
leaving them uncovered.
"""
struct DerivedInput
    quantity::Type
    wrt::Symbol
    of::Function
    then::Function
    function DerivedInput(
        quantity::Type, wrt::Symbol; of::Function=(v, x) -> v, then::Function=identity
    )
        wrt in (:T, :β) || throw(
            ArgumentError(
                "DerivedInput: `wrt` must be :T or :β (the temperature axis); got " *
                ":$(wrt). Field derivatives need a model-parameter mechanism that " *
                "does not exist yet — see the docstring.",
            ),
        )
        return new(quantity, wrt, of, then)
    end
end

"""
    ∂(quantity, wrt)

Terse constructor for [`DerivedInput`](@ref): `∂(FreeEnergy, :T)` is `dF/dT`.
"""
function ∂(quantity::Type, wrt::Symbol; of::Function=(v, x) -> v, then::Function=identity)
    return DerivedInput(quantity, wrt; of=of, then=then)
end

"""
    ResponseEdge

A relation with at least one [`DerivedInput`](@ref) slot, materialized on every
hub implementing its quantity slots.  Shares the `name` / `sweep` / `finite_N` /
`exclusions` / `notes` / `references` vocabulary with the other edge kinds.

The tolerance is `max(default_rtol(backend), rtol_floor)`.  The backend part is
the accuracy of the DERIVATIVE; `rtol_floor` is the accuracy of the QUANTITIES
being related, which no amount of exact differentiation can improve on.

Measured motivation: IsingSquare's `SpecificHeat` and `ThermalEntropy` are
Onsager closed forms evaluated by quadrature, and `C` vs `T·dS/dT` agree to
2.5e-5 — with AD *and* with finite differences, to the same figure, so this is
the accuracy of the fetches, not of the derivative.  Holding that hub to AD's
1e-6 reported a quadrature residue as a physics failure.  An identity cannot be
verified more tightly than the values it relates.
"""
struct ResponseEdge
    name::Symbol
    relation::AbstractRelation
    subject::Symbol
    derived::NamedTuple
    rtol_floor::Float64
    sweep::NamedTuple
    finite_N::Int
    exclusions::Vector{Pair{Any,String}}
    notes::String
    references::Vector{String}
end

"""
    RESPONSES :: Vector{ResponseEdge}

Every declared response edge, in declaration order.
"""
const RESPONSES = ResponseEdge[]

# The quantity type behind a slot, or `nothing` if the slot is a field/untyped.
_quantity_slot(T) = (T isa Type && T <: AbstractQuantity) ? T : nothing
_field_slot(T) = (T isa Type && T <: AbstractQAtlas.AbstractField) ? T : nothing

"""
    response!(name; relation, derived, sweep, finite_N, exclusions, notes, references)

Declare a response edge.  The subject — the quantity slot the check solves for
and compares against `fetch` — is derived, not declared: a relation with exactly
one quantity slot has no ambiguity, and one with several is rejected rather than
guessed at.

Also rejected at declaration time: a duplicate name; a `derived` key that is not
an untyped slot of the relation; and any untyped slot left unsupplied (which
would make the check unrunnable in a way only visible at run time).
"""
function response!(
    name::Symbol;
    relation,
    derived::NamedTuple,
    rtol_floor::Real=0.0,
    sweep::NamedTuple=NamedTuple(),
    finite_N::Int=8,
    exclusions::AbstractVector=Pair{Any,String}[],
    notes::AbstractString="",
    references::AbstractVector{<:AbstractString}=String[],
)
    any(e -> e.name === name, RESPONSES) &&
        throw(ArgumentError("response!: :$(name) already declared"))
    rel = relation isa Type ? relation() : relation
    rel isa AbstractRelation ||
        throw(ArgumentError("response!: `relation` must be an AbstractQAtlas relation"))

    slots = variable_slots(rel)
    qslots = [(n, T) for (n, T) in slots if _quantity_slot(T) !== nothing]
    untyped = [n for (n, T) in slots if T === nothing]

    length(qslots) == 1 || throw(
        ArgumentError(
            "response!: :$(name) — $(nameof(typeof(rel))) has $(length(qslots)) quantity " *
            "slots; the subject to compare against `fetch` would be ambiguous. Only " *
            "single-subject relations are supported.",
        ),
    )
    for k in keys(derived)
        k in untyped || throw(
            ArgumentError(
                "response!: :$(name) — `$(k)` is not an untyped slot of " *
                "$(nameof(typeof(rel))); its untyped slots are $(untyped).",
            ),
        )
        derived[k] isa DerivedInput ||
            throw(ArgumentError("response!: :$(name) — `$(k)` must be a DerivedInput (∂)"))
    end
    missing_slots = setdiff(untyped, collect(keys(derived)))
    isempty(missing_slots) || throw(
        ArgumentError(
            "response!: :$(name) — untyped slot(s) $(missing_slots) are unsupplied, so " *
            "the check could never run. Supply them with ∂(...) or drop the edge.",
        ),
    )

    push!(
        RESPONSES,
        ResponseEdge(
            name,
            rel,
            first(qslots)[1],
            derived,
            Float64(rtol_floor),
            sweep,
            finite_N,
            Pair{Any,String}[p for p in exclusions],
            String(notes),
            String[r for r in references],
        ),
    )
    return nothing
end

"""
    @response(:name, relation = R, derived = (dF_dT = ∂(FreeEnergy, :T),), ...)

Keyword-macro front end for [`response!`](@ref), matching `@identity` / `@bound`.
"""
macro response(name, kwargs...)
    kws = [esc(k) for k in kwargs]
    return :(response!($(esc(name)); $(kws...)))
end

# ──────────────────────────────────────────────────────────────────────
# Generator — the :response kind of generated_checks()
# ──────────────────────────────────────────────────────────────────────

# The fetch kwargs and evaluation point for differentiating along `wrt` at a
# sweep point that carries `beta`.  β and T are one axis reparameterized, which
# is why only these two are allowed.
function _diff_axis(wrt::Symbol, point::NamedTuple)
    β = getfield(point, :beta)
    wrt === :β && return (x -> merge(point, (; beta=x)), float(β))
    return (x -> merge(point, (; beta=1 / x)), 1 / float(β))
end

function response_checks()
    out = GeneratedCheck[]
    backend = preferred_backend()
    for e in RESPONSES
        # The derivative's accuracy AND the fetched values' accuracy both bound
        # what this check can assert; the looser of the two wins.
        rtol = max(default_rtol(backend), e.rtol_floor)
        slots = variable_slots(e.relation)
        subject_T = only(T for (n, T) in slots if n === e.subject)
        for hub in _implemented_hubs(Type[subject_T])
            hub_id = string(
                "response/", e.name, "/", _kgshort(hub.model), "/", _kgshort(hub.bc)
            )
            reason = _exclusion_reason(e, hub.model, hub.bc)
            if reason !== nothing
                _push_excluded_check!(out, :response, hub_id, reason)
                continue
            end
            for point in _sweep_points(e.sweep)
                id = hub_id * _point_suffix(point)
                model_T, bc_T, rel = hub.model, hub.bc, e.relation
                runner = function ()
                    m = model_T()
                    bc = _bc_instance(bc_T; finite_N=e.finite_N)
                    args = Dict{Symbol,Any}()
                    for (n, T) in slots
                        n === e.subject && continue
                        F = _field_slot(T)
                        F === nothing && continue
                        # Field slots are the state point, not a fetch.
                        args[n] = if F === AbstractQAtlas.Temperature
                            1 / point.beta
                        else
                            point.beta
                        end
                    end
                    # Cross-check every derived input against an independent
                    # differentiation method.  A disagreement is evidence about
                    # the METHOD (an iterative solve AD differentiates through, a
                    # non-differentiable point) and must not be reported as a
                    # failed physical identity — so it becomes a visible skip.
                    for (n, d) in pairs(e.derived)
                        kw, x0 = _diff_axis(d.wrt, point)
                        g =
                            x -> d.of(
                                fetch(m, _quantity_instance(d.quantity), bc; kw(x)...),
                                x,
                            )
                        # A fetch can be non-differentiable by the chosen backend even
                        # when its signature says `::Real`: some implementations narrow to
                        # Float64 internally, and an AD dual then dies on a
                        # `MethodError(Float64, Dual)`.  (Measured: AKLT2D and Cluster1D
                        # do this for Energy{:per_site}.)  That is a fact about the
                        # METHOD, exactly like a cross-check disagreement, so it becomes a
                        # visible skip rather than an :error.
                        #
                        # Deliberately NOT a fallback to finite differences: this edge's
                        # rtol was computed for `backend` before the run, so silently
                        # substituting a lower-accuracy method would report an FD result
                        # at an AD tolerance.  Scoped to the derivative alone — a throw
                        # from the subject fetch or from `solve` is a real error and still
                        # surfaces as one.
                        local val, trusted, why
                        try
                            val, _, trusted, why = derivative_agreement(g, x0; primary=backend)
                        catch err
                            return _skip_outcome(
                                "$(n): not differentiable by $(nameof(typeof(backend))) " *
                                "— $(sprint(showerror, err))",
                            )
                        end
                        trusted || return _skip_outcome("$(n): $(why)")
                        args[n] = d.then(val)
                    end
                    expected = solve(rel, Val(e.subject); args...)
                    actual = fetch(m, _quantity_instance(subject_T), bc; point...)
                    return _outcome(
                        actual,
                        expected;
                        rtol=rtol,
                        atol=1e-10,
                        detail="derivative via $(nameof(typeof(backend)))",
                    )
                end
                push!(
                    out,
                    GeneratedCheck(
                        :response,
                        id,
                        string(
                            "response :",
                            e.name,
                            " (",
                            nameof(typeof(rel)),
                            ") on ",
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
    end
    return out
end

register_check_generator!(:response, response_checks)
register_edge_store!(:response, RESPONSES; location_of=e -> "response :$(e.name)")
