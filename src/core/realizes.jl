# core/realizes.jl — the model ↔ universality-class correspondence.
#
# A concrete `Model` *realizes* a `Universality{class}` at a critical point /
# scaling regime: its long-distance physics there is governed by that class.
# This is the relation that makes "which models are in the Ising class?"
# answerable, and the data behind the `by/universality` atlas view.
#
# Orthogonal to `REGISTRY` (which records how to *compute* a quantity): this
# records *membership* (model → class) and the regime where it holds.

"""
    Realization

One `model realizes class` row of the `REALIZES` correspondence: the
concrete `model` flows to `Universality{class}` in the stated `regime`
(e.g. a quantum critical point), resting on `references`.
"""
struct Realization
    model::Type
    class::Symbol
    regime::String
    at::Union{Function,Nothing}                  # predicate: is this instance at this critical locus?
    example::Union{AbstractQAtlasModel,Nothing}  # a representative critical instance on the locus
    references::Vector{String}
end

"""
    REALIZES :: Vector{Realization}

The model ↔ universality-class correspondence, populated at include-time by
[`realizes!`](@ref) / [`@realizes`](@ref).  Query with [`realizations`](@ref)
(by model) and [`realized_by`](@ref) (by class).
"""
const REALIZES = Realization[]

"""
    realizes!(model_T, class; regime, at=nothing, example=nothing, references=String[])

Record that `model_T` realizes `Universality{class}` in `regime`.  `class` is a
`Symbol` naming a universality class (`:Ising`, `:XY`, `:Heisenberg`, …).

`at` is an optional predicate `model_instance -> Bool` marking the *critical
locus* (a point, line, or surface in parameter space) where the model realizes
the class — multiple rows for one model must have mutually-exclusive `at`
predicates (a critical point belongs to exactly one class). `example` is a
representative critical instance on that locus (used to verify mutual exclusion
and to probe universal behaviour). Both are needed for the universal-quantity
delegation / verification to engage.
"""
function realizes!(
    model_T::Type,
    class::Symbol;
    regime::AbstractString,
    at::Union{Function,Nothing}=nothing,
    example::Union{AbstractQAtlasModel,Nothing}=nothing,
    references::AbstractVector{<:AbstractString}=String[],
)
    push!(
        REALIZES,
        Realization(
            model_T, class, String(regime), at, example, String[r for r in references]
        ),
    )
    return nothing
end

"""
    @realizes Model :class regime="…" references=[…]

Macro sugar around [`realizes!`](@ref): the positional `Model` is spliced as a
type and `:class` as the class symbol; the remaining `key=value` pairs are
forwarded as keyword arguments.
"""
macro realizes(model_T, class, kwargs...)
    return _forward_kw_macro(realizes!, :realizes, (model_T, class), kwargs)
end

"""
    realizations(model) -> Vector{NamedTuple}

The universality classes `model` realizes: `(class, regime, references)` rows.
"""
function realizations(model)
    m_T = _as_type(model)
    return [
        (class=r.class, regime=r.regime, references=r.references) for
        r in REALIZES if r.model === m_T
    ]
end

"""
    realized_by(class::Symbol) -> Vector{NamedTuple}

The concrete models realizing `Universality{class}`: `(model, regime,
references)` rows — the membership list behind the `by/universality` view.
"""
function realized_by(class::Symbol)
    return [
        (model=r.model, regime=r.regime, references=r.references) for
        r in REALIZES if r.class === class
    ]
end

"""
    realized_class(model) -> Union{Symbol,Nothing}

The universality class a `model` *instance* realizes at its current parameters:
the `class` of the unique [`@realizes`](@ref) row whose `at` predicate holds for
`model`, or `nothing` if the instance sits on no registered critical locus.
Errors if more than one row matches (non-exclusive `at` predicates — a coherence
violation; a critical point belongs to exactly one class).
"""
function realized_class(model::AbstractQAtlasModel)
    m_T = typeof(model)
    hits = Symbol[]
    for r in REALIZES
        r.model === m_T && r.at !== nothing && r.at(model) && push!(hits, r.class)
    end
    isempty(hits) && return nothing
    length(hits) == 1 || error(
        "realized_class($(m_T)): instance matches multiple classes $(hits) — the " *
        "@realizes `at` predicates for $(m_T) are not mutually exclusive.",
    )
    return only(hits)
end

function fetch(
    model::AbstractQAtlasModel,
    ::UniversalityClass,
    bc::BoundaryCondition;
    kwargs...,
)
    m_T = typeof(model)
    entries = [r for r in REALIZES if r.model === m_T]
    if isempty(entries)
        error("Model of type $(m_T) has no registered realizations.")
    end

    class = realized_class(model)
    if class !== nothing
        return Universality(class)
    end

    has_locus_predicates = any(r -> r.at !== nothing, entries)
    if has_locus_predicates
        throw(DomainError(model, "Model instance is not at any critical locus (off-critical/gapped)."))
    end

    if length(entries) == 1
        return Universality(entries[1].class)
    else
        throw(DomainError(model, "Model has multiple realization entries and no critical locus is active/specified."))
    end
end

function fetch(
    u::Universality{C},
    ::UniversalityClass,
    bc::BoundaryCondition;
    kwargs...,
) where {C}
    return u
end

