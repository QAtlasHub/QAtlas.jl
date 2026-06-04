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

One `model realizes class` row of the [`REALIZES`](@ref) correspondence: the
concrete `model` flows to `Universality{class}` in the stated `regime`
(e.g. a quantum critical point), resting on `references`.
"""
struct Realization
    model::Type
    class::Symbol
    regime::String
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
    realizes!(model_T, class; regime, references=String[])

Record that `model_T` realizes `Universality{class}` in `regime`.  `class` is a
`Symbol` naming a universality class (`:Ising`, `:XY`, `:Heisenberg`, …).
"""
function realizes!(
    model_T::Type,
    class::Symbol;
    regime::AbstractString,
    references::AbstractVector{<:AbstractString}=String[],
)
    push!(
        REALIZES, Realization(model_T, class, String(regime), String[r for r in references])
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
    kw = map(kwargs) do k
        k isa Expr && k.head === :(=) || error("@realizes: expected key=value, got $k")
        return Expr(:kw, k.args[1], esc(k.args[2]))
    end
    return Expr(:call, realizes!, Expr(:parameters, kw...), esc(model_T), esc(class))
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
