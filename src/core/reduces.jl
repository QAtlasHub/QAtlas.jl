# core/reduces.jl — the model → model reduction correspondence.
#
# A concrete `source` model *reduces to* another concrete `target` model in a
# stated `regime` (a parameter limit or special point): there its Hamiltonian
# becomes the target's, so it may delegate quantities to the target.  This is
# the edge that makes a model→model `method=:delegation` legitimate (vs a
# model→class delegation, which is recorded by `REALIZES`).
#
# Orthogonal to `REALIZES` (model ↔ universality class) and to `REGISTRY` (how
# to compute a quantity): this records model ↔ model membership-under-limit.

"""
    Reduction

One `source reduces to target` row of the `REDUCES` correspondence:
the concrete `source` model becomes the concrete `target` model in the stated
`regime` (a limit / special point), resting on `references`.  This is what
makes a model→model delegation coherent — see [`@reduces`](@ref).
"""
struct Reduction
    source::Type
    target::Type
    regime::String
    references::Vector{String}
end

"""
    REDUCES :: Vector{Reduction}

The model ↔ model reduction correspondence, populated at include-time by
[`reduces!`](@ref) / [`@reduces`](@ref).  Query with [`reductions`](@ref) (by
source) and [`reduced_from`](@ref) (by target).
"""
const REDUCES = Reduction[]

"""
    reduces!(source_T, target_T; regime, references=String[])

Record that `source_T` reduces to `target_T` in `regime`.  Both arguments are
concrete model types.
"""
function reduces!(
    source_T::Type,
    target_T::Type;
    regime::AbstractString,
    references::AbstractVector{<:AbstractString}=String[],
)
    isempty(strip(regime)) && throw(
        ArgumentError(
            "reduces!: regime must be a non-empty description of the limit / special point",
        ),
    )
    push!(
        REDUCES,
        Reduction(source_T, target_T, String(regime), String[r for r in references]),
    )
    return nothing
end

"""
    @reduces Source Target regime="…" references=[…]

Macro sugar around [`reduces!`](@ref): the positional `Source` and `Target` are
spliced as model types; the remaining `key=value` pairs are forwarded as
keyword arguments.

```julia
@reduces MixedFieldIsing1D TFIM regime="longitudinal field h_z = 0"
```
"""
macro reduces(source_T, target_T, kwargs...)
    return _forward_kw_macro(reduces!, :reduces, (source_T, target_T), kwargs)
end

"""
    reductions(model) -> Vector{NamedTuple}

The models `model` reduces to: `(target, regime, references)` rows.
"""
function reductions(model)
    m_T = _as_type(model)
    return [
        (target=r.target, regime=r.regime, references=r.references) for
        r in REDUCES if r.source === m_T
    ]
end

"""
    reduced_from(model) -> Vector{NamedTuple}

The concrete models that reduce to `model`: `(source, regime, references)`
rows — the inverse of [`reductions`](@ref).
"""
function reduced_from(model)
    m_T = _as_type(model)
    return [
        (source=r.source, regime=r.regime, references=r.references) for
        r in REDUCES if r.target === m_T
    ]
end
