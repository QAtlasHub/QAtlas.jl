# core/about.jl — the one-line "what is this model" card.
#
# A `ModelCard` is the Wikipedia-first-sentence of a model: a one-sentence
# `summary` plus the `hamiltonian` in LaTeX, so a reader landing on the model
# page immediately sees *what the model is* and *what it solves*.  Orthogonal
# to `REGISTRY` (how to compute a quantity) and `REALIZES` (class membership):
# this records the human-facing description.
#
# Cards are optional and authored incrementally in `src/about_registry.jl`.
# Where no card exists, the docs generator falls back to the model's struct
# docstring — so every model gets *something*, and `@about` is the curated,
# math-rendered upgrade.

"""
    ModelCard

One `@about` row: the human-facing description of a `model` — a one-sentence
`summary`, the `hamiltonian` as a LaTeX string (rendered as display math; may
be empty), and optional `references` (bibkeys).  See [`@about`](@ref) /
[`about`](@ref).
"""
struct ModelCard
    model::Type
    summary::String
    hamiltonian::String
    references::Vector{String}
end

"""
    ABOUT :: Vector{ModelCard}

Module-level store of model description cards, populated at include-time by
[`about!`](@ref) / [`@about`](@ref) from `src/about_registry.jl`.  Query with
[`about`](@ref).
"""
const ABOUT = ModelCard[]

"""
    about!(model_T; summary, hamiltonian="", references=String[])

Record a [`ModelCard`](@ref) for `model_T`.  `summary` is a one-sentence
description (required, may contain inline `\$…\$` math); `hamiltonian` is a
LaTeX string rendered as display math on the model page (optional).  Usually
called via the [`@about`](@ref) macro.
"""
function about!(
    model_T::Type;
    summary::AbstractString,
    hamiltonian::AbstractString="",
    references::AbstractVector{<:AbstractString}=String[],
)
    isempty(strip(summary)) &&
        throw(ArgumentError("about!: summary must be a non-empty one-sentence description"))
    push!(
        ABOUT,
        ModelCard(
            model_T, String(summary), String(hamiltonian), String[r for r in references]
        ),
    )
    return nothing
end

"""
    @about Model summary="…" hamiltonian=raw"…" references=[…]

Macro sugar around [`about!`](@ref): the positional `Model` is spliced as a
type; the remaining `key=value` pairs are forwarded as keyword arguments.  Use
`raw"…"` for the `hamiltonian` so LaTeX backslashes survive.

```julia
@about TFIM summary="The 1D transverse-field Ising model, the canonical solvable quantum phase transition." \\
    hamiltonian=raw"H = -J\\sum_i \\sigma^z_i \\sigma^z_{i+1} - h\\sum_i \\sigma^x_i"
```
"""
macro about(model_T, kwargs...)
    return _forward_kw_macro(about!, :about, (model_T,), kwargs)
end

"""
    about(model) -> NamedTuple | nothing

The description card for `model` (instance or type): `(summary, hamiltonian,
references)`, or `nothing` if no [`@about`](@ref) card was authored.
"""
function about(model)
    m_T = _as_type(model)
    for c in ABOUT
        c.model === m_T &&
            return (summary=c.summary, hamiltonian=c.hamiltonian, references=c.references)
    end
    return nothing
end
