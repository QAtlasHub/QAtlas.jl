
"""
    AbstractQAtlasModel

Abstract parent type for every QAtlas model.  Concrete subtypes carry
their physics parameters as typed fields (e.g.
`struct TFIM <: AbstractQAtlasModel; J::Float64; h::Float64; end`).

The older `Model{S}(params::Dict)` phantom-typed wrapper is still an
`AbstractQAtlasModel` (via the deprecated alias below) but new models
must use concrete structs.
"""
abstract type AbstractQAtlasModel end

"""
    const AbstractModel = AbstractQAtlasModel

Backward-compatible alias.  Existing downstream code dispatches on
`::AbstractModel`; new code should use `::AbstractQAtlasModel` directly
or — preferably — a concrete model struct.
"""
const AbstractModel = AbstractQAtlasModel

"""
    Model{M} <: AbstractQAtlasModel  (deprecated)

Phantom-typed Dict wrapper kept for backward compatibility.  The
`Model(:TFIM; J=1.0, h=1.0)` constructor below still works but is
routed through the Symbol-dispatch deprecation shim in
`src/deprecate/legacy_fetch.jl`.  Prefer concrete model structs for
new code.
"""
struct Model{M} <: AbstractQAtlasModel
    params::Dict{Symbol,Any}
end
function Model(name::Symbol; kwargs...)
    canon = canonicalize_model(Val(name))
    return Model{canon}(Dict{Symbol,Any}(kwargs))
end

# `BoundaryCondition` / `Infinite` / `OBC` / `PBC` and the `_bc_size` helper are
# defined once in AbstractQAtlas (the single source of truth) and imported +
# re-exported at the top of `QAtlas.jl` (#734).  Concrete-model `fetch` methods
# dispatch on the same shared BC types, so behaviour is unchanged.

"""
    AbstractQuantity

Abstract parent type for quantities.  New code defines concrete structs
(e.g. `struct MagnetizationX <: AbstractQuantity end`) so dispatch is
static and naming is explicit (axis, entropy variant, …).  The older
`Quantity{S}` phantom-type wrapper is retained for legacy symbol-based
dispatch; see the `Quantity(::Symbol)` shim below.
"""
abstract type AbstractQuantity end

"""
    Quantity{Q} <: AbstractQuantity  (deprecated)

Phantom-typed wrapper kept for the legacy symbol API.  New code should
use concrete quantity structs such as `Energy()`, `MagnetizationX()`,
`ZZCorrelation(; mode=:static)`.
"""
struct Quantity{Q} <: AbstractQuantity end
function Quantity(q::Symbol)
    canon = canonicalize_quantity(Val(q))
    return Quantity{canon}()
end
Quantity(q::AbstractString) = Quantity(Symbol(q))

"""
    fetch(model, quantity, bc; kwargs...)

Return the stored / computed value of `quantity` for `model` under
boundary condition `bc`.  The canonical signature takes a concrete
model struct + concrete quantity struct + BC; a legacy
`fetch(::Symbol, ::Symbol, bc; kwargs...)` shim is also provided in
`src/deprecate/legacy_fetch.jl` for backward compatibility.

Each `(model, quantity, bc)` triple must be implemented as a separate
method; this top-level definition throws an informative error for
un-implemented triples.
"""
function fetch(
    model::AbstractQAtlasModel, quantity::AbstractQuantity, bc::BoundaryCondition; kwargs...
)
    return error(
        "QAtlas: no fetch method for model=$(typeof(model)), " *
        "quantity=$(typeof(quantity)), bc=$(typeof(bc)). " *
        "Define `fetch(::$(typeof(model)), ::$(typeof(quantity)), ::$(typeof(bc)); ...)` " *
        "in src/models/... to register the implementation.",
    )
end

# Note: the legacy `fetch(::Symbol, ::Symbol, bc; kwargs...)` shim has
# been moved to `src/deprecate/legacy_fetch.jl` so the deprecation
# surface is concentrated in one place.  See src/deprecate/README.md.
