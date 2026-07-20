
# `AbstractQAtlasModel` (the root model type) is defined once in AbstractQAtlas
# and imported + re-exported at the top of `QAtlas.jl` (#734).  The deprecation
# surface below (`AbstractModel` alias, `Model{M}` wrapper) stays in QAtlas and
# subtypes the shared root.

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

# `AbstractQuantity` (the root quantity type) is defined once in AbstractQAtlas
# and imported + re-exported at the top of `QAtlas.jl` (#734).  The legacy
# `Quantity{Q}` phantom wrapper below stays in QAtlas and subtypes the shared
# root.

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

# `fetch` is AbstractQAtlas's generic function — QAtlas `import`s it (see
# src/QAtlas.jl) and implements one method per `(model, quantity, bc)` triple in
# src/models/... (#734).  The generic fallback that errors on an un-implemented
# triple now lives in AbstractQAtlas; QAtlas no longer defines its own, so the
# whole ecosystem shares one `fetch` generic (the AbstractFFTs→FFTW idiom).

# Note: the legacy `fetch(::Symbol, ::Symbol, bc; kwargs...)` shim has
# been moved to `src/deprecate/legacy_fetch.jl` so the deprecation
# surface is concentrated in one place.  See src/deprecate/README.md.
