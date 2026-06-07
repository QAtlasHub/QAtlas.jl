
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

"""
    BoundaryCondition

Abstract parent type.  The three concrete subtypes carry system-size
information where applicable, so `fetch` can read it from the
BC instead of `kwargs`.

- `Infinite` — thermodynamic limit; no size.
- `PBC(N::Int)` — periodic boundary conditions at finite `N`.
- `OBC(N::Int)` — open boundary conditions at finite `N`.

For backward compatibility, the zero-argument constructors `PBC()` and
`OBC()` exist and set `N = 0`, which signals "caller will pass `N` via
kwargs" — legacy fetch methods still look at `kwargs[:N]`.  New fetch
methods read `bc.N` directly.
"""
abstract type BoundaryCondition end

"""
    Infinite()

Thermodynamic-limit boundary condition — no finite size.
"""
struct Infinite <: BoundaryCondition end

"""
    OBC(N::Int)
    OBC(; N::Int = 0)

Open boundary condition.  `N` is the chain length.  `N = 0` is a legacy
sentinel meaning "size unspecified — caller passes it via kwargs";
`fetch` methods that accept `OBC(0)` must look up `kwargs[:N]`.
"""
struct OBC <: BoundaryCondition
    N::Int
end
OBC(; N::Int=0) = OBC(N)

"""
    PBC(N::Int)
    PBC(; N::Int = 0)

Periodic boundary condition.  See [`OBC`](@ref) for the `N = 0`
sentinel.
"""
struct PBC <: BoundaryCondition
    N::Int
end
PBC(; N::Int=0) = PBC(N)

"""
    _bc_size(bc::BoundaryCondition, kwargs) -> Int

Return the effective system size for `bc`.  Prefers `bc.N` when it is
positive; otherwise looks up `kwargs[:N]`; otherwise throws.  Legacy
fetch methods can use this helper to accept both `OBC(N=24)` and
`OBC(); N=24` call forms.
"""
function _bc_size end

function _bc_size(::Infinite, kwargs)
    return error("_bc_size called on Infinite; call a size-free fetch method instead")
end
function _bc_size(bc::Union{OBC,PBC}, kwargs)
    bc.N > 0 && return bc.N
    haskey(kwargs, :N) && return Int(kwargs[:N])
    return error("$(typeof(bc)): N unspecified — pass via OBC(N=...) or kwargs N=...")
end

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
