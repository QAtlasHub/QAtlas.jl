# ─────────────────────────────────────────────────────────────────────────────
# Universality{C} — parametric type for universality classes
#
# Each universality class is identified by a Symbol parameter C
# (e.g., :Ising, :XY, :Heisenberg, :Potts3, :Percolation, :KPZ).
# The spatial dimension d is passed as a keyword argument to `fetch`.
#
# Standard critical exponents {α, β, γ, δ, ν, η} are returned as a
# NamedTuple. For exact (analytically known) values, the fields are
# Rational{Int}. For numerical estimates, the fields are Float64 with
# corresponding `_err` fields indicating the uncertainty.
#
# Non-equilibrium classes (e.g., KPZ) have different exponent sets
# and use separate dispatch tags (GrowthExponents instead of
# CriticalExponents).
# ─────────────────────────────────────────────────────────────────────────────

"""
    Universality{C}

Parametric dispatch tag for universality classes. `C` is a `Symbol`
identifying the class (`:Ising`, `:XY`, `:Heisenberg`, `:Potts3`,
`:Potts4`, `:Percolation`, `:KPZ`, etc.).

Use with [`CriticalExponents`](@ref) (equilibrium) or
[`GrowthExponents`](@ref) (KPZ-type) and a `d` keyword to select
the spatial dimension:

```julia
QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2)   # exact Rational
QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=3)   # numerical + _err
QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=4)   # mean-field
```
"""
struct Universality{C} <: AbstractQAtlasModel end
Universality(name::Symbol) = Universality{name}()

"""
    CriticalExponents() <: AbstractQuantity

Standard set of equilibrium critical exponents
{α, β, γ, δ, ν, η} of a universality class. Returns a `NamedTuple`.

For exact values: fields are `Rational{Int}`.
For numerical estimates: fields are `Float64` with corresponding
`_err` fields (e.g., `β_err`) giving the uncertainty.
"""
struct CriticalExponents <: AbstractQuantity end

"""
    GrowthExponents() <: AbstractQuantity

KPZ-type growth / roughness / dynamic exponents.  Returns
`(β_growth, α_rough, z)` instead of the equilibrium set.
"""
struct GrowthExponents <: AbstractQuantity end
