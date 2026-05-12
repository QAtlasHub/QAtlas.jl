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

# ─────────────────────────────────────────────────────────────────────────────
# CFT Casimir / finite-size ground-state energy correction (Cardy 1986)
#
# For a 1+1D conformal field theory with central charge c and CFT
# velocity v on a system of size L the ground-state energy admits the
# expansion
#
#   E_0^PBC(L) = L · ε_∞               - π c v / (6  L) + O(L^{-2})
#   E_0^OBC(L) = L · ε_∞ + ε_surf      - π c v / (24 L) + O(L^{-2})
#
# Only the universal 1/L correction is exposed here.  Sign convention:
# the correction is *negative*; the system gains energy from the
# Casimir-like vacuum-mode shift on a finite cylinder / strip.
#
# Provenance: each per-class central charge below cites a primary
# reference for the 1+1D CFT, distinct from (but consistent with) the
# 2D classical universality data already shipped under
# `CriticalExponents`.  The PBC-to-OBC ratio of 4 is a *kinematic*
# consequence of the conformal map between cylinder and strip.
#
# Phase 2 (TODO, issue #150 follow-up): expose the primary tower
# (h, h̄) for each minimal model via a separate `ConformalTower`
# quantity.  Not implemented in this commit.
#
# References:
#   J. Cardy, Nucl. Phys. B 270, 186 (1986).
#   H. W. J. Blöte, J. L. Cardy, M. P. Nightingale,
#     Phys. Rev. Lett. 56, 742 (1986).
#   I. Affleck, Phys. Rev. Lett. 56, 746 (1986).
# ─────────────────────────────────────────────────────────────────────────────

raw"""
    _universality_central_charge(::Universality{C}) -> Real

Return the central charge `c` of the 1+1D CFT associated with
universality class `C`.  Used internally by
[`CasimirEnergyCorrection`](@ref) fetch.

Supported classes:

| `C`           | `c`     | Reference                                            |
|---------------|---------|------------------------------------------------------|
| `:Ising`      | `1//2`  | BPZ minimal model M(3,4); Cardy 1986                 |
| `:Potts3`     | `4//5`  | M(5,6) minimal model; Dotsenko–Fateev 1984           |
| `:Potts4`     | `1//1`  | Free-boson radius limit; Saleur 1987                 |
| `:XY`         | `1//1`  | Compact free boson (1+1D Luttinger liquid)           |
| `:Heisenberg` | `1//1`  | SU(2)_1 WZW (Affleck–Haldane); 1+1D AFM chain        |

Other classes raise `ErrorException`.  The 2D classical critical
percolation CFT is c = 0 (logarithmic, non-unitary) and is *not*
described by the simple Cardy formula in the same form, so it is
deliberately excluded here.  KPZ is non-equilibrium and has no CFT
central charge.  Mean-field is above the upper critical dimension and
has no 1+1D CFT representation.
"""
_universality_central_charge(::Universality{:Ising}) = 1 // 2
_universality_central_charge(::Universality{:Potts3}) = 4 // 5
_universality_central_charge(::Universality{:Potts4}) = 1 // 1
_universality_central_charge(::Universality{:XY}) = 1 // 1
_universality_central_charge(::Universality{:Heisenberg}) = 1 // 1

function _universality_central_charge(::Universality{C}) where {C}
    return error(
        "QAtlas Universality{:$C}: no 1+1D CFT central charge is " *
        "registered for this universality class.  CasimirEnergyCorrection " *
        "requires `c` from a unitary 1+1D CFT; classes such as :KPZ " *
        "(non-equilibrium), :Percolation (non-unitary, c = 0 logarithmic), " *
        "and :MeanField (no 1+1D CFT) are not supported.  If `C` should " *
        "have a 1+1D CFT entry, add a method to " *
        "`QAtlas._universality_central_charge` and document the source.",
    )
end

raw"""
    fetch(::Universality{C}, ::CasimirEnergyCorrection, ::PBC; L, v) -> Real

Return the universal Cardy 1/L correction
``-\pi c v / (6 L)`` at periodic boundary conditions, where `c` is
the central charge of the 1+1D CFT for class `C` (see
[`_universality_central_charge`](@ref)) and `v` is the CFT velocity
supplied by the caller.

`L` and `v` must be positive.  The return type is `Rational` when both
`c` is rational and `v` is rational/integer; otherwise `Float64`.
"""
function fetch(
    u::Universality{C}, ::CasimirEnergyCorrection, ::PBC; L::Real, v::Real, kwargs...
) where {C}
    L > 0 || throw(ArgumentError("CasimirEnergyCorrection: L must be positive; got $L"))
    v > 0 || throw(ArgumentError("CasimirEnergyCorrection: v must be positive; got $v"))
    c = _universality_central_charge(u)
    return -π * c * v / (6 * L)
end

raw"""
    fetch(::Universality{C}, ::CasimirEnergyCorrection, ::OBC; L, v) -> Real

Return the universal Cardy 1/L correction
``-\pi c v / (24 L)`` at open boundary conditions.

The PBC : OBC ratio of the 1/L term is exactly 4, a kinematic
consequence of the conformal map between the cylinder (PBC) and the
strip (OBC).
"""
function fetch(
    u::Universality{C}, ::CasimirEnergyCorrection, ::OBC; L::Real, v::Real, kwargs...
) where {C}
    L > 0 || throw(ArgumentError("CasimirEnergyCorrection: L must be positive; got $L"))
    v > 0 || throw(ArgumentError("CasimirEnergyCorrection: v must be positive; got $v"))
    c = _universality_central_charge(u)
    return -π * c * v / (24 * L)
end
