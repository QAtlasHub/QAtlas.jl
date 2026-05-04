# ─────────────────────────────────────────────────────────────────────────────
# KPZ (Kardar–Parisi–Zhang) universality class
#
# The KPZ equation is a non-equilibrium growth model. Its scaling
# exponents differ from equilibrium critical exponents and are accessed
# via the `GrowthExponents` tag rather than `CriticalExponents`.
#
# References (exact, d=1):
#   M. Kardar, G. Parisi, Y.-C. Zhang, Phys. Rev. Lett. 56, 889 (1986).
#   T. Sasamoto, H. Spohn, Nucl. Phys. B 834, 523 (2010).
#
# References (numerical, d≥2):
#   A. Pagnani, G. Parisi, Phys. Rev. E 92, 010101(R) (2015)  — d = 2.
#   J. Kelling, G. Ódor, Phys. Rev. E 84, 061150 (2011)       — d = 3.
# ─────────────────────────────────────────────────────────────────────────────

"""
    KPZ1D <: AbstractQAtlasModel

Dispatch tag for the **1+1-dimensional** KPZ (Kardar–Parisi–Zhang)
universality class. Acts as a convenience wrapper over
`Universality(:KPZ)` with the dimension pinned to `d = 1`:

```julia
QAtlas.fetch(KPZ1D(), CriticalExponents())
# ≡ QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d = 1)
```

The struct is kept as a distinct nominal type (rather than
`const KPZ1D = Universality{:KPZ}`) because the dimension is fixed
here — a type-level alias would defeat that fix. Subtyped to
`AbstractQAtlasModel` so dispatch composes with the canonical
`fetch(::AbstractQAtlasModel, ::AbstractQuantity, ::BoundaryCondition)`
signature.

See also: [`Universality`](@ref), [`GrowthExponents`](@ref).
"""
struct KPZ1D <: AbstractQAtlasModel end
function fetch(::KPZ1D, ::CriticalExponents; kwargs...)
    return fetch(Universality(:KPZ), GrowthExponents(); d=1, kwargs...)
end

"""
    fetch(::Universality{:KPZ}, ::GrowthExponents; d) -> NamedTuple

Scaling exponents of the KPZ universality class.

- **d = 1** (1+1D): all three exponents are exact `Rational{Int}`.
  Scaling relations: ``α + z = 2`` (Galilean invariance),
  ``β = α / z``.

- **d = 2** (2+1D): numerical estimates from Pagnani–Parisi 2015
  (large-scale RSOS simulations).
  Returns Float64 fields with `_err` companions:
  `β_growth = 0.2415(15)`, `α_rough = 0.393(5)`, `z = 1.613(9)`.

- **d = 3** (3+1D): numerical estimates from Kelling–Ódor 2011
  (octahedron model on GPU).
  Returns Float64 fields with `_err` companions:
  `β_growth ≈ 0.18`, `α_rough ≈ 0.31`, `z ≈ 1.51`. Galilean
  invariance ``α + z = 2`` is *not* strictly satisfied by these
  estimates — the discrepancy reflects estimation uncertainty (and
  ongoing debate about whether KPZ has an upper critical dimension)
  rather than a violation of the symmetry. Treat the d=3 entry as
  best-numerical, not as a sharp reference.

Higher dimensions (d ≥ 4) are not implemented; the upper critical
dimension of KPZ remains an open problem.
"""
function fetch(::Universality{:KPZ}, ::GrowthExponents; d::Int, kwargs...)
    if d == 1
        return (β_growth=1 // 3, α_rough=1 // 2, z=3 // 2)
    elseif d == 2
        # Pagnani & Parisi, PRE 92, 010101(R) (2015).  Errors are the
        # quoted statistical 1-σ from RSOS large-N extrapolation.
        return (
            β_growth=0.2415,
            β_growth_err=0.0015,
            α_rough=0.393,
            α_rough_err=0.005,
            z=1.613,
            z_err=0.009,
        )
    elseif d == 3
        # Kelling & Ódor, PRE 84, 061150 (2011).  The paper quotes
        # rough estimates rather than tight 1-σ; we adopt 0.01 as a
        # conservative uncertainty consistent with the cross-method
        # spread reported in the literature (cf. Halpin-Healy 2013).
        return (
            β_growth=0.18,
            β_growth_err=0.01,
            α_rough=0.31,
            α_rough_err=0.01,
            z=1.51,
            z_err=0.01,
        )
    end
    return error("KPZ growth exponents: d=$d not supported (only d ∈ {1, 2, 3}).")
end
