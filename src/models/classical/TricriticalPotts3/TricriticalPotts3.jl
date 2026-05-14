# ─────────────────────────────────────────────────────────────────────────────
# TricriticalPotts3 — tricritical 3-state Potts CFT M(6, 7).
#
# The tricritical 3-state Potts model is the multicritical fixed
# point of the 2-D q=3 Potts model with dilution / vacancies tuned
# to the tricritical line.  It is conformally realised by the
# unitary Virasoro minimal model M(p, p') with (p, p') = (6, 7),
# giving the central charge
#
#     c(6, 7) = 1 - 6 (6 − 7)² / (6 · 7) = 1 - 6/42 = 1 - 1/7 = 6/7.
#
# Andrews-Baxter-Forrester (1984) constructed it as an exactly
# solvable RSOS lattice model; Huse (1984) identified the tricritical
# 3-state Potts realisation.
#
# This Phase-1 entry registers `CentralCharge` only, delegating to
# the already-implemented [`MinimalModel(6, 7)`](@ref).  Exact
# critical exponents (ν = 5/9, β = 1/12, ...) and the ABF / Coulomb-gas
# operator content are tracked as Phase 2 and follow the same
# delegation pattern once the relevant quantity types
# (`CriticalExponents`, `PrimaryFields`) are wired through.
#
# References:
#   - G. E. Andrews, R. J. Baxter, P. J. Forrester,
#     J. Stat. Phys. 35, 193 (1984).
#   - D. A. Huse, Phys. Rev. B 30, 3908 (1984).
# ─────────────────────────────────────────────────────────────────────────────

"""
    TricriticalPotts3() <: AbstractQAtlasModel

Tricritical 3-state Potts CFT — the unitary Virasoro minimal model
M(6, 7), realised on the lattice by Andrews-Baxter-Forrester (1984)
RSOS models and by the q = 3 Potts model with dilution
(Huse 1984).

Quantities registered:

| Quantity                       | BC         | Method                              |
| ------------------------------ | ---------- | ----------------------------------- |
| [`CentralCharge`](@ref)        | `Infinite` | delegated to MinimalModel(7, 6)     |
| [`ConformalWeights`](@ref)     | `Infinite` | delegated to MinimalModel(7, 6)     |
| [`PrimaryFields`](@ref)        | `Infinite` | delegated to MinimalModel(7, 6)     |

# References

- G. E. Andrews, R. J. Baxter, P. J. Forrester,
  *J. Stat. Phys.* **35**, 193 (1984).
- D. A. Huse, *Phys. Rev. B* **30**, 3908 (1984).
"""
struct TricriticalPotts3 <: AbstractQAtlasModel end

# ═══════════════════════════════════════════════════════════════════════════════
# Central charge via MinimalModel(6, 7) delegation
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::TricriticalPotts3, ::CentralCharge, ::Infinite; kwargs...)
        -> Rational{Int}

Central charge of the tricritical 3-state Potts CFT, delegated to
[`MinimalModel(6, 7)`](@ref):

    c = 6/7.

Returned as an exact `Rational{Int}` (`6//7`).

# References

- G. E. Andrews, R. J. Baxter, P. J. Forrester,
  *J. Stat. Phys.* **35**, 193 (1984).
- D. A. Huse, *Phys. Rev. B* **30**, 3908 (1984).
"""
function fetch(::TricriticalPotts3, ::CentralCharge, ::Infinite; kwargs...)
    # MinimalModel(p, p_prime) requires p > p_prime; the tricritical
    # 3-state Potts CFT is M(7, 6) in this convention (= M(6, 7) in
    # literature where the order is reversed).  c is symmetric in the
    # two arguments so the result is the same: 6/7.
    return QAtlas.fetch(QAtlas.MinimalModel(7, 6), CentralCharge())
end

# ═══════════════════════════════════════════════════════════════════════════════
# Conformal weights via MinimalModel(7, 6) delegation (Phase 2)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::TricriticalPotts3, ::ConformalWeights, ::Infinite; r::Integer, s::Integer, kwargs...)
        -> Rational{Int}

Conformal weight `h_{r,s}` from the Kac table of the tricritical
3-state Potts CFT, delegated to [`MinimalModel(7, 6)`](@ref):

    h_{r,s}(p, p_prime) = ((p r - p_prime s)² - (p - p_prime)²) / (4 p p_prime),

with `(p, p_prime) = (7, 6)` in QAtlas' `p > p_prime` convention.

Index ranges follow [`MinimalModel`](@ref): `1 ≤ r ≤ p_prime - 1 = 5`
and `1 ≤ s ≤ p - 1 = 6`.  Returns `Rational{Int}`.

Examples:

- `h_{1,1} = 0`   (identity, lowest weight)
- `h_{1,2} = 1/7` (energy operator ε)
- `h_{2,1} = 3/8` (M(7,6) Kac (2,1))

# References

- A. A. Belavin, A. M. Polyakov, A. B. Zamolodchikov,
  *Nucl. Phys. B* **241**, 333 (1984).
- G. E. Andrews, R. J. Baxter, P. J. Forrester,
  *J. Stat. Phys.* **35**, 193 (1984).
"""
function fetch(
    ::TricriticalPotts3, ::ConformalWeights, ::Infinite; r::Integer, s::Integer, kwargs...
)
    return QAtlas.fetch(QAtlas.MinimalModel(7, 6), ConformalWeights(); r=r, s=s)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Primary fields via MinimalModel(7, 6) delegation (Phase 2)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::TricriticalPotts3, ::PrimaryFields, ::Infinite; kwargs...) -> Vector{NamedTuple}

Kac-table primary fields of the tricritical 3-state Potts CFT,
delegated to [`MinimalModel(7, 6)`](@ref).  Returns the
`((p_prime - 1) * (p - 1) / 2) = 15` independent primaries with
their `(r, s)` labels and `h` weights.

# References

- A. A. Belavin, A. M. Polyakov, A. B. Zamolodchikov,
  *Nucl. Phys. B* **241**, 333 (1984).
"""
function fetch(::TricriticalPotts3, ::PrimaryFields, ::Infinite; kwargs...)
    return QAtlas.fetch(QAtlas.MinimalModel(7, 6), PrimaryFields())
end
