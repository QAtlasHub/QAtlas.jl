# ─────────────────────────────────────────────────────────────────────────────
# TightBinding1D — 1D non-interacting spinless-fermion chain.
#
# Hamiltonian (lattice spacing a = 1):
#
#     H = -t Σᵢ (c†ᵢ c_{i+1} + h.c.) - μ Σᵢ c†ᵢ cᵢ
#
# • t > 0 — nearest-neighbour hopping amplitude (real, strictly positive)
# • μ ∈ ℝ — chemical potential (arbitrary real)
#
# Single-band momentum-space dispersion (PBC ring, thermodynamic limit):
#
#     ε(k) = -2t cos(k) - μ,    k ∈ [-π, π).
#
# Ground state fills every Bloch level with ε(k) ≤ 0, i.e.
# cos(k) ≥ -μ/(2t).  This yields three regimes:
#
#   μ ≤ -2t            → empty band  (band bottom -2t already above 0)
#   |μ| < 2t           → partial filling, k_F = arccos(-μ/(2t)),
#                        gapless metal at the chemical potential
#   μ ≥ 2t             → full band (band top +2t already below 0)
#
# The mass gap (single-particle gap to the chemical potential) is then
#
#     Δ = max(0, |μ| - 2t),
#
# and at partial filling the Fermi velocity is
#
#     v_F = |∂ε/∂k|_{k_F} = 2t sin(k_F) = 2t √(1 - μ²/(4t²)).
#
# The ground-state energy density (per site, thermodynamic limit) is
#
#     E/N = (1/(2π)) ∫_{filled} ε(k) dk
#         = ⎧ 0                                       (μ ≤ -2t, empty)
#           ⎨ -(2t/π) sin(k_F) - (μ/π) k_F            (|μ| < 2t)
#           ⎩ -μ                                      (μ ≥ 2t, full)
#
# All three closed-form expressions are textbook (Ashcroft & Mermin,
# *Solid State Physics* (1976), chapter 9; equivalently any introductory
# many-body or condensed-matter text).  They serve as a clean reference
# point for free-fermion / Luttinger-liquid identities and as the
# non-interacting (U = 0) limit of `Hubbard1D` for spinless half-filling.
#
# References:
#   - N. W. Ashcroft, N. D. Mermin, *Solid State Physics*,
#     Harcourt College Publishers (1976), chapter 9.
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Fermion bilinears c†c
#   Observable:  Fermion (number n = c†c, bilinear ⟨c†_i c_j⟩); derived spin observables follow spin S = σ/2
#   Reference:   docs/src/conventions.md §Fermion convention

"""
    TightBinding1D(; t::Real = 1.0, μ::Real = 0.0) <: AbstractQAtlasModel

1D non-interacting spinless-fermion chain.

    H = -t Σᵢ (c†ᵢ c_{i+1} + h.c.) - μ Σᵢ c†ᵢ cᵢ

with `t > 0` (strict) and arbitrary real `μ`.

Quantities registered (Phase 1, all closed-form):

| Quantity                       | BC         | Formula                                  |
| ------------------------------ | ---------- | ---------------------------------------- |
| [`Energy`](@ref)`{:per_site}`  | `Infinite` | partial-filling integral / band edges    |
| [`MassGap`](@ref)              | `Infinite` | `max(0, |μ| - 2t)`                       |
| [`FermiVelocity`](@ref)        | `Infinite` | `2t √(1 - μ²/(4t²))` for `|μ| < 2t`, else 0 |

# Examples

```julia
julia> using QAtlas

julia> fetch(TightBinding1D(), Energy{:per_site}(), Infinite())   # half-filling
-0.6366197723675814

julia> fetch(TightBinding1D(; μ = 3.0), MassGap(), Infinite())     # band insulator
1.0

julia> fetch(TightBinding1D(; μ = 1.0), FermiVelocity(), Infinite())
1.7320508075688772
```

# References

- N. W. Ashcroft, N. D. Mermin, *Solid State Physics* (1976), chapter 9.
"""
struct TightBinding1D <: AbstractQAtlasModel
    t::Float64
    μ::Float64
    function TightBinding1D(t::Real, μ::Real)
        t > 0 || throw(DomainError(t, "TightBinding1D requires hopping t > 0; got t = $t."))
        return new(Float64(t), Float64(μ))
    end
end
TightBinding1D(; t::Real=1.0, μ::Real=0.0) = TightBinding1D(t, μ)

# Native energy granularity — the closed-form integral lives in the
# thermodynamic limit, so the per-site value is the primary one.
QAtlas.native_energy_granularity(::TightBinding1D, ::Infinite) = :per_site

# ═══════════════════════════════════════════════════════════════════════════════
# Energy{:per_site} at Infinite — closed-form ground-state energy density
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::TightBinding1D, ::Energy{:per_site}, ::Infinite;
          t=m.t, μ=m.μ, kwargs...) -> Float64

Ground-state energy per site in the thermodynamic limit:

    E/N = ⎧ 0                                       (μ ≤ -2t, empty band)
          ⎨ -(2t/π) sin(k_F) - (μ/π) k_F            (|μ| < 2t, k_F = arccos(-μ/(2t)))
          ⎩ -μ                                      (μ ≥ 2t,  full  band)

with `t > 0` enforced (DomainError otherwise).
"""
function fetch(
    m::TightBinding1D, ::Energy{:per_site}, ::Infinite; t::Real=m.t, μ::Real=m.μ, kwargs...
)
    t > 0 || throw(DomainError(t, "TightBinding1D Energy requires t > 0; got t = $t."))
    if μ ≤ -2t
        return 0.0                                 # empty band
    elseif μ ≥ 2t
        return -float(μ)                           # full band
    else
        k_F = acos(-μ / (2t))                      # partial filling
        return -(2t / π) * sin(k_F) - (μ / π) * k_F
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# MassGap at Infinite — closed form Δ = max(0, |μ| - 2t)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::TightBinding1D, ::MassGap, ::Infinite;
          t=m.t, μ=m.μ, kwargs...) -> Float64

Single-particle gap to the chemical potential:

    Δ = max(0, |μ| - 2t).

Gapless (metallic) for `|μ| ≤ 2t` (Fermi surface inside the band);
band insulator for `|μ| > 2t`.  At `|μ| = 2t` the band edge touches
the chemical potential exactly (Lifshitz transition), Δ = 0.
"""
function fetch(
    m::TightBinding1D, ::MassGap, ::Infinite; t::Real=m.t, μ::Real=m.μ, kwargs...
)
    t > 0 || throw(DomainError(t, "TightBinding1D MassGap requires t > 0; got t = $t."))
    return max(0.0, abs(float(μ)) - 2 * float(t))
end

# ═══════════════════════════════════════════════════════════════════════════════
# FermiVelocity at Infinite — closed form v_F = 2t sin(k_F)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::TightBinding1D, ::FermiVelocity, ::Infinite;
          t=m.t, μ=m.μ, kwargs...) -> Float64

Fermi velocity in the metallic regime (`|μ| < 2t`):

    v_F = |∂ε/∂k|_{k_F} = 2t sin(k_F) = 2t √(1 - μ²/(4t²)),

with `k_F = arccos(-μ/(2t))`.  In the gapped phase (`|μ| ≥ 2t`)
there is no Fermi surface and `v_F = 0` is returned by convention
(see [`MassGap`](@ref) for the insulating regime).
"""
function fetch(
    m::TightBinding1D, ::FermiVelocity, ::Infinite; t::Real=m.t, μ::Real=m.μ, kwargs...
)
    t > 0 ||
        throw(DomainError(t, "TightBinding1D FermiVelocity requires t > 0; got t = $t."))
    if abs(float(μ)) ≥ 2 * float(t)
        # Gapped phase: no Fermi surface, v_F = 0 by convention.
        return 0.0
    end
    return 2 * float(t) * sqrt(1 - (float(μ))^2 / (4 * float(t)^2))
end
