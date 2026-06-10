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
# ═══════════════════════════════════════════════════════════════════════════════
# Finite-temperature thermodynamics — free-fermion BZ integrals
#
# Per-site grand-canonical quantities at fixed (β, μ).  Since the model
# Hamiltonian already absorbs μ into the single-particle dispersion
# ε(k) = -2t cos(k) - μ, the per-site "free energy" implemented below is
# the grand potential per site
#
#     ω(β; t, μ) = -(1/(2πβ)) ∫_{-π}^{π} log(1 + e^{-β ε(k)}) dk
#                = -(1/(πβ))  ∫_{0}^{π}  log(1 + e^{-β ε(k)}) dk
#
# (the k ↔ -k symmetry of cos folds the BZ in half).  The associated
# internal energy density u and entropy density s satisfy
#
#     u = (1/π) ∫_0^π ε(k) n_F(βε(k)) dk,
#     s = β (u - ω),
#     c_μ = (β²/π) ∫_0^π ε(k)² n_F(βε(k)) (1 - n_F(βε(k))) dk
#         = (∂u/∂T)_μ.
#
# All three integrals are textbook (e.g. Mahan §1.3, Coleman §2.4) and
# evaluated by `QuadGK.quadgk` with `rtol = 1e-10`.  Naming convention
# matches the existing `Energy{:per_site}` at T = 0: that quantity is
# the β → ∞ limit of `u` here, not of `ω`.  `Energy{:per_site}` at
# finite T is *not* registered here to avoid colliding with the T = 0
# closed form already on the file.
# ═══════════════════════════════════════════════════════════════════════════════

@inline _tb1d_dispersion(k::Real, t::Real, μ::Real) = -2 * t * cos(k) - μ

# Numerically stable log(1 + exp(y)).
@inline function _tb1d_log1pexp(y::Real)
    return y > 0 ? y + log1p(exp(-y)) : log1p(exp(y))
end

# Numerically stable Fermi-Dirac occupation n_F(βε) = 1/(1 + e^{βε}).
@inline function _tb1d_nF(βε::Real)
    return βε > 0 ? exp(-βε) / (1 + exp(-βε)) : 1 / (1 + exp(βε))
end

function _tb1d_thermo_infinite(quantity::Symbol, t::Real, μ::Real, β::Real)
    if quantity === :free_energy
        integrand = k -> _tb1d_log1pexp(-β * _tb1d_dispersion(k, t, μ))
        val, _ = quadgk(integrand, 0.0, π; rtol=1e-10)
        return -val / (π * β)
    elseif quantity === :entropy
        integrand = k -> begin
            εk = _tb1d_dispersion(k, t, μ)
            y = β * εk
            _tb1d_log1pexp(-y) + y * _tb1d_nF(y)
        end
        val, _ = quadgk(integrand, 0.0, π; rtol=1e-10)
        return val / π
    elseif quantity === :specific_heat
        integrand = k -> begin
            εk = _tb1d_dispersion(k, t, μ)
            n = _tb1d_nF(β * εk)
            εk^2 * n * (1 - n)
        end
        val, _ = quadgk(integrand, 0.0, π; rtol=1e-10)
        return β^2 * val / π
    else
        error("Unknown TightBinding1D thermal quantity: $quantity")
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# FreeEnergy at Infinite — grand-potential density
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::TightBinding1D, ::FreeEnergy, ::Infinite;
          beta::Real, t=m.t, μ=m.μ, kwargs...) -> Float64

Per-site grand-potential density of the 1D non-interacting tight-binding
chain at inverse temperature `β > 0`:

    ω(β; t, μ) = -(1/(πβ)) ∫_0^π log(1 + e^{-β ε(k)}) dk,
    ε(k) = -2t cos k - μ.

Evaluated by `QuadGK.quadgk` with `rtol = 1e-10`.  The T = 0
`Energy`{`:per_site`} closed form is the β → ∞ limit of the
*internal* energy `u(β)`, not of `ω(β)`; the two coincide only when
`μ` is outside the band (`|μ| > 2t`) and `u_kin → 0`.

# References
- G. D. Mahan, *Many-Particle Physics* (3rd ed., 2000), §1.3.
- P. Coleman, *Introduction to Many-Body Physics* (2015), §2.4.
"""
function fetch(
    m::TightBinding1D,
    ::FreeEnergy,
    ::Infinite;
    beta::Real,
    t::Real=m.t,
    μ::Real=m.μ,
    kwargs...,
)
    t > 0 || throw(DomainError(t, "TightBinding1D FreeEnergy requires t > 0; got t = $t."))
    beta > 0 ||
        throw(DomainError(beta, "TightBinding1D FreeEnergy requires β > 0; got β = $beta."))
    return _tb1d_thermo_infinite(:free_energy, t, μ, beta)
end

# ═══════════════════════════════════════════════════════════════════════════════
# ThermalEntropy at Infinite — per-site Gibbs entropy
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::TightBinding1D, ::ThermalEntropy, ::Infinite;
          beta::Real, t=m.t, μ=m.μ, kwargs...) -> Float64

Per-site Gibbs entropy

    s(β; t, μ) = (1/π) ∫_0^π { log(1 + e^{-βε}) + βε · n_F(βε) } dk,
    ε(k) = -2t cos k - μ,   n_F(x) = 1/(1 + e^x).

Equivalent to the thermodynamic identity `s = β (u − ω)`.  In the
high-T limit (`β → 0⁺`) each mode is half-occupied and `s → log 2`
per site; in the metallic ground-state limit (`β → ∞`, `|μ| < 2t`)
the Sommerfeld expansion gives `s ~ (π/3) v_F⁻¹ T` per site.

# References
- G. D. Mahan, *Many-Particle Physics* (3rd ed., 2000), §1.3.
"""
function fetch(
    m::TightBinding1D,
    ::ThermalEntropy,
    ::Infinite;
    beta::Real,
    t::Real=m.t,
    μ::Real=m.μ,
    kwargs...,
)
    t > 0 ||
        throw(DomainError(t, "TightBinding1D ThermalEntropy requires t > 0; got t = $t."))
    beta > 0 || throw(
        DomainError(beta, "TightBinding1D ThermalEntropy requires β > 0; got β = $beta."),
    )
    return _tb1d_thermo_infinite(:entropy, t, μ, beta)
end

# ═══════════════════════════════════════════════════════════════════════════════
# SpecificHeat at Infinite — per-site c_μ = (∂u/∂T)_μ
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::TightBinding1D, ::SpecificHeat, ::Infinite;
          beta::Real, t=m.t, μ=m.μ, kwargs...) -> Float64

Per-site specific heat at fixed chemical potential

    c_μ(β; t, μ) = (β²/π) ∫_0^π ε(k)² n_F(βε) (1 − n_F(βε)) dk,
    ε(k) = -2t cos k - μ.

Both T → 0 (metal: linear in T, Sommerfeld; insulator: Arrhenius)
and T → ∞ (~ β² · ⟨ε²⟩ / 4 → 0) limits are returned correctly.

# References
- G. D. Mahan, *Many-Particle Physics* (3rd ed., 2000), §1.3.
"""
function fetch(
    m::TightBinding1D,
    ::SpecificHeat,
    ::Infinite;
    beta::Real,
    t::Real=m.t,
    μ::Real=m.μ,
    kwargs...,
)
    t > 0 ||
        throw(DomainError(t, "TightBinding1D SpecificHeat requires t > 0; got t = $t."))
    beta > 0 || throw(
        DomainError(beta, "TightBinding1D SpecificHeat requires β > 0; got β = $beta.")
    )
    return _tb1d_thermo_infinite(:specific_heat, t, μ, beta)
end

function _tb1d_nmr_relaxation_infinite(t::Real, μ::Real, β::Real, η::Real)
    integrand_outer =
        k1 -> begin
            ε1 = _tb1d_dispersion(k1, t, μ)
            f1 = _tb1d_nF(β * ε1)
            integrand_inner = k2 -> begin
                ε2 = _tb1d_dispersion(k2, t, μ)
                f2 = _tb1d_nF(β * ε2)
                lorentz = η / (π * ((ε1 - ε2)^2 + η^2))
                return f1 * (1.0 - f2) * lorentz
            end
            val_inner, _ = quadgk(integrand_inner, 0.0, π; rtol=1e-6)
            return val_inner
        end
    val_outer, _ = quadgk(integrand_outer, 0.0, π; rtol=1e-6)
    return val_outer / π^2
end

# ═══════════════════════════════════════════════════════════════════════════════
# NMRSpinRelaxationRate at Infinite — per-site 1/T_1(β, η)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::TightBinding1D, ::NMRSpinRelaxationRate, ::Infinite;
          beta::Real, eta::Real=0.1, t=m.t, μ=m.μ, kwargs...) -> Float64

Per-site NMR spin-lattice relaxation rate \$1/T_1\$ of the 1D non-interacting tight-binding
chain at inverse temperature `β > 0`, in the Korringa golden-rule (contact-hyperfine,
\$ω → 0\$) form — the \$q\$-summed dynamic structure factor \$S(q, ω→0)\$ broadened by a
Lorentzian of width `η > 0`:

    1/T_1(β, η) = 1/π³ ∫_0^π dk₁ ∫_0^π dk₂ f(ε(k₁)) (1 - f(ε(k₂))) η / ((ε(k₁) - ε(k₂))² + η²)

where `ε(k) = -2t cos k - μ` and `f` is the Fermi function.

!!! note "η is a physical regulator, not a removable artifact"
    In a 1D band the energy-conserving (\$η → 0⁺\$) δ-function is cut off by the
    inverse-square-root van Hove singularities at the band edges, so the bare Korringa
    rate is **η-dependent**: this quantity is registered as `status=:approx` and is
    reported at finite `η` (default `0.1`). It is *not* the canonical 3D metallic
    Korringa law `1/T_1 ∝ T`; see Korringa (1950) for that limit. The value is verified
    against an independent discrete N-mode \$k\$-sum and the \$β → 0\$ ¼-factorisation
    `1/T_1 → ¼ · (η-broadened phase space)`.
"""
function fetch(
    m::TightBinding1D,
    ::NMRSpinRelaxationRate,
    ::Infinite;
    beta::Real,
    eta::Real=0.1,
    t::Real=m.t,
    μ::Real=m.μ,
    kwargs...,
)
    t > 0 || throw(
        DomainError(t, "TightBinding1D NMRSpinRelaxationRate requires t > 0; got t = \$t."),
    )
    beta > 0 || throw(
        DomainError(
            beta, "TightBinding1D NMRSpinRelaxationRate requires β > 0; got β = \$beta."
        ),
    )
    eta > 0 || throw(
        DomainError(
            eta, "TightBinding1D NMRSpinRelaxationRate requires η > 0; got η = \$eta."
        ),
    )
    return _tb1d_nmr_relaxation_infinite(t, μ, beta, eta)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Finite-size utilities (OBC/PBC)
# ═══════════════════════════════════════════════════════════════════════════════

function _tb1d_obc_spectrum(N::Int, t::Real, μ::Real)
    # ε_n = -2t cos(nπ / (N+1)) - μ, for n = 1,...,N
    return [-2 * t * cos(n * π / (N + 1)) - μ for n in 1:N]
end

function _tb1d_pbc_spectrum(N::Int, t::Real, μ::Real)
    # ε_n = -2t cos(2nπ / N) - μ, for n = 0,...,N-1
    return [-2 * t * cos(2 * n * π / N) - μ for n in 0:(N - 1)]
end

function _tb1d_energy_finite(eigenvalues::AbstractVector{<:Real})
    # T = 0 ground state: sum of all occupied level energies (ε ≤ 0)
    return sum(ε -> ε ≤ 0 ? ε : 0.0, eigenvalues)
end

function _tb1d_thermo_finite(
    quantity::Symbol, eigenvalues::AbstractVector{<:Real}, β::Real; kwargs...
)
    N = length(eigenvalues)
    if quantity === :free_energy
        return -sum(ε -> _tb1d_log1pexp(-β * ε), eigenvalues) / (N * β)
    elseif quantity === :entropy
        return sum(eigenvalues) do ε
            y = β * ε
            _tb1d_log1pexp(-y) + y * _tb1d_nF(y)
        end / N
    elseif quantity === :specific_heat
        return sum(eigenvalues) do ε
            y = β * ε
            n = _tb1d_nF(y)
            ε^2 * n * (1 - n)
        end * (β^2 / N)
    elseif quantity === :nmr_relaxation
        eta = Float64(get(kwargs, :eta, 0.1))
        eta > 0 || throw(
            DomainError(
                eta,
                "TightBinding1D NMRSpinRelaxationRate requires η > 0; got η = $eta.",
            ),
        )
        s = 0.0
        for ε1 in eigenvalues
            f1 = _tb1d_nF(β * ε1)
            for ε2 in eigenvalues
                f2 = _tb1d_nF(β * ε2)
                lorentz = eta / (π * ((ε1 - ε2)^2 + eta^2))
                s += f1 * (1.0 - f2) * lorentz
            end
        end
        return s / N^2
    else
        error("Unknown TightBinding1D finite quantity: $quantity")
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch dispatch for OBC and PBC
# ═══════════════════════════════════════════════════════════════════════════════

# ── OBC / PBC energy ──
function fetch(
    m::TightBinding1D,
    ::Energy{:total},
    bc::Union{OBC,PBC};
    t::Real=m.t,
    μ::Real=m.μ,
    kwargs...,
)
    t > 0 || throw(DomainError(t, "TightBinding1D Energy requires t > 0; got t = $t."))
    N = _bc_size(bc, kwargs)
    eigenvalues = bc isa OBC ? _tb1d_obc_spectrum(N, t, μ) : _tb1d_pbc_spectrum(N, t, μ)
    return _tb1d_energy_finite(eigenvalues)
end

function fetch(
    m::TightBinding1D,
    ::Energy{:per_site},
    bc::Union{OBC,PBC};
    t::Real=m.t,
    μ::Real=m.μ,
    kwargs...,
)
    t > 0 || throw(DomainError(t, "TightBinding1D Energy requires t > 0; got t = $t."))
    N = _bc_size(bc, kwargs)
    eigenvalues = bc isa OBC ? _tb1d_obc_spectrum(N, t, μ) : _tb1d_pbc_spectrum(N, t, μ)
    return _tb1d_energy_finite(eigenvalues) / N
end

# ── OBC / PBC one-particle gap ──
function fetch(
    m::TightBinding1D, ::MassGap, bc::Union{OBC,PBC}; t::Real=m.t, μ::Real=m.μ, kwargs...
)
    t > 0 || throw(DomainError(t, "TightBinding1D MassGap requires t > 0; got t = $t."))
    N = _bc_size(bc, kwargs)
    eigenvalues = bc isa OBC ? _tb1d_obc_spectrum(N, t, μ) : _tb1d_pbc_spectrum(N, t, μ)
    return minimum(abs.(eigenvalues))
end

# ── OBC / PBC thermal quantities ──
const _TB1D_THERMAL_QUANTITIES = (
    (FreeEnergy, :free_energy),
    (ThermalEntropy, :entropy),
    (SpecificHeat, :specific_heat),
    (NMRSpinRelaxationRate, :nmr_relaxation),
)

for (QTy, qsym) in _TB1D_THERMAL_QUANTITIES
    @eval begin
        function fetch(
            m::TightBinding1D,
            ::$QTy,
            bc::Union{OBC,PBC};
            beta::Real,
            t::Real=m.t,
            μ::Real=m.μ,
            kwargs...,
        )
            t > 0 || throw(
                DomainError(
                    t, "TightBinding1D $($(string(QTy))) requires t > 0; got t = $t."
                ),
            )
            beta > 0 || throw(
                DomainError(
                    beta,
                    "TightBinding1D $($(string(QTy))) requires β > 0; got β = $beta.",
                ),
            )
            N = _bc_size(bc, kwargs)
            eigenvalues =
                bc isa OBC ? _tb1d_obc_spectrum(N, t, μ) : _tb1d_pbc_spectrum(N, t, μ)
            return _tb1d_thermo_finite($(QuoteNode(qsym)), eigenvalues, beta; kwargs...)
        end
    end
end
