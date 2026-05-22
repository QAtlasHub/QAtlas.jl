# ─────────────────────────────────────────────────────────────────────────────
# TightBindingV1D — 1D spinless fermion chain with nearest-neighbour interaction
# (the "t-V model"; Jordan-Wigner-equivalent to the spin-1/2 XXZ chain).
#
# Hamiltonian:
#
#     H = -t  Σ_i (c†_i c_{i+1} + h.c.)
#         + V  Σ_i n_i n_{i+1}
#         - μ  Σ_i n_i ,
#
#     t > 0,   V ∈ ℝ,   μ ∈ ℝ.
#
# At V = 0 the chain reduces to the free-fermion tight-binding model with
# closed-form spectrum
#
#     ε(k) = -2t cos(k) - μ ,
#
# from which
#
#     MassGap         = max(0, |μ| - 2t)              (insulating for |μ| ≥ 2t)
#     FermiVelocity   = 2t · sin(arccos(-μ/(2t)))     (for |μ| < 2t)
#
# follow immediately.  At V ≠ 0 the model is, via the Jordan-Wigner
# transformation, equivalent to the XXZ chain with V/t setting the
# anisotropy Δ_XXZ; the gap and Fermi velocity then require the
# Yang-Yang (1966) Bethe-ansatz result.  Phase 1 of this file exposes
# only the V = 0 free-fermion point; V ≠ 0 raises `DomainError` and
# will be addressed in Phase 2 via the JW ↔ XXZ1D mapping.
#
# NOTE: We deliberately do NOT delegate to `TightBinding1D` (PR #316,
# not yet merged on `main`); the free-fermion closed forms are
# implemented inline here.
#
# References:
#   - C. N. Yang, C. P. Yang, "One-dimensional chain of anisotropic
#     spin-spin interactions", Phys. Rev. 150, 321 (1966).
#   - N. W. Ashcroft, N. D. Mermin, *Solid State Physics* (1976),
#     Chapter 9 (tight-binding bands).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Fermion bilinears c†c
#   Observable:  Fermion (number n = c†c, bilinear ⟨c†_i c_j⟩); derived spin observables follow spin S = σ/2
#   Reference:   docs/src/conventions.md §Fermion convention

"""
    TightBindingV1D(t::Real, V::Real, μ::Real)
    TightBindingV1D(; t=1.0, V=0.0, μ=0.0)

1D spinless-fermion t-V chain (NN hopping `t > 0`, NN density-density
interaction `V`, chemical potential `μ`).  At `V = 0` the model is the
free-fermion tight-binding chain; at `V ≠ 0` it is Jordan-Wigner
equivalent to the spin-1/2 XXZ chain (Phase 2).

The default keyword constructor lands at the closed-form free-fermion
point `(t, V, μ) = (1, 0, 0)`.
"""
struct TightBindingV1D <: AbstractQAtlasModel
    t::Float64
    V::Float64
    μ::Float64
    function TightBindingV1D(t::Real, V::Real, μ::Real)
        t > 0 || throw(DomainError(t, "TightBindingV1D requires t > 0; got t = $t."))
        return new(Float64(t), Float64(V), Float64(μ))
    end
end
TightBindingV1D(; t::Real=1.0, V::Real=0.0, μ::Real=0.0) = TightBindingV1D(t, V, μ)

# ═══════════════════════════════════════════════════════════════════════════════
# Energy granularity dispatch (see src/core/quantities.jl)
# ═══════════════════════════════════════════════════════════════════════════════

# Required so generic <-> per-site Energy conversions know that
# fetch(model::TightBindingV1D, ::Energy{:per_site}, ::Infinite) is the
# natively-implemented granularity; without this entry the verify()
# harness throws MethodError at the granularity-dispatch step (root
# cause of CI shard s07 fail across PRs #383-#387 and #447).
QAtlas.native_energy_granularity(::TightBindingV1D, ::Infinite) = :per_site

# ═══════════════════════════════════════════════════════════════════════════════
# Mass gap — V = 0 free-fermion closed form
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::TightBindingV1D, ::MassGap, ::Infinite;
          t=m.t, V=m.V, μ=m.μ, kwargs...) -> Float64

Single-particle gap of the tight-binding t-V chain at `V = 0`:

    Δ = max(0, |μ| - 2t) .

Gapless inside the band (|μ| < 2t), zero at the band edges (|μ| = 2t),
linear-in-|μ| insulator gap outside (|μ| > 2t).

`V ≠ 0` raises `DomainError` (Phase 2 via JW ↔ XXZ1D).
"""
function fetch(
    m::TightBindingV1D,
    ::MassGap,
    ::Infinite;
    t::Real=m.t,
    V::Real=m.V,
    μ::Real=m.μ,
    kwargs...,
)
    t > 0 || throw(DomainError(t, "TightBindingV1D MassGap requires t > 0; got t = $t."))
    if !iszero(V)
        throw(
            DomainError(
                V,
                "TightBindingV1D MassGap: V ≠ 0 (Jordan-Wigner-equivalent to interacting XXZ, " *
                "Yang-Yang 1966) deferred to Phase 2. Got V = $V.",
            ),
        )
    end
    return max(0.0, abs(μ) - 2t)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Fermi velocity — V = 0 free-fermion closed form
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::TightBindingV1D, ::FermiVelocity, ::Infinite;
          t=m.t, V=m.V, μ=m.μ, kwargs...) -> Float64

Fermi velocity of the V = 0 tight-binding t-V chain,

    v_F = 2t · sin(k_F),    k_F = arccos(-μ/(2t)) ,

valid in the metallic regime `|μ| < 2t`.  For `|μ| ≥ 2t` there is no
Fermi surface (band insulator); a `DomainError` is raised.

`V ≠ 0` raises `DomainError` (Phase 2 via JW ↔ XXZ1D).
"""
function fetch(
    m::TightBindingV1D,
    ::FermiVelocity,
    ::Infinite;
    t::Real=m.t,
    V::Real=m.V,
    μ::Real=m.μ,
    kwargs...,
)
    t > 0 ||
        throw(DomainError(t, "TightBindingV1D FermiVelocity requires t > 0; got t = $t."))
    if !iszero(V)
        throw(
            DomainError(
                V, "TightBindingV1D FermiVelocity: V ≠ 0 deferred to Phase 2. Got V = $V."
            ),
        )
    end
    abs(μ) < 2t || throw(
        DomainError(
            μ,
            "TightBindingV1D FermiVelocity: no Fermi surface for |μ| ≥ 2t (insulating regime). " *
            "Got μ = $μ, t = $t.",
        ),
    )
    return 2t * sin(acos(-μ / (2t)))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy per site — V = 0 free-fermion closed form
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::TightBindingV1D, ::Energy{:per_site}, ::Infinite;
          t=m.t, V=m.V, μ=m.μ, kwargs...) -> Float64

Ground-state energy density of the V = 0 tight-binding chain at T = 0:

    e₀ = -2t/π · sin(k_F) - μ/π · k_F ,   k_F = arccos(-μ/(2t))   (|μ| < 2t)
    e₀ = 0                                                          (μ ≤ -2t, empty band)
    e₀ = -μ                                                         (μ ≥  2t, filled band)

Obtained from `e₀ = (1/2π) ∫_{-k_F}^{k_F} ε(k) dk` with ε(k) = -2t cos(k) - μ.

`V ≠ 0` raises `DomainError` (Phase 2 via JW ↔ XXZ1D).

References:
  - G. D. Mahan, *Many-Particle Physics* (3rd ed., 2000), Chapter 1.
  - N. W. Ashcroft, N. D. Mermin, *Solid State Physics* (1976), Ch. 9.
"""
function fetch(
    m::TightBindingV1D,
    ::Energy{:per_site},
    ::Infinite;
    t::Real=m.t,
    V::Real=m.V,
    μ::Real=m.μ,
    kwargs...,
)
    t > 0 || throw(
        DomainError(t, "TightBindingV1D Energy{:per_site} requires t > 0; got t = $t.")
    )
    if !iszero(V)
        throw(
            DomainError(
                V,
                "TightBindingV1D Energy{:per_site}: V ≠ 0 (JW-equivalent to interacting XXZ, " *
                "Yang-Yang 1966) deferred to Phase 2. Got V = $V.",
            ),
        )
    end
    if μ <= -2t
        return 0.0
    elseif μ >= 2t
        return -float(μ)
    else
        k_F = acos(-μ / (2t))
        return -(2t / pi) * sin(k_F) - (μ / pi) * k_F
    end
end
# ═══════════════════════════════════════════════════════════════════════════════
# Finite-temperature thermodynamics — V = 0 free-fermion BZ integrals
#
# Identical integrands to the standalone `TightBinding1D` finite-T block
# (Mahan §1.3); we re-implement inline rather than delegating, matching
# the convention adopted for the T = 0 quantities in this file (see
# header note re PR #316).  V ≠ 0 raises `DomainError` and will be
# addressed in Phase 2 via the JW ↔ XXZ1D mapping.
# ═══════════════════════════════════════════════════════════════════════════════

@inline _tbv1d_dispersion(k::Real, t::Real, μ::Real) = -2 * t * cos(k) - μ

@inline function _tbv1d_log1pexp(y::Real)
    return y > 0 ? y + log1p(exp(-y)) : log1p(exp(y))
end

@inline function _tbv1d_nF(βε::Real)
    return βε > 0 ? exp(-βε) / (1 + exp(-βε)) : 1 / (1 + exp(βε))
end

function _tbv1d_thermo_infinite(quantity::Symbol, t::Real, μ::Real, β::Real)
    if quantity === :free_energy
        integrand = k -> _tbv1d_log1pexp(-β * _tbv1d_dispersion(k, t, μ))
        val, _ = quadgk(integrand, 0.0, π; rtol=1e-10)
        return -val / (π * β)
    elseif quantity === :entropy
        integrand = k -> begin
            εk = _tbv1d_dispersion(k, t, μ)
            y = β * εk
            _tbv1d_log1pexp(-y) + y * _tbv1d_nF(y)
        end
        val, _ = quadgk(integrand, 0.0, π; rtol=1e-10)
        return val / π
    elseif quantity === :specific_heat
        integrand = k -> begin
            εk = _tbv1d_dispersion(k, t, μ)
            n = _tbv1d_nF(β * εk)
            εk^2 * n * (1 - n)
        end
        val, _ = quadgk(integrand, 0.0, π; rtol=1e-10)
        return β^2 * val / π
    else
        error("Unknown TightBindingV1D thermal quantity: $quantity")
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# FreeEnergy at Infinite — V = 0 free-fermion grand-potential density
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::TightBindingV1D, ::FreeEnergy, ::Infinite;
          beta::Real, t=m.t, V=m.V, μ=m.μ, kwargs...) -> Float64

Per-site grand-potential density of the spinless-fermion t-V chain at
inverse temperature `β > 0`, currently implemented only at `V = 0`:

    ω(β; t, μ) = -(1/(πβ)) ∫_0^π log(1 + e^{-β ε(k)}) dk,
    ε(k) = -2t cos k - μ.

`V ≠ 0` raises `DomainError` (Phase 2 via JW ↔ XXZ1D).
"""
function fetch(
    m::TightBindingV1D,
    ::FreeEnergy,
    ::Infinite;
    beta::Real,
    t::Real=m.t,
    V::Real=m.V,
    μ::Real=m.μ,
    kwargs...,
)
    t > 0 || throw(DomainError(t, "TightBindingV1D FreeEnergy requires t > 0; got t = $t."))
    beta > 0 || throw(
        DomainError(beta, "TightBindingV1D FreeEnergy requires β > 0; got β = $beta.")
    )
    if !iszero(V)
        throw(
            DomainError(
                V,
                "TightBindingV1D FreeEnergy: V ≠ 0 (Jordan-Wigner-equivalent to interacting XXZ, " *
                "Yang-Yang 1966) deferred to Phase 2. Got V = $V.",
            ),
        )
    end
    return _tbv1d_thermo_infinite(:free_energy, t, μ, beta)
end

# ═══════════════════════════════════════════════════════════════════════════════
# ThermalEntropy at Infinite — V = 0 Gibbs entropy density
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::TightBindingV1D, ::ThermalEntropy, ::Infinite;
          beta::Real, t=m.t, V=m.V, μ=m.μ, kwargs...) -> Float64

Per-site Gibbs entropy of the V = 0 tight-binding t-V chain:

    s(β; t, μ) = (1/π) ∫_0^π { log(1 + e^{-βε}) + βε · n_F(βε) } dk,
    ε(k) = -2t cos k - μ.

`V ≠ 0` raises `DomainError` (Phase 2 via JW ↔ XXZ1D).
"""
function fetch(
    m::TightBindingV1D,
    ::ThermalEntropy,
    ::Infinite;
    beta::Real,
    t::Real=m.t,
    V::Real=m.V,
    μ::Real=m.μ,
    kwargs...,
)
    t > 0 ||
        throw(DomainError(t, "TightBindingV1D ThermalEntropy requires t > 0; got t = $t."))
    beta > 0 || throw(
        DomainError(beta, "TightBindingV1D ThermalEntropy requires β > 0; got β = $beta."),
    )
    if !iszero(V)
        throw(
            DomainError(
                V, "TightBindingV1D ThermalEntropy: V ≠ 0 deferred to Phase 2. Got V = $V."
            ),
        )
    end
    return _tbv1d_thermo_infinite(:entropy, t, μ, beta)
end

# ═══════════════════════════════════════════════════════════════════════════════
# SpecificHeat at Infinite — V = 0 free-fermion c_μ
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::TightBindingV1D, ::SpecificHeat, ::Infinite;
          beta::Real, t=m.t, V=m.V, μ=m.μ, kwargs...) -> Float64

Per-site specific heat at fixed chemical potential of the V = 0
tight-binding t-V chain:

    c_μ(β; t, μ) = (β²/π) ∫_0^π ε(k)² n_F(βε) (1 − n_F(βε)) dk,
    ε(k) = -2t cos k - μ.

`V ≠ 0` raises `DomainError` (Phase 2 via JW ↔ XXZ1D).
"""
function fetch(
    m::TightBindingV1D,
    ::SpecificHeat,
    ::Infinite;
    beta::Real,
    t::Real=m.t,
    V::Real=m.V,
    μ::Real=m.μ,
    kwargs...,
)
    t > 0 ||
        throw(DomainError(t, "TightBindingV1D SpecificHeat requires t > 0; got t = $t."))
    beta > 0 || throw(
        DomainError(beta, "TightBindingV1D SpecificHeat requires β > 0; got β = $beta.")
    )
    if !iszero(V)
        throw(
            DomainError(
                V, "TightBindingV1D SpecificHeat: V ≠ 0 deferred to Phase 2. Got V = $V."
            ),
        )
    end
    return _tbv1d_thermo_infinite(:specific_heat, t, μ, beta)
end
