# ─────────────────────────────────────────────────────────────────────────────
# Haldane-Shastry — Infinite() finite-T observables via c = 1 CFT
#
# Stopgap path for issue #524 (semion gas finite-T). The Haldane-Shastry
# chain flows to the same c = 1 SU(2)_1 WZW CFT as the SU(2)-symmetric
# Heisenberg chain, with sound velocity v_s = π J / 2. The leading-order
# CFT thermodynamics are therefore identical to the Heisenberg case
# (issue #521 Path B, PR #526) modulo the ground-state energy density:
#
#     v_s = π J / 2
#     f(T) = e_0 - π T² / (6 v_s) = -π² J / 24 - T² / (3 J)
#     s(T) = π T / (3 v_s) = 2 T / (3 J)
#     c(T) = π T / (3 v_s) = 2 T / (3 J)
#
# Validity window (matches PR #526)
# =================================
#
# The Eggert-Affleck-Takahashi log corrections that apply to Heisenberg
# do *not* apply identically to HS — HS has no marginally irrelevant
# operator (it sits exactly on the SU(2)_1 fixed point), so the
# leading-order CFT is accurate to higher T than for Heisenberg.
# Nevertheless we keep the same conservative β > 5/J floor here
# until the full semion gas (Haldane 1991, BHR 2008) replaces this
# under #524, at which point the validity will be β > 0.
#
# References
# ==========
#
#   - F. D. M. Haldane, Phys. Rev. Lett. 60, 635 (1988) — model + e_0
#   - F. D. M. Haldane, Phys. Rev. Lett. 66, 1529 (1991) — Yangian,
#     spinon dispersion confirming v_s = π J / 2
#   - I. Affleck, Phys. Rev. Lett. 56, 746 (1986) — c=1 CFT form
#   - B. A. Bernevig, F. D. M. Haldane, N. Regnault, J. Phys. A 41,
#     304005 (2008) — semion gas thermodynamics (issue #524)
# ─────────────────────────────────────────────────────────────────────────────

const _HS_CFT_BETA_MIN = 5.0  # in units of 1/J; mirrors Heisenberg1D LO-CFT floor

"""
    _haldane_shastry_cft_freeenergy(J::Real, beta::Real) -> Float64

Leading-order c = 1 CFT free-energy density for HS:
`f = -π² J / 24 - T² / (3 J)`, with v_s = π J / 2.
"""
function _haldane_shastry_cft_freeenergy(J::Real, beta::Real)
    e0 = -π^2 * J / 24
    T = 1 / beta
    v_s = π * J / 2
    return e0 - π * T^2 / (6 * v_s)
end

"""
    _haldane_shastry_cft_entropy(J::Real, beta::Real) -> Float64

Leading-order c = 1 CFT entropy density: `s = 2 T / (3 J)`.
"""
function _haldane_shastry_cft_entropy(J::Real, beta::Real)
    T = 1 / beta
    v_s = π * J / 2
    return π * T / (3 * v_s)
end

"""
    _haldane_shastry_cft_specific_heat(J::Real, beta::Real) -> Float64

Leading-order c = 1 CFT specific heat: `c_v = 2 T / (3 J)`. Coincides
with the entropy at LO CFT.
"""
function _haldane_shastry_cft_specific_heat(J::Real, beta::Real)
    T = 1 / beta
    v_s = π * J / 2
    return π * T / (3 * v_s)
end

"""
    _haldane_shastry_cft_validity_warn(quantity::Symbol, J::Real, beta::Real)

NaN-return + warning above the LO CFT validity floor.
"""
function _haldane_shastry_cft_validity_warn(quantity::Symbol, J::Real, beta::Real)
    @warn (
        "HaldaneShastry " *
        String(quantity) *
        " at Infinite() uses a c=1 " *
        "CFT low-T expansion that is only validated for β > $(_HS_CFT_BETA_MIN)/J. " *
        "At β = $(beta) (T = $(round(1/beta; digits=3))) the LO term may carry > 5% " *
        "systematic error. The full Haldane semion gas (issue #524) will replace " *
        "this stopgap with the full β > 0 result. Returning NaN."
    )
    return NaN
end

# ── Dispatches ──────────────────────────────────────────────────────────────

"""
    fetch(::HaldaneShastry, ::FreeEnergy, ::Infinite; beta, kwargs...)

Per-site free energy of the infinite HS chain via leading-order c = 1
CFT: `f = -π² J / 24 - T² / (3 J)`. Valid β > 5/J; NaN+warn outside.
"""
function fetch(m::HaldaneShastry, ::FreeEnergy, ::Infinite; beta::Real, kwargs...)
    isempty(kwargs) || @warn(
        "fetch(HaldaneShastry, FreeEnergy, Infinite) received unrecognized kwargs; they are ignored.",
        kwargs=collect(keys(kwargs))
    )
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    if beta * m.J ≤ _HS_CFT_BETA_MIN
        return _haldane_shastry_cft_validity_warn(:FreeEnergy, m.J, beta)
    end
    return _haldane_shastry_cft_freeenergy(m.J, beta)
end

"""
    fetch(::HaldaneShastry, ::ThermalEntropy, ::Infinite; beta, kwargs...)

Per-site entropy via leading-order c = 1 CFT: `s = 2 T / (3 J)`.
"""
function fetch(m::HaldaneShastry, ::ThermalEntropy, ::Infinite; beta::Real, kwargs...)
    isempty(kwargs) || @warn(
        "fetch(HaldaneShastry, ThermalEntropy, Infinite) received unrecognized kwargs; they are ignored.",
        kwargs=collect(keys(kwargs))
    )
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    if beta * m.J ≤ _HS_CFT_BETA_MIN
        return _haldane_shastry_cft_validity_warn(:ThermalEntropy, m.J, beta)
    end
    return _haldane_shastry_cft_entropy(m.J, beta)
end

"""
    fetch(::HaldaneShastry, ::SpecificHeat, ::Infinite; beta, kwargs...)

Per-site specific heat via leading-order c = 1 CFT: `c_v = 2 T / (3 J)`.
Coincides with the entropy at LO.
"""
function fetch(m::HaldaneShastry, ::SpecificHeat, ::Infinite; beta::Real, kwargs...)
    isempty(kwargs) || @warn(
        "fetch(HaldaneShastry, SpecificHeat, Infinite) received unrecognized kwargs; they are ignored.",
        kwargs=collect(keys(kwargs))
    )
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    if beta * m.J ≤ _HS_CFT_BETA_MIN
        return _haldane_shastry_cft_validity_warn(:SpecificHeat, m.J, beta)
    end
    return _haldane_shastry_cft_specific_heat(m.J, beta)
end

# ─────────────────────────────────────────────────────────────────────────────
# Calabrese-Cardy entanglement at Infinite via Universality(:Heisenberg)
#
# The Haldane-Shastry chain is gapless with linear dispersion (free
# spinons of the SU(2)_1 WZW model) and shares the c = 1 free-boson
# universality class with the SU(2)-symmetric Heisenberg chain. Its
# single-interval von Neumann and Renyi entanglement entropies in the
# thermodynamic limit are therefore given by the Calabrese-Cardy
# closed forms with c = 1, identical to Heisenberg1D.
#
# Reference: Calabrese-Cardy J. Stat. Mech. P06002 (2004) §4.
# Issue: #580 Phase 1, extending the universality-layer delegation
# pattern introduced for Heisenberg1D.
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::HaldaneShastry, ::VonNeumannEntropy{:equilibrium}, ::Infinite;
          ℓ::Int, beta::Real = Inf, kwargs...) -> Float64

Single-interval von Neumann entanglement entropy of the
Haldane-Shastry chain in the thermodynamic limit, delegated to the
c = 1 Calabrese-Cardy form via `Universality(:Heisenberg)`.

- `beta = Inf` (default): T = 0 ground state, `S = (1/3) log ℓ`.
- `beta < Inf`: thermal state, `S = (1/3) log[(β/π) sinh(π ℓ / β)]`.
"""
function fetch(
    ::HaldaneShastry,
    ::VonNeumannEntropy{:equilibrium},
    ::Infinite;
    ℓ::Int,
    beta::Real=Inf,
    kwargs...,
)
    return fetch(
        Universality(:Heisenberg),
        VonNeumannEntropy(),
        Infinite();
        ℓ=ℓ,
        beta=beta,
        kwargs...,
    )
end

"""
    fetch(::HaldaneShastry, q::RenyiEntropy, ::Infinite;
          ℓ::Int, beta::Real = Inf, kwargs...) -> Float64

Single-interval Renyi-α entanglement entropy delegated to
`Universality(:Heisenberg)` with the standard
`c -> c · (1 + 1/α) / 2` substitution.
"""
function fetch(
    ::HaldaneShastry, q::RenyiEntropy, ::Infinite; ℓ::Int, beta::Real=Inf, kwargs...
)
    return fetch(Universality(:Heisenberg), q, Infinite(); ℓ=ℓ, beta=beta, kwargs...)
end
