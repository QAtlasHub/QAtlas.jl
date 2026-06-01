# ─────────────────────────────────────────────────────────────────────────────
# Hubbard1D — Infinite() finite-T observables via regime-based delegation
#
# Phase-2A stopgap for issue #523 (Juttner-Klumper-Suzuki QTM NLIE).
# The full 4-NLIE QTM solver is large enough to warrant its own PR;
# this file gives Hubbard1D some Infinite() finite-T coverage right now
# by delegating to known limits and refusing intermediate regimes with
# NaN-and-warn.
#
# Regimes covered (half filling mu = U/2 only)
# ============================================
#
# (A)  U = 0 exact: free spinless fermion * 2 species
#        f_Hubbard(beta, t, 0) = 2 * f_TightBinding1D(beta, t, mu = 0)
#
# (B)  Very high T (beta * max(t, U) <= 0.05): atomic limit dominates,
#      every site has 4 equally weighted states
#        f -> -T ln 4 + O(beta * max(t, U))
#
# (C)  Strong coupling + low T (U/t >= 10 AND beta * J_eff >= 5,
#      J_eff = 4 t^2 / U): Anderson superexchange + Heisenberg c=1 CFT
#        f = e_0_LiebWu(t, U) - T^2 / (3 J_eff)
#      The Lieb-Wu integral gives the exact GS density (incl. higher
#      order in t/U); the CFT excess uses v_s = pi J_eff / 2.
#
# (D)  Otherwise: NaN + warn, pointing to #523 for the full QTM NLIE.
#
# References
# ==========
#
#   - E. H. Lieb, F. Y. Wu, Phys. Rev. Lett. 20, 1445 (1968) — GS BAE
#   - P. W. Anderson, Phys. Rev. 115, 2 (1959) — strong-coupling
#     superexchange J_eff = 4 t^2 / U
#   - G. Juttner, A. Klumper, J. Suzuki, Nucl. Phys. B 522, 471 (1998)
#     [arXiv:cond-mat/9711310] — full QTM NLIE; the U -> infinity
#     analytic limit (Sec. 7.1) confirms the Heisenberg AFM mapping
# ─────────────────────────────────────────────────────────────────────────────

const _HUBBARD_HIGH_T_CUTOFF = 0.05   # beta * max(t, U) below this -> high-T limit
const _HUBBARD_STRONG_U_RATIO = 10.0  # U/t above this -> strong-coupling regime
const _HUBBARD_CFT_BETA_MIN = 5.0     # beta * J_eff above this -> CFT regime

"""
    _hubbard1d_assert_half_filling(m::Hubbard1D)

Half-filling guard for Phase-2A stopgap (mu = U/2 within 1e-12). Off
half-filling raises `DomainError`. Doped finite-T is deferred to a
later phase together with the JKS NLIE.
"""
function _hubbard1d_assert_half_filling(m::Hubbard1D)
    if !isapprox(m.μ, m.U / 2; atol=1e-12)
        throw(
            DomainError(
                m.μ,
                "Hubbard1D finite-T (Phase-2A) is half-filling only; need μ = U/2 (= $(m.U/2)).",
            ),
        )
    end
    return nothing
end

"""
    _hubbard1d_thermal_regime(m::Hubbard1D, beta::Real) -> Symbol

Classify (m, β) into :u_zero, :high_T, :strong_low_T, or :intermediate
for delegation by the FreeEnergy/Entropy/SpecificHeat dispatches.
"""
function _hubbard1d_thermal_regime(m::Hubbard1D, beta::Real)
    if iszero(m.U)
        return :u_zero
    end
    if beta * max(m.t, m.U) ≤ _HUBBARD_HIGH_T_CUTOFF
        return :high_T
    end
    J_eff = 4 * m.t^2 / m.U
    if m.U / m.t ≥ _HUBBARD_STRONG_U_RATIO && beta * J_eff ≥ _HUBBARD_CFT_BETA_MIN
        return :strong_low_T
    end
    return :intermediate
end

"""
    _hubbard1d_intermediate_warn(quantity::Symbol, m::Hubbard1D, beta::Real)

NaN-return + warn for the intermediate regime (no delegate). Names #523
as the proper extension.
"""
function _hubbard1d_intermediate_warn(quantity::Symbol, m::Hubbard1D, beta::Real)
    J_eff_str = m.U > 0 ? "J_eff = $(round(4 * m.t^2 / m.U; digits=3))" : ""
    @warn (
        "Hubbard1D " *
        String(quantity) *
        " at Infinite(): (U/t = $(m.U/m.t), " *
        "β · t = $(beta * m.t)) lies in the intermediate regime that requires " *
        "the full Juttner-Klumper-Suzuki QTM NLIE (issue #523). Phase-2A stopgap " *
        "covers U = 0 exactly, very-high T (β · max(t,U) ≤ 0.05), and strong-" *
        "coupling low T (U/t ≥ 10 with β · J_eff ≥ 5" *
        (isempty(J_eff_str) ? "" : ", " * J_eff_str) *
        "). Returning NaN."
    )
    return NaN
end

# ── FreeEnergy ──────────────────────────────────────────────────────────────

"""
    fetch(m::Hubbard1D, ::FreeEnergy, ::Infinite; beta, kwargs...)

Per-site Helmholtz free energy of the half-filled infinite Hubbard1D
chain. Phase-2A stopgap: dispatches to (A) free-fermion × 2 at U = 0,
(B) -T ln 4 at very high T, (C) Lieb-Wu + Heisenberg CFT at strong
coupling + low T, (D) NaN + warn otherwise.
"""
function fetch(m::Hubbard1D, ::FreeEnergy, ::Infinite; beta::Real, kwargs...)
    isempty(kwargs) || @warn(
        "fetch(Hubbard1D, FreeEnergy, Infinite) received unrecognized kwargs; they are ignored.",
        kwargs=collect(keys(kwargs))
    )
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    m.t > 0 || throw(DomainError(m.t, "t must be > 0"))
    m.U ≥ 0 || throw(DomainError(m.U, "U must be ≥ 0"))
    _hubbard1d_assert_half_filling(m)

    regime = _hubbard1d_thermal_regime(m, beta)
    if regime === :u_zero
        tb = TightBinding1D(; t=m.t, μ=0.0)
        return 2 * QAtlas.fetch(tb, FreeEnergy(), Infinite(); beta=beta)
    elseif regime === :high_T
        return -log(4) / beta
    elseif regime === :strong_low_T
        e0 = QAtlas.fetch(m, GroundStateEnergyDensity(), Infinite())
        J_eff = 4 * m.t^2 / m.U
        T = 1 / beta
        return e0 - T^2 / (3 * J_eff)
    else
        return _hubbard1d_intermediate_warn(:FreeEnergy, m, beta)
    end
end

# ── ThermalEntropy ──────────────────────────────────────────────────────────

"""
    fetch(m::Hubbard1D, ::ThermalEntropy, ::Infinite; beta, kwargs...)

Per-site Gibbs entropy density. Same regime classification as
FreeEnergy: 2 × free-fermion at U = 0, ln 4 at very high T, 2T/(3 J_eff)
at strong coupling + low T, NaN + warn otherwise.
"""
function fetch(m::Hubbard1D, ::ThermalEntropy, ::Infinite; beta::Real, kwargs...)
    isempty(kwargs) || @warn(
        "fetch(Hubbard1D, ThermalEntropy, Infinite) received unrecognized kwargs; they are ignored.",
        kwargs=collect(keys(kwargs))
    )
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    m.t > 0 || throw(DomainError(m.t, "t must be > 0"))
    m.U ≥ 0 || throw(DomainError(m.U, "U must be ≥ 0"))
    _hubbard1d_assert_half_filling(m)

    regime = _hubbard1d_thermal_regime(m, beta)
    if regime === :u_zero
        tb = TightBinding1D(; t=m.t, μ=0.0)
        return 2 * QAtlas.fetch(tb, ThermalEntropy(), Infinite(); beta=beta)
    elseif regime === :high_T
        return log(4)
    elseif regime === :strong_low_T
        J_eff = 4 * m.t^2 / m.U
        T = 1 / beta
        return 2 * T / (3 * J_eff)
    else
        return _hubbard1d_intermediate_warn(:ThermalEntropy, m, beta)
    end
end

# ── SpecificHeat ────────────────────────────────────────────────────────────

"""
    fetch(m::Hubbard1D, ::SpecificHeat, ::Infinite; beta, kwargs...)

Per-site heat capacity density. 2 × free-fermion at U = 0; vanishes at
leading order in the very-high-T regime; 2T/(3 J_eff) at strong coupling
+ low T; NaN + warn otherwise.
"""
function fetch(m::Hubbard1D, ::SpecificHeat, ::Infinite; beta::Real, kwargs...)
    isempty(kwargs) || @warn(
        "fetch(Hubbard1D, SpecificHeat, Infinite) received unrecognized kwargs; they are ignored.",
        kwargs=collect(keys(kwargs))
    )
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    m.t > 0 || throw(DomainError(m.t, "t must be > 0"))
    m.U ≥ 0 || throw(DomainError(m.U, "U must be ≥ 0"))
    _hubbard1d_assert_half_filling(m)

    regime = _hubbard1d_thermal_regime(m, beta)
    if regime === :u_zero
        tb = TightBinding1D(; t=m.t, μ=0.0)
        return 2 * QAtlas.fetch(tb, SpecificHeat(), Infinite(); beta=beta)
    elseif regime === :high_T
        # Leading order: c_v -> 0 (entropy already saturated at ln 4).
        # Sub-leading O((β t)²) plus O((β U)²) corrections are not captured.
        return 0.0
    elseif regime === :strong_low_T
        J_eff = 4 * m.t^2 / m.U
        T = 1 / beta
        return 2 * T / (3 * J_eff)
    else
        return _hubbard1d_intermediate_warn(:SpecificHeat, m, beta)
    end
end
