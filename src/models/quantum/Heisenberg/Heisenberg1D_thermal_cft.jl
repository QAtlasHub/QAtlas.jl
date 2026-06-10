# ─────────────────────────────────────────────────────────────────────────────
# Heisenberg1D — Infinite() finite-T observables via c = 1 CFT low-T expansion
#
# Path B of issue #521 (stopgap before the full Klümper NLIE Δ → 1 limit).
# This file fills the FreeEnergy / ThermalEntropy / SpecificHeat gap at
# the SU(2)-symmetric point so external benchmarks have an analytic
# Infinite() reference inside the CFT-valid window.
#
# Physics
# =======
#
# The S = 1/2 antiferromagnetic Heisenberg chain
#
#     H = J · Σᵢ Sᵢ · Sᵢ₊₁
#
# is critical with central charge c = 1 and sound velocity
#
#     v_s = π J / 2          (Δ → 1 limit of v_s = π J sin γ / (2 γ))
#
# Leading-order CFT then gives per-site
#
#     f(T) ≈ e₀ − π T² / (6 v_s) = e₀ − T² / (3 J)
#     s(T) ≈ π T / (3 v_s)        =  2 T / (3 J)
#     c(T) ≈ π T / (3 v_s)        =  2 T / (3 J)
#
# with e₀ = J (1/4 − ln 2) the Hulthén ground-state energy density.
#
# Validity
# ========
#
# Eggert-Affleck-Takahashi (1994) log corrections modify each leading
# term by a factor ≈ 1 + 1 / (2 ln(T₀/T)) with T₀ ≈ 7.7 J. The plain
# LO expression here is accurate to ≲ 5 % for T ≲ J/5 (β ≥ 5/J) and
# rapidly degrades above. We refuse to return a value outside the
# validity window:
#
#     β > 5/J  → CFT value
#     β ≤ 5/J  → NaN + warning naming the Klümper-Δ→1 NLIE path
#                (issue #521) as the proper extension
#
# References
# ==========
#
#   - I. Affleck, Phys. Rev. Lett. 56, 746 (1986)
#   - H. W. J. Blöte, J. L. Cardy, M. P. Nightingale, Phys. Rev. Lett.
#     56, 742 (1986)
#   - S. Lukyanov, A. Zamolodchikov, Nucl. Phys. B 493, 571 (1997)
#     — sound velocity at the isotropic point
#   - S. Eggert, I. Affleck, M. Takahashi, Phys. Rev. Lett. 73, 332
#     (1994) — multiplicative log corrections
#   - A. Klümper, Z. Phys. B 91, 507 (1993) — Δ → 1 limit of XXZ NLIE
# ─────────────────────────────────────────────────────────────────────────────

const _HEIS_CFT_BETA_MIN = 5.0  # in units of 1/J; below this β the LO CFT degrades > 5 %

"""
    _heisenberg1d_cft_freeenergy(J::Real, beta::Real) -> Float64

Leading-order c = 1 CFT free-energy density for the SU(2) Heisenberg chain:
`f = e₀ − T² / (3J)`, with `e₀ = J(1/4 − ln 2)` and `v_s = π J / 2`.
"""
function _heisenberg1d_cft_freeenergy(J::Real, beta::Real)
    e0 = J * (0.25 - log(2))
    T = 1 / beta
    v_s = π * J / 2
    return e0 - π * T^2 / (6 * v_s)
end

"""
    _heisenberg1d_cft_entropy(J::Real, beta::Real) -> Float64

Leading-order c = 1 CFT entropy density: `s = 2T / (3J) = π T / (3 v_s)`.
"""
function _heisenberg1d_cft_entropy(J::Real, beta::Real)
    T = 1 / beta
    v_s = π * J / 2
    return π * T / (3 * v_s)
end

"""
    _heisenberg1d_cft_specific_heat(J::Real, beta::Real) -> Float64

Leading-order c = 1 CFT specific-heat density: `c_v = 2T / (3J)`.
Equals `s(T)` at LO CFT (`c_v = T ∂s/∂T = s` for linear-in-T entropy).
"""
function _heisenberg1d_cft_specific_heat(J::Real, beta::Real)
    T = 1 / beta
    v_s = π * J / 2
    return π * T / (3 * v_s)
end

"""
    _heisenberg1d_cft_validity_warn(quantity::Symbol, J::Real, beta::Real)

Emit a `@warn` naming Klümper-Δ→1 (#521 Path A) as the proper extension
and return NaN. Called when β is below the LO-CFT validity floor.
"""
function _heisenberg1d_cft_validity_warn(quantity::Symbol, J::Real, beta::Real)
    @warn (
        "Heisenberg1D " *
        String(quantity) *
        " at Infinite() uses a c=1 " *
        "CFT low-T expansion that is only accurate for β > $(_HEIS_CFT_BETA_MIN)/J. " *
        "At β = $(beta) (T = $(round(1/beta; digits=3))) the LO term has > 5% " *
        "systematic error from EAT log corrections. The full Klümper NLIE " *
        "Δ → 1 limit (issue #521 Path A) will replace this. Returning NaN."
    )
    return NaN
end

# ── Dispatches ──────────────────────────────────────────────────────────────

"""
    fetch(::Heisenberg1D, ::FreeEnergy, ::Infinite; beta, J=1.0)

Per-site free energy of the infinite spin-1/2 Heisenberg AF chain via
leading-order c = 1 CFT. Returns `e₀ - π T² / (6 v_s)` for β > 5/J;
otherwise NaN + warn. `β ≤ 0` raises `DomainError`.
"""
function fetch(::Heisenberg1D, ::FreeEnergy, ::Infinite; beta::Real, J::Real=1.0, kwargs...)
    isempty(kwargs) || @warn(
        "fetch(Heisenberg1D, FreeEnergy, Infinite) received unrecognized kwargs; they are ignored.",
        kwargs=collect(keys(kwargs))
    )
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    J > 0 || throw(DomainError(J, "J must be > 0"))
    if beta * J ≤ _HEIS_CFT_BETA_MIN
        return _heisenberg1d_cft_validity_warn(:FreeEnergy, J, beta)
    end
    return _heisenberg1d_cft_freeenergy(J, beta)
end

"""
    fetch(::Heisenberg1D, ::ThermalEntropy, ::Infinite; beta, J=1.0)

Per-site entropy via leading-order c = 1 CFT: `s = π T / (3 v_s) = 2T / (3J)`.
"""
function fetch(
    ::Heisenberg1D, ::ThermalEntropy, ::Infinite; beta::Real, J::Real=1.0, kwargs...
)
    isempty(kwargs) || @warn(
        "fetch(Heisenberg1D, ThermalEntropy, Infinite) received unrecognized kwargs; they are ignored.",
        kwargs=collect(keys(kwargs))
    )
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    J > 0 || throw(DomainError(J, "J must be > 0"))
    if beta * J ≤ _HEIS_CFT_BETA_MIN
        return _heisenberg1d_cft_validity_warn(:ThermalEntropy, J, beta)
    end
    return _heisenberg1d_cft_entropy(J, beta)
end

"""
    fetch(::Heisenberg1D, ::SpecificHeat, ::Infinite; beta, J=1.0)

Per-site heat capacity via leading-order c = 1 CFT: `c_v = π T / (3 v_s) = 2T / (3J)`.
"""
function fetch(
    ::Heisenberg1D, ::SpecificHeat, ::Infinite; beta::Real, J::Real=1.0, kwargs...
)
    isempty(kwargs) || @warn(
        "fetch(Heisenberg1D, SpecificHeat, Infinite) received unrecognized kwargs; they are ignored.",
        kwargs=collect(keys(kwargs))
    )
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    J > 0 || throw(DomainError(J, "J must be > 0"))
    if beta * J ≤ _HEIS_CFT_BETA_MIN
        return _heisenberg1d_cft_validity_warn(:SpecificHeat, J, beta)
    end
    return _heisenberg1d_cft_specific_heat(J, beta)
end

# -----------------------------------------------------------------------------
# Adjacent-interval quantum-information triple at the gapless point
# (Calabrese-Cardy 2009 MI / CC-Tonni 2012 LN / CC 2005 saturation)
# -----------------------------------------------------------------------------

"""
    fetch(::Heisenberg1D, ::MutualInformation, ::Infinite;
          ℓ_A::Real, ℓ_B::Real, beta::Real=Inf, kwargs...) -> Float64

Calabrese-Cardy mutual information of two adjacent intervals,
delegated to `Universality(:Heisenberg)` (c = 1). The chain is gapless
SU(2)-invariant Luttinger liquid, so the Universality formula applies
directly.
"""
function fetch(
    ::Heisenberg1D,
    ::MutualInformation,
    ::Infinite;
    ℓ_A::Real,
    ℓ_B::Real,
    beta::Real=Inf,
    kwargs...,
)
    return fetch(
        Universality(:Heisenberg),
        MutualInformation(),
        Infinite();
        ℓ_A=ℓ_A,
        ℓ_B=ℓ_B,
        beta=beta,
        kwargs...,
    )
end

"""
    fetch(::Heisenberg1D, ::LogarithmicNegativity, ::Infinite;
          ℓ_A::Real, ℓ_B::Real, kwargs...) -> Float64

Calabrese-Cardy-Tonni 2012 logarithmic negativity of two adjacent
intervals, delegated to `Universality(:Heisenberg)`:

    E = (1/4) log[ℓ_A * ℓ_B / (ℓ_A + ℓ_B)].
"""
function fetch(
    ::Heisenberg1D, ::LogarithmicNegativity, ::Infinite; ℓ_A::Real, ℓ_B::Real, kwargs...
)
    return fetch(
        Universality(:Heisenberg),
        LogarithmicNegativity(),
        Infinite();
        ℓ_A=ℓ_A,
        ℓ_B=ℓ_B,
        kwargs...,
    )
end

"""
    fetch(::Heisenberg1D, ::EntanglementSaturationDensity, ::Infinite;
          beta_eff::Real, kwargs...) -> Float64

Long-time post-quench saturation of half-system entanglement entropy
per unit length, delegated to `Universality(:Heisenberg)`:

    S_A(infty) / L = pi / (6 beta_eff).
"""
function fetch(
    ::Heisenberg1D, ::EntanglementSaturationDensity, ::Infinite; beta_eff::Real, kwargs...
)
    return fetch(
        Universality(:Heisenberg),
        EntanglementSaturationDensity(),
        Infinite();
        beta_eff=beta_eff,
        kwargs...,
    )
end

# -----------------------------------------------------------------------------
# Lieb-Robinson velocity + entanglement growth slope (des Cloizeaux-Pearson)
# -----------------------------------------------------------------------------

"""
    fetch(m::Heisenberg1D, ::LiebRobinsonVelocity, ::Infinite;
          J = m.J, kwargs...) -> Float64

Maximum group velocity of the des Cloizeaux-Pearson spinon
dispersion of the spin-1/2 Heisenberg chain. With epsilon(k) =
(pi/2) |J sin k| the spinon velocity is

    v_s = pi J / 2,

attained at k = 0 and k = pi. This is the standard CC quasi-particle
velocity that enters the entanglement-spreading formulas.
"""
function fetch(m::Heisenberg1D, ::LiebRobinsonVelocity, ::Infinite; J::Real=1.0, kwargs...)
    return π * abs(J) / 2
end

"""
    fetch(m::Heisenberg1D, ::EntanglementGrowthSlope, ::Infinite;
          beta_eff::Real, kwargs...) -> Float64

Calabrese-Cardy 2005 linear-growth slope of post-quench half-system
entanglement entropy for the gapless Heisenberg chain (c = 1, v = pi J / 2):

    dS_A / dt = pi c v / (3 beta_eff) = pi^2 J / (6 beta_eff).

Delegates to `Universality(:Heisenberg)`.
"""
function fetch(
    m::Heisenberg1D, ::EntanglementGrowthSlope, ::Infinite; beta_eff::Real, kwargs...
)
    return fetch(
        Universality(:Heisenberg),
        EntanglementGrowthSlope(),
        Infinite();
        v=fetch(m, LiebRobinsonVelocity(), Infinite()),
        beta_eff=beta_eff,
        kwargs...,
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Conformal tower of states (SU(2) symmetric point)
# ─────────────────────────────────────────────────────────────────────────────

raw"""
    fetch(model::Heisenberg1D, q::ConformalTower, bc::PBC; J::Real=1.0, kwargs...) -> Vector{NamedTuple}

Conformal tower of states excitation energies of the critical spin-1/2 antiferromagnetic
Heisenberg chain (SU(2)_1 WZW universality, c=1). Delegates to `Universality(:Heisenberg)`
with the des Cloizeaux-Pearson sound velocity `v = π J / 2` and system size `L` extracted
from `bc`.

The two SU(2)_1 primary representations and their scaling dimensions are:
- j = 0 (vacuum):  Δ = 0,   degeneracy = 1
- j = 1/2 (spin):  Δ = 1/4, degeneracy = 4  [(2j+1)² = 4 for both chiral sectors]

(j = 1 is not a WZW primary for k = 1; the constraint is j ≤ k/2 = 1/2.)

The first descendant level at Δ = 1 carries 9 states from 3 left × 3 right SU(2) currents
acting on the j = 0 vacuum.

Only PBC is implemented. OBC raises `ErrorException`.

# Arguments
- `bc::PBC`: periodic boundary condition; system size `L` is read from `bc.N`.
- `J::Real=1.0`: exchange coupling (des Cloizeaux-Pearson velocity v = π|J|/2).

# Returns
`Vector{NamedTuple{(:energy, :dimension, :degeneracy), Tuple{Float64, Float64, Int}}}` —
see `fetch(::Universality{:Heisenberg}, ::ConformalTower, ...)` for full documentation.

# References
- I. Affleck, *Phys. Rev. Lett.* **56**, 746 (1986). — SU(2)_1 WZW spectrum in spin chains.
- J. Cardy, *Nucl. Phys. B* **270**, 186 (1986). — operator content of 1+1D CFTs.
- J. des Cloizeaux, J. J. Pearson, *Phys. Rev.* **128**, 2131 (1962). — spinon velocity v = πJ/2.
"""
function fetch(model::Heisenberg1D, q::ConformalTower, bc::PBC; J::Real=1.0, kwargs...)
    L = _bc_size(bc, kwargs)
    v = π * abs(J) / 2
    return fetch(Universality(:Heisenberg), q, bc; L=L, v=v, J=J, kwargs...)
end

function fetch(model::Heisenberg1D, ::ConformalTower, bc::OBC; J::Real=1.0, kwargs...)
    return error("Heisenberg1D ConformalTower is only implemented for PBC boundary conditions.")
end
