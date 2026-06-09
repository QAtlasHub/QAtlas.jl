# ──────────────────────────────────────────────────────────────────────────────
# Kitaev1D — free-fermion (BdG) finite-temperature thermodynamics.
#
# The Kitaev chain is a free-fermion (p-wave BdG) system with bulk quasiparticle
# dispersion E(k) = √((2t cos k + μ)² + 4Δ² sin²k) (`_kitaev1d_dispersion`, in
# Kitaev1D.jl).  The thermodynamic-limit potentials per site are the standard
# BdG free-fermion integrals (cf. TFIM, the Jordan–Wigner dual):
#
#    f  = -(1/2π) ∫_{-π}^{π} [ E/2 + β⁻¹ ln(1 + e^{-βE}) ] dk      (FreeEnergy)
#    ε  = -(1/2π) ∫_{-π}^{π} (E/2) tanh(βE/2) dk                   (Energy, Kitaev1D.jl)
#   c_v = (β²/8π) ∫_{-π}^{π} E² sech²(βE/2) dk = β² Var(E)/L       (energy FDT form)
#    s  = β(ε − f)
#
# Limits pin the conventions: high-T s → ln 2 (a spinless mode has two states);
# T → 0 ε → ε₀ (ground state, β → ∞), s → 0.  src keeps ForwardDiff in [extras],
# so these are explicit closed-form integrals (QuadGK); the test cross-checks
# c_v == -β² ∂ε/∂β (AutoDiff) and against `fd_free_fermion_thermo` (#676).
# ──────────────────────────────────────────────────────────────────────────────

using QuadGK: quadgk

@inline function _kitaev1d_require_beta(beta::Real)
    return beta > 0 || throw(
        DomainError(beta, "Kitaev1D finite-T thermodynamics require β > 0; got β = $beta.")
    )
end

"""
    fetch(model::Kitaev1D, ::FreeEnergy, ::Infinite; beta) -> Float64

Helmholtz free energy per site `f(β) = -β⁻¹ log Z/L` of the infinite Kitaev
chain (BdG free fermions):

    f(β) = -(1/2π) ∫_{-π}^{π} [ E(k)/2 + β⁻¹ ln(1 + e^{-βE(k)}) ] dk.

`β → ∞` recovers the ground-state energy per site.
"""
function fetch(model::Kitaev1D, ::FreeEnergy, ::Infinite; beta::Real, kwargs...)
    _kitaev1d_require_beta(beta)
    μ, t, Δ = model.μ, model.t, model.Δ
    result, _ = quadgk(-π, π; rtol=1e-10) do k
        E = _kitaev1d_dispersion(k, μ, t, Δ)
        return -(E / 2 + log1p(exp(-beta * E)) / beta)
    end
    return result / (2π)
end

"""
    fetch(model::Kitaev1D, ::SpecificHeat, ::Infinite; beta) -> Float64

Specific heat per site of the infinite Kitaev chain, the energy
fluctuation–dissipation form for free BdG quasiparticles:

    c_v(β) = (β²/8π) ∫_{-π}^{π} E(k)² sech²(βE(k)/2) dk
           = β² ∫_{-π}^{π} (dk/2π) E(k)² f_k (1 - f_k),   f_k = 1/(e^{βE}+1).
"""
function fetch(model::Kitaev1D, ::SpecificHeat, ::Infinite; beta::Real, kwargs...)
    _kitaev1d_require_beta(beta)
    μ, t, Δ = model.μ, model.t, model.Δ
    result, _ = quadgk(-π, π; rtol=1e-10) do k
        E = _kitaev1d_dispersion(k, μ, t, Δ)
        return E^2 * sech(beta * E / 2)^2
    end
    return beta^2 * result / (8π)
end

"""
    fetch(model::Kitaev1D, ::ThermalEntropy, ::Infinite; beta) -> Float64

Entropy per site `s(β) = β(ε − f)` of the infinite Kitaev chain.  Bounded
between 0 (`T → 0`) and `ln 2` (`T → ∞`: two states per spinless mode).
"""
function fetch(model::Kitaev1D, ::ThermalEntropy, ::Infinite; beta::Real, kwargs...)
    _kitaev1d_require_beta(beta)
    ε = fetch(model, Energy(:per_site), Infinite(); beta=beta)
    f = fetch(model, FreeEnergy(), Infinite(); beta=beta)
    return beta * (ε - f)
end
