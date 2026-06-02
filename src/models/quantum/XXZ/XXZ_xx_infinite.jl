# ─────────────────────────────────────────────────────────────────────────────
# XXZ chain at Δ = 0 (XX / free fermion) — thermodynamic-limit finite-T
# observables.
#
# After the Jordan-Wigner transformation, the XX chain in the *spin*
# convention used throughout XXZ.jl
#
#     H = J Σᵢ (Sˣᵢ Sˣᵢ₊₁ + Sʸᵢ Sʸᵢ₊₁),    Sᵅ = σᵅ / 2
#
# maps (after the cancellation of nearest-neighbour JW strings) to the
# tight-binding fermion chain
#
#     H_JW = (J/2) Σᵢ (c†ᵢ cᵢ₊₁ + h.c.)
#
# whose single-particle dispersion is
#
#     ε(k) = J cos(k),     k ∈ [-π, π]
#
# at half filling and zero chemical potential.  We work with the
# equivalent (k → π + k) form
#
#     ε(k) = -J cos(k)
#
# so that the band minimum sits at k = 0; this is purely a relabelling
# and leaves all thermal observables invariant since they depend on
# ε only through the even functionals
# `log(2 cosh(βε/2))`, `ε · tanh(βε/2)`, `(βε/2)² sech²(βε/2)`.
#
# Per-site free-fermion thermodynamics in the thermodynamic limit
# (1/(2π) ∫_{-π}^{π} = 1/π ∫_0^π by k ↔ -k symmetry):
#
#     f(β)  = -(1/(πβ))  ∫₀^π log(2 cosh(β ε(k) / 2)) dk
#     e(β)  = -(1/(2π))  ∫₀^π ε(k) · tanh(β ε(k) / 2) dk
#     s(β)  =  (1/π)     ∫₀^π [ log(2 cosh(βε/2)) − (βε/2) tanh(βε/2) ] dk
#           ≡  β (e(β) − f(β))
#     C(β)  =  (1/π)     ∫₀^π (β ε(k) / 2)² sech²(β ε(k) / 2) dk
#
# **Why does `e` carry 1/(2π) while `f, s, C` carry 1/π?**  All four
# observables are defined as `(1/(2π)) ∫_{-π}^{π} ⋯ dk`; collapsing
# to [0, π] by the k ↔ -k symmetry of the integrand contributes a
# factor 2 to the overall coefficient.  The remaining factor comes
# from the algebraic identity
#     log(1 + e^{-βε}) = -βε/2 + log(2 cosh(βε/2))
# combined with ∫_{-π}^{π} ε(k) dk = 0 (true for ε(k) = -J cos k since
# cos is odd over a full period of length 2π — equivalently ε is even
# in k but its mean over the BZ is zero).  The linear -βε/2 piece
# integrates to zero, so for `f` (and inherited `s, C`) the
# integrand collapses cleanly to `log(2 cosh(βε/2))` — already
# carrying an implicit factor of 2 from the cosh.  For `e` the
# integrand `ε · n_F(ε)` does not benefit from this doubling: the
# similar identity `ε · n_F(ε) = ε/2 - (ε/2) tanh(βε/2)` produces a
# linear ε/2 piece that *also* integrates to zero, but no extra factor
# of 2 appears in the surviving `-(ε/2) tanh` term.  Hence the k →
# [0, π] reduction of `e` retains the prefactor at 1/(2π).
# A direct algebraic check: at β → ∞,
#     e(∞) = -(1/(2π)) ∫₀^π ε(k) sgn(ε(k)) dk
#          = -(1/(2π)) ∫₀^π |ε(k)| dk
#          = -(J/(2π)) · ∫₀^π |cos k| dk
#          = -(J/(2π)) · 2  =  -J/π,
# matching `_xxz1d_energy_free_fermion(J) = -J/π` from XXZ.jl.
#
# Limit checks:
#
#   • β → ∞:  e(∞) = -J/π          ← XX ground-state energy density
#                                   (Hulthén / Yang-Yang, see XXZ.jl)
#             f(∞) → e(∞)           (entropy → 0)
#             s(∞) → 0
#             C(∞) → 0 (frozen-out)
#   • β → 0:  e(0) → 0             (each integrand picks up only an
#                                    O(β) contribution from the
#                                    high-T expansion of `tanh`)
#             f(0) → -(log 2)/β   (free spin per site)
#             s(0) → log 2          (per-site fermion / spin entropy)
#             C(0) → 0
#
# References
#   - G. D. Mahan, *Many-Particle Physics* (3rd ed.), §1.3.
#   - P. Coleman, *Introduction to Many-Body Physics*, §2.4.
#   - M. Takahashi, *Thermodynamics of One-Dimensional Solvable Models*
#     (1999), §4 — XX point as Δ = 0 / γ = π/2 limit of XXZ.
#
# This file implements only the Δ = 0 special case.  General-Δ TBA
# (issue #108) is a separate workstream that will live alongside the
# existing OBC ED methods.
# ─────────────────────────────────────────────────────────────────────────────

using QuadGK: quadgk

# ── Numerically robust log(2 cosh x) ──────────────────────────────────────────
@inline function _xx_logcosh2(x::Real)
    a = abs(x)
    return a + log1p(exp(-2 * a))
end

# ── Single-particle dispersion (spin convention; band minimum at k = 0) ──────
@inline _xx_dispersion(k::Real, J::Real) = -J * cos(k)

# ─────────────────────────────────────────────────────────────────────────────
# Internal kernel: per-site thermodynamic potential of the XX chain at finite β.
#
# `quantity` ∈ (:free_energy, :energy, :entropy, :specific_heat).
#
# Integrals are evaluated by adaptive Gauss-Kronrod quadrature on [0, π];
# the symmetry k ↔ -k is folded into the prefactor (1/π for the symmetric
# integrands, 1/(2π) for `e`; see header for the derivation).  Each branch
# returns the per-site value in the spin-S convention used by XXZ.jl
# (so the β → ∞ limit of `:energy` reproduces -J/π exactly).
# ─────────────────────────────────────────────────────────────────────────────
function _xx_thermo_infinite(quantity::Symbol, J::Real, β::Real)
    if quantity === :free_energy
        # f(β) = -(1/(πβ)) ∫₀^π log(2 cosh(βε(k)/2)) dk
        integrand = k -> _xx_logcosh2(β * _xx_dispersion(k, J) / 2)
        val, _ = quadgk(integrand, 0.0, π; rtol=1e-10)
        return -val / (π * β)
    elseif quantity === :energy
        # e(β) = -(1/(2π)) ∫₀^π ε(k) tanh(β ε(k) / 2) dk
        integrand = k -> begin
            εk = _xx_dispersion(k, J)
            εk * tanh(β * εk / 2)
        end
        val, _ = quadgk(integrand, 0.0, π; rtol=1e-10)
        return -val / (2 * π)
    elseif quantity === :entropy
        # s(β) = (1/π) ∫₀^π [ log(2 cosh(βε/2)) - (βε/2) tanh(βε/2) ] dk
        integrand = k -> begin
            εk = _xx_dispersion(k, J)
            x = β * εk / 2
            _xx_logcosh2(x) - x * tanh(x)
        end
        val, _ = quadgk(integrand, 0.0, π; rtol=1e-10)
        return val / π
    elseif quantity === :specific_heat
        # C(β) = (1/π) ∫₀^π (βε/2)² sech²(βε/2) dk
        integrand = k -> begin
            εk = _xx_dispersion(k, J)
            x = β * εk / 2
            x^2 * sech(x)^2
        end
        val, _ = quadgk(integrand, 0.0, π; rtol=1e-10)
        return val / π
    else
        error("Unknown XX thermal quantity: $quantity")
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# fetch dispatch — Δ = 0 only.
#
# These methods are guarded by `_xx_is_free_fermion(model)`
# (`isapprox(Δ, 0; atol=1e-12)`); for other Δ they emit a warn+NaN
# placeholder.  Issue #108 will replace those NaNs with a TBA
# implementation.
# ─────────────────────────────────────────────────────────────────────────────

@inline _xx_is_free_fermion(model::XXZ1D) = isapprox(model.Δ, 0.0; atol=1e-12)

function _xx_warn_general_delta(quantity::Symbol, Δ::Real)
    @warn (
        "XXZ1D thermal observable at Infinite() for general Δ requires the " *
        "thermal Bethe ansatz (issue #108); only Δ = 0 (XX / free fermion) " *
        "is implemented. Returning NaN."
    ) quantity = quantity Δ = Δ
    return NaN
end

# ── Energy{:per_site}, Infinite — finite-T (β provided) and ground state ─────
#
# Replaces the `Energy{:per_site}, Infinite` method previously defined
# in XXZ.jl by branching on whether `beta` was passed:
#
#   • no `beta` kwarg ⇒ closed-form ground-state value (canonical
#                       three points; warn+NaN otherwise) — identical
#                       to the XXZ.jl implementation it supersedes.
#   • `beta` kwarg     ⇒ free-fermion finite-T integral at Δ = 0;
#                        warn+NaN at other Δ (issue #108).

"""
    fetch(model::XXZ1D, ::Energy{:per_site}, ::Infinite; [beta]) -> Float64

Per-site energy of the infinite XXZ chain.

* Without `beta`: ground-state energy density.  Closed-form values at
  the three canonical points `Δ ∈ {-1, 0, 1}`; warns and returns
  `NaN` otherwise (general-Δ Bethe ansatz is a follow-up).
* With `beta`: thermal energy density `⟨H⟩_β / N`.  At `Δ = 0`
  evaluated by Gauss-Kronrod quadrature of the free-fermion integral

      e(β) = -(1/(2π)) ∫₀^π ε(k) tanh(β ε(k) / 2) dk,   ε(k) = -J cos k.

  The β → ∞ limit reproduces `-J/π` (matching XXZ.jl's
  `_xxz1d_energy_free_fermion`).  At general Δ the thermal Bethe
  ansatz (issue #108) is required; a warning is emitted and `NaN`
  returned.
"""
function fetch(model::XXZ1D, ::Energy{:per_site}, ::Infinite; kwargs...)
    if haskey(kwargs, :beta)
        β = kwargs[:beta]
        if _xx_is_free_fermion(model)
            return _xx_thermo_infinite(:energy, model.J, β)
        end
        return _xx_warn_general_delta(:energy, model.Δ)
    end
    # Ground state — same logic as the original XXZ.jl method.
    J, Δ = model.J, model.Δ
    if isapprox(Δ, 0.0; atol=1e-12)
        return _xxz1d_energy_free_fermion(J)
    elseif isapprox(Δ, 1.0; atol=1e-12)
        return _xxz1d_energy_heisenberg_af(J)
    elseif isapprox(Δ, -1.0; atol=1e-12)
        return _xxz1d_energy_heisenberg_fm(J)
    elseif -1.0 < Δ < 1.0
        return _xxz1d_energy_yang_yang(J, Δ)
    else
        @warn "XXZ1D Energy: gapped regime |Δ| > 1 not yet implemented; " *
            "use OBC dense ED at small N for a finite-size reference." Δ = Δ
        return NaN
    end
end

# ── FreeEnergy, ThermalEntropy, SpecificHeat at Infinite ────────────────────

"""
    fetch(model::XXZ1D, ::FreeEnergy, ::Infinite; beta::Real, kwargs...)

Per-site Helmholtz free energy of the infinite XXZ chain.
Currently only `Δ = 0` (XX / free fermion) is implemented:

    f(β) = -(1/(πβ)) ∫₀^π log(2 cosh(β ε(k) / 2)) dk,   ε(k) = -J cos k.

For other Δ the thermal Bethe ansatz (issue #108) is required; this
method emits a warning and returns `NaN`.

References: Mahan, *Many-Particle Physics*, §1.3; Coleman, *Introduction
to Many-Body Physics*, §2.4; Takahashi (1999), §4.
"""
function fetch(model::XXZ1D, ::FreeEnergy, ::Infinite; beta::Real, kwargs...)
    isempty(kwargs) || @warn(
        "fetch(XXZ1D, FreeEnergy, Infinite) received unrecognized kwargs; they are ignored. " *
            "Magnetic field h ≠ 0 is not yet wired through the Klümper NLIE (issue #521).",
        kwargs=collect(keys(kwargs))
    )
    if _xx_is_free_fermion(model)
        return _xx_thermo_infinite(:free_energy, model.J, beta)
    end
    if -1 < model.Δ < 1
        e0 = fetch(model, Energy{:per_site}(), Infinite())
        excess = _xxz_klumper_free_energy_excess(model, beta)
        return e0 + excess
    end
    return _xx_warn_general_delta(:free_energy, model.Δ)
end

"""
    fetch(model::XXZ1D, ::ThermalEntropy, ::Infinite; beta::Real, kwargs...)

Per-site Gibbs entropy of the infinite XXZ chain.  Implemented at
Δ = 0 only via

    s(β) = (1/π) ∫₀^π [ log(2 cosh(βε/2)) - (βε/2) tanh(βε/2) ] dk,
    ε(k) = -J cos k,

equivalent to `s(β) = β (e(β) - f(β))`.  In the high-T limit
(`β → 0`) this saturates to `log 2` per site.  Returns `NaN`
(with a warning) for general Δ pending issue #108.
"""
function fetch(model::XXZ1D, ::ThermalEntropy, ::Infinite; beta::Real, kwargs...)
    isempty(kwargs) || @warn(
        "fetch(XXZ1D, ThermalEntropy, Infinite) received unrecognized kwargs; they are ignored.",
        kwargs=collect(keys(kwargs))
    )
    if _xx_is_free_fermion(model)
        return _xx_thermo_infinite(:entropy, model.J, beta)
    end
    if -1 < model.Δ < 1
        return _xxz_klumper_entropy(model, beta)
    end
    return _xx_warn_general_delta(:entropy, model.Δ)
end

"""
    fetch(model::XXZ1D, ::SpecificHeat, ::Infinite; beta::Real, kwargs...)

Per-site heat capacity of the infinite XXZ chain.  At Δ = 0,

    C(β) = (1/π) ∫₀^π (β ε / 2)² sech²(β ε / 2) dk,   ε(k) = -J cos k.

For `-1 < Δ < 1` (Δ ≠ 0) routes through the Klümper NLIE (issue #521).
See `_xxz_klumper_specific_heat` for the finite-difference details.
"""
function fetch(model::XXZ1D, ::SpecificHeat, ::Infinite; beta::Real, kwargs...)
    isempty(kwargs) || @warn(
        "fetch(XXZ1D, SpecificHeat, Infinite) received unrecognized kwargs; they are ignored.",
        kwargs=collect(keys(kwargs))
    )
    if _xx_is_free_fermion(model)
        return _xx_thermo_infinite(:specific_heat, model.J, beta)
    end
    if -1 < model.Δ < 1
        return _xxz_klumper_specific_heat(model, beta)
    end
    return _xx_warn_general_delta(:specific_heat, model.Δ)
end
