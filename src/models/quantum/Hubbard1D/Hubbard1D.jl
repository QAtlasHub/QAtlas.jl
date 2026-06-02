# ─────────────────────────────────────────────────────────────────────────────
# Hubbard1D — 1D Hubbard model, Lieb–Wu Bethe ansatz half-filling closed form.
#
# Hamiltonian:
#
#   H = -t Σ_{i, σ} (c†_{i,σ} c_{i+1,σ} + h.c.)
#       + U Σ_i n_{i,↑} n_{i,↓}
#       - μ Σ_i n_i
#
#   (spin-½ fermions; t > 0 hopping, U > 0 on-site repulsion, μ chemical
#   potential — half filling corresponds to μ = U/2 by particle-hole
#   symmetry.)
#
# Phase 1 (this file) implements the Lieb–Wu (1968) closed-form
# integrals at half filling only:
#
#   E₀/N = -4 t ∫_0^∞ dω  J₀(ω) J₁(ω) / [ ω (1 + exp(ω U / 2t)) ]
#
#   Δ_c  = (16 t² / U) ∫_1^∞ dω  √(ω² - 1) / sinh(2π t ω / U)
#
#   Δ_s  = 0                     (gapless spinons; rigorous Lieb–Wu)
#
# Asymptotic limits used as test anchors:
#
#   U/t → 0   (free fermion):    E₀/N → -4 t / π,        Δ_c → 0
#   U/t → ∞   (Heisenberg AFM):  E₀/N → -4 t² log 2 / U, Δ_c → U - 4t + 8 t² log 2 / U
#
# Phase 2 (deferred): general filling (μ ≠ U/2) and finite-T QTM —
# both require the full Lieb–Wu coupled integral equations.
#
# References:
#
#   - E. H. Lieb, F. Y. Wu, "Absence of Mott transition in an exact
#     solution of the short-range, one-band model in one dimension",
#     Phys. Rev. Lett. 20, 1445 (1968).
#   - F. H. L. Essler, H. Frahm, F. Göhmann, A. Klümper, V. E. Korepin,
#     "The One-Dimensional Hubbard Model", Cambridge University Press
#     (2005) — the canonical textbook.
#   - E. H. Lieb, F. Y. Wu, "The one-dimensional Hubbard model: a
#     reminiscence", Physica A 321, 1 (2003).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Fermion bilinears c†c
#   Observable:  Fermion (number n = c†c, bilinear ⟨c†_i c_j⟩); derived spin observables follow spin S = σ/2
#   Reference:   docs/src/conventions.md §Fermion convention

using QuadGK: quadgk
using SpecialFunctions: besselj0, besselj1

"""
    Hubbard1D(; t::Real = 1.0, U::Real = 4.0, μ::Real = 2.0) <: AbstractQAtlasModel

1D Hubbard model

    H = -t Σ_{i, σ} (c†_{i,σ} c_{i+1,σ} + h.c.)
        + U Σ_i n_{i,↑} n_{i,↓}
        - μ Σ_i n_i

Convention: `t > 0` hopping, `U > 0` on-site repulsion.  Half filling
corresponds to `μ = U/2` (particle-hole-symmetric point).

# Phase 1 scope (this release)

Phase 1 implements the Lieb–Wu (1968) closed-form integrals at **half
filling only** (`μ = U/2`).  The supported `Infinite()` quantities are

- [`GroundStateEnergyDensity`](@ref) — `E₀/N` via the Lieb–Wu integral,
- [`ChargeGap`](@ref) — Mott gap `Δ_c` via the Lieb–Wu second integral,
- [`SpinGap`](@ref) — `0` exactly (rigorous gapless-spinon result).

Off half filling raises a `DomainError` — the general-filling
(coupled-integral-equation) and finite-temperature (QTM) surfaces are
deferred to Phase 2.
"""
struct Hubbard1D <: AbstractQAtlasModel
    t::Float64
    U::Float64
    μ::Float64
end
function Hubbard1D(; t::Real=1.0, U::Real=4.0, μ::Real=2.0)
    return Hubbard1D(Float64(t), Float64(U), Float64(μ))
end

# ─── Half-filling guard ───────────────────────────────────────────────

"""
    _hubbard1d_check_half_filling(model::Hubbard1D)

Throw a `DomainError` if `model.μ` is not the half-filling chemical
potential `U/2`.  Phase 1 only exposes the Lieb–Wu closed forms at the
particle-hole-symmetric point; general filling requires the coupled
integral equations and is tracked as Phase 2 follow-up.
"""
function _hubbard1d_check_half_filling(model::Hubbard1D)
    half = model.U / 2
    if !isapprox(model.μ, half; atol=1e-12, rtol=1e-12)
        throw(
            DomainError(
                model.μ,
                "Hubbard1D Phase 1: only half-filling (μ = U/2) implemented; " *
                "got μ=$(model.μ) vs U/2=$(half)",
            ),
        )
    end
    return nothing
end

# ─── Lieb–Wu integrals (half filling) ─────────────────────────────────

"""
    _hubbard1d_e0_integrand(ω, t, U)

Integrand of the Lieb–Wu ground-state-energy integral at half filling,

    f(ω) = J₀(ω) J₁(ω) / [ ω (1 + exp(ω U / 2t)) ].

The integrand is regular at `ω = 0` (`J₁(ω)/ω → 1/2` as `ω → 0`); we
evaluate the limit explicitly to avoid 0/0 at the lower endpoint that
some `quadgk` adaptive subdivisions produce.
"""
function _hubbard1d_e0_integrand(ω::Float64, t::Float64, U::Float64)::Float64
    if ω == 0.0
        # J₀(0) = 1, J₁(ω)/ω → 1/2, denominator = 1 + 1 = 2 → integrand → 1/4.
        return 0.25
    end
    return besselj0(ω) * besselj1(ω) / (ω * (1.0 + exp(ω * U / (2.0 * t))))
end

"""
    _hubbard1d_e0(t, U) -> Float64

Lieb–Wu (1968) ground-state energy per site at half filling:

    E₀/N = -4 t ∫_0^∞ dω  J₀(ω) J₁(ω) / [ ω (1 + exp(ω U / 2t)) ].
"""
function _hubbard1d_e0(t::Float64, U::Float64)::Float64
    integral, _ = quadgk(
        ω -> _hubbard1d_e0_integrand(ω, t, U), 0.0, Inf; rtol=1e-12, atol=1e-14
    )
    return -4.0 * t * integral
end

"""
    _hubbard1d_charge_gap(t, U) -> Float64

Lieb–Wu (1968) Mott (charge) gap at half filling:

    Δ_c = (16 t² / U) ∫_1^∞ dω  √(ω² - 1) / sinh(2π t ω / U).

The integrand vanishes at `ω = 1` like `√(ω-1)` and decays
exponentially at large `ω`; `quadgk` resolves both endpoints with the
default Gauss–Kronrod adaptive subdivision.
"""
function _hubbard1d_charge_gap(t::Float64, U::Float64)::Float64
    g(ω) = sqrt(ω^2 - 1.0) / sinh(2.0 * π * t * ω / U)
    integral, _ = quadgk(g, 1.0, Inf; rtol=1e-12, atol=1e-14)
    return (16.0 * t^2 / U) * integral
end

# ─── fetch methods ────────────────────────────────────────────────────

"""
    fetch(model::Hubbard1D, ::GroundStateEnergyDensity, ::Infinite) -> Float64

Lieb–Wu ground-state energy density `E₀/N` at half filling.  Currently
only implemented for `μ = U/2`; off-half-filling raises a
`DomainError`.

# Asymptotic limits

- `U/t → 0`:  `E₀/N → -4t/π`  (free 1D fermion).
- `U/t → ∞`: `E₀/N → -4 t² log 2 / U`  (Heisenberg AFM reduction).

Both are exercised in `test/standalone/test_hubbard1d.jl`.

# References

- Lieb–Wu, *PRL* **20**, 1445 (1968).
- Essler et al., *The One-Dimensional Hubbard Model* (Cambridge, 2005).
"""
function fetch(model::Hubbard1D, ::GroundStateEnergyDensity, ::Infinite; kwargs...)
    _hubbard1d_check_half_filling(model)
    return _hubbard1d_e0(model.t, model.U)
end

"""
    fetch(model::Hubbard1D, ::ChargeGap, ::Infinite) -> Float64

Lieb–Wu Mott (charge) gap at half filling:

    Δ_c = (16 t² / U) ∫_1^∞ dω  √(ω² - 1) / sinh(2π t ω / U).

Strictly positive for any `U > 0` (no Mott transition in 1D — the
chain is insulating at half filling for arbitrarily small `U`, the
celebrated Lieb–Wu result).

# Asymptotic limits

- `U → 0`: `Δ_c → 0` (exponentially small, ∝ exp(-2π t / U)).
- `U → ∞`: `Δ_c → U - 4t + 8 t² log 2 / U`.

# References

- Lieb–Wu, *PRL* **20**, 1445 (1968).
"""
function fetch(model::Hubbard1D, ::ChargeGap, ::Infinite; kwargs...)
    _hubbard1d_check_half_filling(model)
    return _hubbard1d_charge_gap(model.t, model.U)
end

"""
    fetch(model::Hubbard1D, ::SpinGap, ::Infinite) -> Float64

Spin gap of the half-filled 1D Hubbard chain.  Returns `0.0` exactly:
the spinon branch is gapless for any `U > 0` (rigorous Lieb–Wu
result; spin-charge separation in the low-energy effective theory).

# References

- Lieb–Wu, *PRL* **20**, 1445 (1968).
- Essler et al., *The One-Dimensional Hubbard Model* (Cambridge, 2005).
"""
function fetch(model::Hubbard1D, ::SpinGap, ::Infinite; kwargs...)
    _hubbard1d_check_half_filling(model)
    return 0.0
end

# ═══════════════════════════════════════════════════════════════════════════════
# Luttinger parameter at U=0 (free-fermion limit, Phase 2)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::Hubbard1D, ::LuttingerParameter, ::Infinite; t=m.t, U=m.U, μ=m.μ) -> Float64

Luttinger-liquid parameter of the 1D Hubbard model in the **free-fermion
limit** U = 0:

    K = 1               (non-interacting spinful fermions; Voit 1995)

For U > 0, the Lieb–Wu Bethe-ansatz solution gives a non-closed-form
expression for both K_ρ (charge) and K_σ (spin); deferred to Phase 2
(or Phase 3).  This entry exposes only the U = 0 free-fermion limit
and throws `DomainError` for any U ≠ 0.

# References

- E. H. Lieb, F. Y. Wu, *Phys. Rev. Lett.* **20**, 1445 (1968).
- J. Voit, *Rep. Prog. Phys.* **58**, 977 (1995) — TLL review for Hubbard.
"""
function fetch(
    m::Hubbard1D,
    ::LuttingerParameter,
    ::Infinite;
    t::Real=m.t,
    U::Real=m.U,
    μ::Real=m.μ,
    kwargs...,
)
    t > 0 ||
        throw(DomainError(t, "Hubbard1D LuttingerParameter requires t > 0; got t = $t."))
    if !iszero(U)
        throw(
            DomainError(
                U,
                "Hubbard1D LuttingerParameter: U ≠ 0 requires Lieb–Wu Bethe-ansatz K_ρ(U), " *
                "K_σ(U) integrals (Voit 1995) — non-closed-form, deferred to Phase 2. Got U = $U.",
            ),
        )
    end
    return 1.0       # free-fermion limit
end

# ═══════════════════════════════════════════════════════════════════════════════
# Free energy at finite T via JKS NLIE (Stage D.2 wire-up of #523)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::Hubbard1D, ::FreeEnergy, ::Infinite; beta, kwargs...) -> Float64

Per-site Helmholtz free energy of the 1D Hubbard chain at finite
temperature 1/beta, computed by the JKS 1998 quantum-transfer-matrix
non-linear integral equations (issue #523).

The solver uses the paper-precise eq (47) NLIE in 3 channels (b, c, c̄)
on a discretised contour, then evaluates the QTM eigenvalue via paper
eq (49) third form. At β → 0 the result reproduces the atomic limit
`atomic_free_energy(β, U, μ)` to within a few percent (Stage C.10
test guards this).

# Keyword arguments

- `beta::Real`: inverse temperature (β > 0 required).
- `H::Real = 0`: external magnetic field (couples to magnetisation).
- `grid_N::Int = 64`, `x_max::Real = 8.0`: discretisation knobs.
- `alpha::Real = m.U / 6`: contour shift in the b channel
  (0 < α < η = U/4).
- `tol::Real = 1e-6`, `maxiter::Int = 40`: Newton convergence knobs.

# Notes

Half-filling (μ = U/2) is the well-tested regime. Other fillings
work as long as the NLIE converges (the Newton solver does not warn
on non-convergence; the result is `NaN` in that case).

# References

- JKS 1998 = Jüttner, Klümper, Suzuki, *Nucl. Phys. B* **522**, 471 (1998),
  arXiv:cond-mat/9711310.
- Paper eqs (23), (47), (48), (49), (54), (55) all implemented per PDF
  (Stage C.22c paper-precise rewrite + Stage C.24 FE evaluator fix).
"""
function fetch(
    m::Hubbard1D,
    ::FreeEnergy,
    ::Infinite;
    beta::Real,
    H::Real=0.0,
    grid_N::Int=128,
    x_max::Real=32.0,
    alpha::Real=m.U / 6,
    tol::Real=1e-6,
    maxiter::Int=40,
    kwargs...,
)
    beta > 0 ||
        throw(DomainError(beta, "Hubbard1D FreeEnergy@Infinite requires β > 0; got β = $(beta)."))
    isempty(kwargs) || @warn(
        "fetch(Hubbard1D, FreeEnergy, Infinite) received unrecognized kwargs; ignored.",
        kwargs=collect(keys(kwargs))
    )
    return Hubbard1DJKSNLIE.hubbard1d_jks_free_energy(
        m.t, m.U, m.μ, beta;
        H=H, grid_N=grid_N, x_max=x_max, alpha=alpha, tol=tol, maxiter=maxiter,
        solver=:full_newton,
    )
end
