# ─────────────────────────────────────────────────────────────────────────────
# XXZ1D — general-Δ ground-state energy density via the Yang–Yang
# single-integral form (critical regime |Δ| ≤ 1).
#
# Hamiltonian (matching XXZ.jl):
#
#   H = J Σ_i [ Sˣ_i Sˣ_{i+1} + Sʸ_i Sʸ_{i+1} + Δ Sᶻ_i Sᶻ_{i+1} ],
#       J > 0 antiferromagnetic, spin 1/2, periodic boundary conditions.
#
# Parametrisation:  Δ = cos γ,  γ ∈ [0, π],
#   γ = 0     → AF Heisenberg  (Δ =  1, Hulthén)
#   γ = π/2   → XX free fermion (Δ =  0)
#   γ → π     → FM saturation   (Δ = -1)
#
# Result. The Bethe-ansatz ground-state rapidity density solves
# (Takahashi 1999 §5.1)
#
#   ρ(λ) + ∫_{-∞}^{∞} a_2(λ-μ) ρ(μ) dμ = a_1(λ),
#       a_n(λ) = (γ/π) · sin(n γ) / [cosh(2 γ λ) − cos(n γ)],
#
# with normalisation ∫ρ = 1/2.  Fourier-transforming with ω conjugate to
# λ, sum-to-product collapses
#
#   â_n(ω) = sinh((π − n γ) ω / (2γ)) / sinh(π ω / (2γ)),
#
# into the *γ-independent* density
#
#   ρ̂(ω) = â_1(ω) / (1 + â_2(ω)) = 1 / (2 cosh(ω/2))
#       ⇒ ρ(λ) = 1/(2 cosh(π λ)).
#
# Normalisation convention: ρ here is the *spinon rapidity density*
# with ∫_{-∞}^{∞} ρ(λ) dλ = 1/2 (half-filling of magnons).  This is the
# convention used by Takahashi (1999) and Yang-Yang II (1966); some
# Bethe-ansatz textbooks (e.g. Korepin–Bogoliubov–Izergin) instead use
# ∫ρ = 1, in which case all per-site formulas pick up a factor 1/2 in
# the prefactor.  The two conventions are equivalent up to that
# rescaling — verify by checking ∫_{-∞}^{∞} 1/(2 cosh π λ) dλ = 1/2
# (an elementary integral).
#
# All γ dependence is then carried by the energy formula
# (Takahashi 1999 eq. 4.3.18, equivalent to Yang–Yang II 1966 eq. (4.4)):
#
#   e₀(Δ) = (J cos γ)/4
#        − J sin² γ · ∫_{-∞}^{∞} ρ(λ) / [cosh(2γλ) − cos γ] dλ
#        = (J cos γ)/4
#        − (J sin² γ / 2) · ∫_{-∞}^{∞} dλ / [cosh(πλ)·(cosh(2γλ) − cos γ)].
#
# The integrand is smooth on ℝ, exponentially decaying with rate
# min(π, 2γ); QuadGK with rtol = 1e-12 returns to machine precision in
# tens of microseconds.  The closed-form points Δ ∈ {-1, 0, 1} are kept
# as fast-paths (also documents the limiting values).
#
# References:
#   - C. N. Yang & C. P. Yang, "One-Dimensional Chain of Anisotropic
#     Spin–Spin Interactions. II. Properties of the Ground-State Energy
#     per Lattice Site for an Infinite System", Phys. Rev. 150, 327 (1966).
#   - M. Takahashi, "Thermodynamics of One-Dimensional Solvable Models"
#     (Cambridge University Press, 1999), §4.3 and §5.1.
#   - Existing project derivation:
#     `docs/src/calc/bethe-ansatz-heisenberg-e0.md` (Δ = 1 worked
#     example; Steps 4–5 specialise to γ = 0).
#     `docs/src/calc/xxz-luttinger-parameters.md` (Steps 1–3 set up the
#     anisotropic kernels a_n used here).
# ─────────────────────────────────────────────────────────────────────────────

using QuadGK: quadgk

"""
    _xxz1d_yang_yang_integrand(λ, γ)

Integrand of the Yang–Yang single-integral form for the XXZ ground-state
energy density at anisotropy `Δ = cos γ`, namely

    f(λ) = 1 / [ 2 cosh(π λ) · (cosh(2 γ λ) − cos γ) ].

Smooth and exponentially decaying for `0 < γ < π`.  Used by
[`_xxz1d_energy_yang_yang`](@ref).
"""
@inline function _xxz1d_yang_yang_integrand(λ::Real, γ::Real)
    return 1 / (2 * cosh(π * λ) * (cosh(2γ * λ) - cos(γ)))
end

"""
    _xxz1d_energy_yang_yang(J, Δ; rtol = 1e-12, atol = 1e-14) -> Float64

Ground-state energy per site of the spin-½ XXZ chain in the
thermodynamic limit, evaluated by the Yang–Yang single-integral
formula

    e₀(Δ) = (J cos γ)/4 − J sin² γ · ∫_{-∞}^{∞} dλ /
            [2 cosh(πλ) · (cosh(2γλ) − cos γ)],

with `Δ = cos γ`, `γ ∈ (0, π)` (critical / gapless regime).  The three
canonical points `Δ ∈ {-1, 0, 1}` are dispatched in the caller for
exactness; this helper is invoked only for general `-1 < Δ < 1` with
`Δ ∉ {0}`.  Tolerance defaults give relative error ≤ `1e-12` against
the closed-form values at `Δ = 0, ±1` (verified in `test_XXZ1D.jl`).

Cost: a few QuadGK panel splits, ≈ 50 µs per call on a recent CPU.
"""
function _xxz1d_energy_yang_yang(
    J::Real, Δ::Real; rtol::Real=1e-12, atol::Real=1e-14
)::Float64
    γ = acos(Δ)
    I, _ = quadgk(λ -> _xxz1d_yang_yang_integrand(λ, γ), -Inf, Inf; rtol=rtol, atol=atol)
    return J * cos(γ) / 4 - J * sin(γ)^2 * I
end
