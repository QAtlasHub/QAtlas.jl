# ─────────────────────────────────────────────────────────────────────────────
# CurieWeissIsing — classical Ising model on the complete graph (mean field).
#
# Hamiltonian (in the thermodynamic-limit normalisation):
#
#     H = -(J/N) Σ_{i<j} σ_i σ_j,   σ_i ∈ {-1, +1}.
#
# With k_B = 1 this is the textbook Curie-Weiss / mean-field Ising
# model.  The saddle-point free energy per site (Landau-Lifshitz §149) is
#
#     f(β, m) = J m² / 2 - β⁻¹ log[2 cosh(β J m)],
#
# whose stationary points are the solutions of the self-consistency equation
#
#     m = tanh(β J m).
#
# For T ≥ T_c = J the only root is m = 0; for T < T_c there is an
# additional pair of nontrivial ±m*(T) roots that spontaneously break
# the ℤ₂ flip symmetry.
#
# The complete-graph geometry is the canonical *upper-critical
# realisation* of the mean-field exponents already exposed in the
# Universality(:MeanField) entry; this file provides the *model side*
# whose thermodynamics matches those exponents exactly.
#
# References:
#   - L. D. Landau, E. M. Lifshitz, *Statistical Physics* §149.
#   - The "Curie-Weiss" naming is a later attribution of the mean-field
#     Ising model after Weiss 1907 and Curie 1907.
# ─────────────────────────────────────────────────────────────────────────────

"""
    CurieWeissIsing(; J::Real = 1.0) <: AbstractQAtlasModel

Classical Ising model on the complete graph (mean field) with
saddle-point Hamiltonian

    H = -(J/N) Σ_{i<j} σ_i σ_j.

The 1/N normalisation makes the energy extensive and the thermodynamic
limit well-defined.  With k_B = 1 the model has critical temperature
`T_c = J` and exhibits the mean-field critical exponents `(β, γ, δ, ν)
= (1/2, 1, 3, 1/2)` already exposed by `Universality(:MeanField)`.

Quantities registered:

| Quantity                          | BC         | Method                |
| --------------------------------- | ---------- | --------------------- |
| [`CriticalTemperature`](@ref)     | `Infinite` | analytic              |
| [`SpontaneousMagnetization`](@ref)| `Infinite` | Newton root solve     |

# References

- L. D. Landau, E. M. Lifshitz, *Statistical Physics* §149.
"""
struct CurieWeissIsing <: AbstractQAtlasModel
    J::Float64
end
CurieWeissIsing(; J::Real=1.0) = CurieWeissIsing(Float64(J))

# ═══════════════════════════════════════════════════════════════════════════════
# Critical temperature
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::CurieWeissIsing, ::CriticalTemperature, ::Infinite; J=m.J) -> Float64

Mean-field critical temperature `T_c = J` (with k_B = 1) for the
Curie-Weiss Ising model.  For `J ≤ 0` the model is paramagnetic /
antiferromagnetic on a complete graph (no ferromagnetic order at any
positive temperature) and `T_c` is returned as `0`.

# References

- L. D. Landau, E. M. Lifshitz, *Statistical Physics* §149.
"""
function fetch(m::CurieWeissIsing, ::CriticalTemperature, ::Infinite; J::Real=m.J)
    return J > 0 ? Float64(J) : 0.0
end

# ═══════════════════════════════════════════════════════════════════════════════
# Spontaneous magnetisation (zero field)
# ═══════════════════════════════════════════════════════════════════════════════

# Newton solver for the self-consistency equation m = tanh(βJ m).
# Returns the nontrivial positive root for βJ > 1 and 0 for βJ ≤ 1
# (paramagnetic phase).
function _curie_weiss_solve_m(βJ::Real; tol::Real=1e-14, maxiter::Int=200)
    if βJ ≤ 1
        return 0.0
    end
    # Initial guess: Landau expansion near T_c gives
    #   m² ≈ 3 (βJ - 1) / (βJ)³  ⇒  m₀ ≈ √(3 (βJ - 1)) / (βJ)^{3/2}.
    m = sqrt(3 * (βJ - 1)) / βJ^1.5
    m = clamp(m, 1e-12, 1.0 - 1e-15)
    for _ in 1:maxiter
        s = tanh(βJ * m)
        # f(m) = tanh(βJ m) - m,  f'(m) = βJ (1 - s²) - 1.
        f = s - m
        fp = βJ * (1 - s^2) - 1
        # Near βJ → 1 from above, fp → 0; protect against degenerate step.
        if abs(fp) < 1e-300
            break
        end
        m_new = clamp(m - f / fp, 0.0, 1.0)
        if abs(m_new - m) ≤ tol * max(1.0, abs(m_new))
            return m_new
        end
        m = m_new
    end
    return m
end

"""
    fetch(::CurieWeissIsing, ::SpontaneousMagnetization, ::Infinite;
          beta::Real, J=m.J) -> Float64

Spontaneous magnetisation `m*(β)` of the Curie-Weiss Ising model at
inverse temperature `β` and zero external field, defined as the
ℤ₂-positive nontrivial root of the self-consistency equation

    m = tanh(β J m).

For `T ≥ T_c = J` (equivalently `β J ≤ 1`) only the trivial `m = 0`
root exists and `0.0` is returned.  For `T < T_c` the nontrivial
root is found by Newton iteration from the Landau-expansion seed
`m₀ ≈ √(3 (βJ - 1)) / (βJ)^{3/2}`, which converges quadratically
across the full ordered phase.

# References

- L. D. Landau, E. M. Lifshitz, *Statistical Physics* §149.
"""
function fetch(
    m::CurieWeissIsing,
    ::SpontaneousMagnetization,
    ::Infinite;
    beta::Real,
    J::Real=m.J,
    kwargs...,
)
    beta > 0 || throw(
        DomainError(
            beta,
            "CurieWeissIsing SpontaneousMagnetization requires β > 0; got β = $beta.",
        ),
    )
    J > 0 || return 0.0     # AFM or zero coupling: no FM order on the complete graph.
    return _curie_weiss_solve_m(beta * J)
end
