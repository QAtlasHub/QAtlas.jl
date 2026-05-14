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
| [`SpontaneousMagnetization`](@ref)| `Infinite` | fixed-point iter      |

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

# Bisection solver for the self-consistency equation m = tanh(beta*J*m).
# Returns the nontrivial positive root for beta*J > 1 and 0 for beta*J <= 1.
#
# Why not Picard / fixed-point iteration: the map m -> tanh(beta*J*m) has
# derivative beta*J*sech^2(beta*J*m_star) at the stable root, and this
# derivative approaches 1 as beta*J -> 1+ (where m_star -> 0).  Near
# criticality the Picard iterate decays as q^k with q -> 1, so machine
# precision needs millions of steps - and at maxiter=400 the iterate sits
# at ~1e-4 relative error, breaking the Landau exponent test (atol = 5e-3
# on m* / sqrt(3 t) at t = 1 - T/T_c = 1e-4).
#
# Bisection on g(m) := m - tanh(beta*J*m) sidesteps the rate degeneracy:
#   g(0) = 0 (unstable trivial root),
#   g has a unique minimum at m_min = atanh(sqrt(1 - 1/beta*J)) / (beta*J)
#     with g(m_min) < 0,
#   g(m) -> 1 - tanh(beta*J) > 0 as m -> 1-.
# So [m_min, 1 - 1e-15] brackets the positive root m_star with g(lo) < 0
# and g(hi) > 0, and ~52 halvings reach 1e-15 relative width.
function _curie_weiss_solve_m(beta_J::Real; tol::Real=1e-14, maxiter::Int=200)
    if beta_J <= 1
        return 0.0
    end
    # Lower bracket: m_min, the unique minimum of g(m) = m - tanh(beta_J*m).
    m_min = atanh(sqrt(1 - 1 / beta_J)) / beta_J
    lo = nextfloat(m_min)
    hi = 1.0 - 1e-15
    g_lo = lo - tanh(beta_J * lo)
    g_hi = hi - tanh(beta_J * hi)
    # In the extreme beta_J -> 1+ limit, floating-point cancellation can
    # push m_min above the true root by an ulp; in that case fall back to
    # a tiny positive seed where g is guaranteed negative.
    if !(g_lo < 0)
        lo = 1e-300
        g_lo = lo - tanh(beta_J * lo)
    end
    if !(g_hi > 0)
        return hi
    end
    for _ in 1:maxiter
        mid = 0.5 * (lo + hi)
        g_mid = mid - tanh(beta_J * mid)
        if g_mid == 0 || (hi - lo) <= tol * max(1.0, hi)
            return mid
        end
        if g_mid < 0
            lo, g_lo = mid, g_mid
        else
            hi, g_hi = mid, g_mid
        end
    end
    return 0.5 * (lo + hi)
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
root is found by bisection of `g(m) := m - tanh(?J m)` on
`[m_min, 1 - eps)` where `m_min = atanh(?(1 - 1/?J)) / ?J`.  Bisection
converges absolutely (no rate degeneracy at criticality) in ~52 steps.

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

# ═══════════════════════════════════════════════════════════════════════════════
# Critical exponents — mean-field (delegate to MeanField universality)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::CurieWeissIsing, ::CriticalExponents, ::Infinite; kwargs...) -> NamedTuple

Mean-field (Landau-Ginzburg) critical exponents of the Curie-Weiss
Ising / complete-graph Ising universality, delegated to
[`MeanField`](@ref):

    α = 0,  β = 1/2,  γ = 1,  δ = 3,  ν = 1/2,  η = 0.

Hyperscaling (Rushbrooke `α + 2β + γ = 2`; Widom `γ = β(δ − 1)`) and
the upper-critical dimension d_c = 4 (above which mean-field is exact)
are encoded in these values.

# References

- L. D. Landau, *Phys. Z. Sowjet.* **11**, 26 (1937).
- H. E. Stanley, *Introduction to Phase Transitions and Critical Phenomena* (1971).
"""
function fetch(::CurieWeissIsing, ::CriticalExponents, ::Infinite; kwargs...)
    return QAtlas.fetch(QAtlas.MeanField(), CriticalExponents())
end
