# ─────────────────────────────────────────────────────────────────────────────
# CurieWeissIsing — classical Ising model on the complete graph (mean field).
#
# Hamiltonian (thermodynamic-limit normalisation, k_B = 1):
#
#     H = -(J/N) Σ_{i<j} σ_i σ_j  -  h Σ_i σ_i,   σ_i ∈ {-1, +1}.
#
# The 1/N normalisation keeps the energy extensive.  The saddle-point free
# energy per site (Landau-Lifshitz §149) is
#
#     f(β, m; J, h) = J m² / 2  -  β⁻¹ log[2 cosh(β (J m + h))],
#
# whose stable stationary point is the unique solution of
#
#     m = tanh(β (J m + h))
#
# on the same side as `h` (for `h = 0`, both ±m* and 0 are stationary;
# the stable ones are ±m*(T) for T < T_c, only m = 0 for T ≥ T_c = J).
#
# All thermodynamic quantities below use the equilibrium magnetisation
# m*(β, J, h) returned by `_curie_weiss_solve_m(βJ, βh)`:
#
#     u   = -J m*²/2 - h m*
#     f   =  J m*²/2 - β⁻¹ log(2 cosh(β(J m* + h)))
#     s   = log(2 cosh(β(J m* + h))) - β(J m* + h) m*           (Gibbs)
#     c_v = β² (J m* + h)² (1 - m*²) / (1 - β J (1 - m*²))
#     χ   = β (1 - m*²) / (1 - β J (1 - m*²))
#
# `SpontaneousMagnetization` is by definition the h → 0⁺ limit and so
# always uses the zero-field branch independent of the model's `h`
# field.  `CriticalTemperature` returns `J` (the h = 0 transition
# temperature); strictly there is no sharp transition at h ≠ 0, but the
# h = 0 value is the natural reference.
#
# References:
#   - L. D. Landau, E. M. Lifshitz, *Statistical Physics* §149.
#   - The "Curie-Weiss" naming postdates Weiss 1907 and Curie 1907.
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: -(J/N)*Σ σσ - h*Σ σ  (FM convention, h ≥ 0 → m* ≥ 0)
#   Observable:  Spin-1/2 (σ ∈ {±1})
#   Reference:   docs/src/conventions.md §Spin convention

"""
    CurieWeissIsing(; J::Real = 1.0, h::Real = 0.0) <: AbstractQAtlasModel

Classical Ising model on the complete graph (mean field) with
saddle-point Hamiltonian

    H = -(J/N) Σ_{i<j} σ_i σ_j  -  h Σ_i σ_i.

The 1/N normalisation makes the energy extensive and the thermodynamic
limit well-defined.  At `h = 0` and `J > 0` the model has zero-field
critical temperature `T_c = J` and mean-field critical exponents
`(β, γ, δ, ν) = (1/2, 1, 3, 1/2)` exposed by `Universality(:MeanField)`.

For `h ≠ 0` there is no sharp transition; the equilibrium magnetisation
`m*(β, J, h)` is the unique stable root of the self-consistency
equation `m = tanh(β(Jm + h))` (same sign as `h`).

Quantities registered:

| Quantity                              | BC         | Notes                                       |
| ------------------------------------- | ---------- | ------------------------------------------- |
| [`CriticalTemperature`](@ref)         | `Infinite` | `J` (h = 0 reference)                       |
| [`SpontaneousMagnetization`](@ref)    | `Infinite` | `h → 0⁺` limit (independent of model.h)     |
| [`CriticalExponents`](@ref)           | `Infinite` | MeanField universality                      |
| [`Energy`](@ref)`{:per_site}`         | `Infinite` | `-Jm²/2 - hm`                               |
| [`FreeEnergy`](@ref)                  | `Infinite` | `Jm²/2 - β⁻¹ log[2 cosh(β(Jm+h))]`         |
| [`ThermalEntropy`](@ref)              | `Infinite` | `log[2cosh(β(Jm+h))] - β(Jm+h)m`           |
| [`SpecificHeat`](@ref)                | `Infinite` | `β²(Jm+h)²(1-m²)/(1-βJ(1-m²))`             |
| [`SusceptibilityZZ`](@ref)            | `Infinite` | `β(1-m²)/(1-βJ(1-m²))`                     |

# References

- L. D. Landau, E. M. Lifshitz, *Statistical Physics* §149.
"""
struct CurieWeissIsing <: AbstractQAtlasModel
    J::Float64
    h::Float64
end
CurieWeissIsing(; J::Real=1.0, h::Real=0.0) = CurieWeissIsing(Float64(J), Float64(h))

# ═══════════════════════════════════════════════════════════════════════════════
# Critical temperature
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::CurieWeissIsing, ::CriticalTemperature, ::Infinite; J=m.J) -> Float64

Mean-field critical temperature `T_c = J` (with k_B = 1) for the
zero-field Curie-Weiss Ising model.  For `J ≤ 0` no ferromagnetic
order at any positive temperature; returns `0`.

At `h ≠ 0` there is no sharp transition — this dispatch still returns
`J` as the natural reference scale (the temperature at which `χ(h=0)`
diverges), independent of the model's `h` field.

# References
- L. D. Landau, E. M. Lifshitz, *Statistical Physics* §149.
"""
function fetch(
    m::CurieWeissIsing, ::CriticalTemperature, ::Infinite; J::Real=m.J, kwargs...
)
    return J > 0 ? Float64(J) : 0.0
end

# ═══════════════════════════════════════════════════════════════════════════════
# Self-consistency-equation (SCE) solver
# ═══════════════════════════════════════════════════════════════════════════════

# Bisection solver for the self-consistency equation
#
#     m = tanh(βJ m + βh)
#
# At βh = 0:  trivial root m = 0 if βJ ≤ 1; otherwise the unique
# nontrivial positive spontaneous root is found by bisection on
# g(m) = m - tanh(βJ m) over [m_min, 1 - eps) where m_min is the unique
# minimum of g (excludes the trivial root).
#
# At βh ≠ 0: a unique stable root has the same sign as h; bisection
# brackets [0, 1 - eps] for h > 0 (resp. [-1 + eps, 0] for h < 0)
# isolate it cleanly (the opposite-sign metastable root and the
# unstable root near 0 lie outside the bracket).
#
# Why not Picard / fixed-point iteration: the map m ↦ tanh(βJ m + βh)
# has derivative βJ sech²(...) which approaches 1 from below at the
# stable root as βJ → 1⁺, so the Picard rate degenerates near
# criticality and needs millions of iterations.  Bisection halves the
# bracket every step regardless.
function _curie_weiss_solve_m(
    beta_J::Real, beta_h::Real=0.0; tol::Real=1e-14, maxiter::Int=200
)
    if iszero(beta_h)
        # ── zero-field branch ────────────────────────────────────────
        if beta_J <= 1
            return 0.0
        end
        m_min = atanh(sqrt(1 - 1 / beta_J)) / beta_J
        lo = nextfloat(m_min)
        hi = 1.0 - 1e-15
        g_lo = lo - tanh(beta_J * lo)
        g_hi = hi - tanh(beta_J * hi)
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
    # ── h ≠ 0 branch: unique stable root has sign(h) ─────────────────
    g(m) = m - tanh(beta_J * m + beta_h)
    if beta_h > 0
        lo, hi = 0.0, 1.0 - 1e-15
    else
        lo, hi = -(1.0 - 1e-15), 0.0
    end
    g_lo, g_hi = g(lo), g(hi)
    if !(g_lo * g_hi < 0)
        # Bracket failed: happens in deep saturation (β|h| or βJ huge) where
        # the true root sits within an ULP of ±1.  Return the saturated
        # endpoint on the same side as h.
        return beta_h > 0 ? hi : lo
    end
    for _ in 1:maxiter
        mid = 0.5 * (lo + hi)
        g_mid = g(mid)
        if g_mid == 0 || (hi - lo) <= tol * max(1.0, abs(hi))
            return mid
        end
        if g_mid * g_lo < 0
            hi, g_hi = mid, g_mid
        else
            lo, g_lo = mid, g_mid
        end
    end
    return 0.5 * (lo + hi)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Spontaneous magnetisation — always the h → 0⁺ limit
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::CurieWeissIsing, ::SpontaneousMagnetization, ::Infinite;
          beta::Real, J=m.J) -> Float64

Spontaneous magnetisation `m*(β) = lim_{h → 0⁺} m(β, J, h)` of the
Curie-Weiss Ising model: the ℤ₂-positive nontrivial root of
`m = tanh(βJm)`.

By definition this dispatch always uses the **zero-field** SCE branch
(the model's `h` field is ignored), since "spontaneous" magnetisation
is the symmetry-broken value that survives in the field-free limit.
For `T ≥ T_c = J` only the trivial `m = 0` root exists and `0.0` is
returned; otherwise bisection on `g(m) = m - tanh(βJm)` over
`[m_min, 1)` converges in ~52 steps (no rate degeneracy at criticality).

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
    J > 0 || return 0.0
    return _curie_weiss_solve_m(beta * J, 0.0)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Critical exponents — delegate to mean-field universality
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::CurieWeissIsing, ::CriticalExponents, ::Infinite; kwargs...) -> NamedTuple

Mean-field (Landau-Ginzburg) critical exponents
`α = 0,  β = 1/2,  γ = 1,  δ = 3,  ν = 1/2,  η = 0` delegated to
[`MeanField`](@ref).  Independent of `h` (universality is a zero-field
concept).

# References
- L. D. Landau, *Phys. Z. Sowjet.* **11**, 26 (1937).
- H. E. Stanley, *Introduction to Phase Transitions and Critical Phenomena* (1971).
"""
function fetch(::CurieWeissIsing, ::CriticalExponents, ::Infinite; kwargs...)
    return QAtlas.fetch(QAtlas.MeanField(), CriticalExponents())
end

# ═══════════════════════════════════════════════════════════════════════════════
# FreeEnergy — saddle-point Helmholtz free energy per site
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::CurieWeissIsing, ::FreeEnergy, ::Infinite;
          beta::Real, J=m.J, h=m.h) -> Float64

Saddle-point free energy per site

    f(β, J, h) = J m*²/2 - β⁻¹ log[2 cosh(β(J m* + h))],

with `m*(β, J, h)` the unique stable SCE root.  Limits:

- `J ≤ 0`: non-interacting single spin in field, `f = -β⁻¹ log(2 cosh(βh))`.
- `J > 0`, `h = 0`, `T > T_c`: paramagnet, `f = -log(2)/β`.
- `J > 0`, `T → 0`: saturated, `f → -J/2 - |h|`.

# References
- L. D. Landau, E. M. Lifshitz, *Statistical Physics* §149.
"""
function fetch(
    m::CurieWeissIsing,
    ::FreeEnergy,
    ::Infinite;
    beta::Real,
    J::Real=m.J,
    h::Real=m.h,
    kwargs...,
)
    beta > 0 || throw(
        DomainError(beta, "CurieWeissIsing FreeEnergy requires β > 0; got β = $beta.")
    )
    if J <= 0
        return -log(2 * cosh(beta * h)) / beta
    end
    mag = _curie_weiss_solve_m(beta * J, beta * h)
    arg = beta * (J * mag + h)
    return J * mag^2 / 2 - log(2 * cosh(arg)) / beta
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy{:per_site} — mean-field internal energy
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::CurieWeissIsing, ::Energy{:per_site}, ::Infinite;
          beta::Real, J=m.J, h=m.h) -> Float64

Internal energy per site

    u(β, J, h) = -J m*² / 2 - h m*.

At `J ≤ 0`: single-spin, `u = -h tanh(βh)`.
At `J > 0`, `h = 0`, `T > T_c`: `u = 0` (paramagnet).
At `J > 0`, `T → 0`, `h ≥ 0`: `u → -J/2 - h` (saturated).

# References
- L. D. Landau, E. M. Lifshitz, *Statistical Physics* §149.
"""
function fetch(
    m::CurieWeissIsing,
    ::Energy{:per_site},
    ::Infinite;
    beta::Real,
    J::Real=m.J,
    h::Real=m.h,
    kwargs...,
)
    beta > 0 || throw(
        DomainError(
            beta, "CurieWeissIsing Energy{:per_site} requires β > 0; got β = $beta."
        ),
    )
    if J <= 0
        return -h * tanh(beta * h)
    end
    mag = _curie_weiss_solve_m(beta * J, beta * h)
    return -J * mag^2 / 2 - h * mag
end

# ═══════════════════════════════════════════════════════════════════════════════
# ThermalEntropy — per-site Gibbs entropy
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::CurieWeissIsing, ::ThermalEntropy, ::Infinite;
          beta::Real, J=m.J, h=m.h) -> Float64

Per-site entropy via the Gibbs identity `s = β(u - f)`:

    s(β, J, h) = log[2 cosh(β(J m* + h))]  -  β(J m* + h) m*.

Bounded between `0` (T → 0, saturated) and `log 2` (T → ∞, h finite).

# References
- L. D. Landau, E. M. Lifshitz, *Statistical Physics* §149.
"""
function fetch(
    m::CurieWeissIsing,
    ::ThermalEntropy,
    ::Infinite;
    beta::Real,
    J::Real=m.J,
    h::Real=m.h,
    kwargs...,
)
    beta > 0 || throw(
        DomainError(beta, "CurieWeissIsing ThermalEntropy requires β > 0; got β = $beta."),
    )
    if J <= 0
        return log(2 * cosh(beta * h)) - beta * h * tanh(beta * h)
    end
    mag = _curie_weiss_solve_m(beta * J, beta * h)
    arg = beta * (J * mag + h)
    return log(2 * cosh(arg)) - arg * mag
end

# ═══════════════════════════════════════════════════════════════════════════════
# SpecificHeat — c_v = (∂u/∂T) at fixed (J, h)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::CurieWeissIsing, ::SpecificHeat, ::Infinite;
          beta::Real, J=m.J, h=m.h) -> Float64

Specific heat per site

    c_v(β, J, h) = β² (J m* + h)² (1 - m*²) / (1 - β J (1 - m*²)).

Derived from `u(β) = -J m*²/2 - h m*` and `dm*/dβ` via implicit
differentiation of the SCE.  At `J > 0`, `h = 0`: zero for `T > T_c`,
jumps to `3/2` at `T_c⁻` (mean-field Landau jump).  At `J ≤ 0`:
non-interacting, `c_v = (βh)² sech²(βh)`.

# References
- L. D. Landau, E. M. Lifshitz, *Statistical Physics* §149.
"""
function fetch(
    m::CurieWeissIsing,
    ::SpecificHeat,
    ::Infinite;
    beta::Real,
    J::Real=m.J,
    h::Real=m.h,
    kwargs...,
)
    beta > 0 || throw(
        DomainError(beta, "CurieWeissIsing SpecificHeat requires β > 0; got β = $beta.")
    )
    if J <= 0
        return (beta * h * sech(beta * h))^2
    end
    mag = _curie_weiss_solve_m(beta * J, beta * h)
    if iszero(mag) && iszero(h)
        return 0.0
    end
    one_m_m2 = 1 - mag^2
    arg = beta * (J * mag + h)
    denom = 1 - beta * J * one_m_m2
    return arg^2 * one_m_m2 / denom
end

# ═══════════════════════════════════════════════════════════════════════════════
# SusceptibilityZZ — longitudinal isothermal susceptibility ∂m/∂h
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::CurieWeissIsing, ::SusceptibilityZZ, ::Infinite;
          beta::Real, J=m.J, h=m.h) -> Float64

Longitudinal isothermal susceptibility per site

    χ(β, J, h) = β (1 - m*²) / (1 - β J (1 - m*²)).

Derived by implicit differentiation of `m* = tanh(β(Jm* + h))` w.r.t. `h`.
Reduces to the Curie-Weiss law `χ = β / (1 - βJ)` above `T_c` (where
`m* = 0`).  Diverges at `T_c⁻` from the m* = 0 side.  At `J ≤ 0`:
non-interacting, `χ = β sech²(βh)`.

# References
- L. D. Landau, E. M. Lifshitz, *Statistical Physics* §149.
"""
function fetch(
    m::CurieWeissIsing,
    ::SusceptibilityZZ,
    ::Infinite;
    beta::Real,
    J::Real=m.J,
    h::Real=m.h,
    kwargs...,
)
    beta > 0 || throw(
        DomainError(
            beta, "CurieWeissIsing SusceptibilityZZ requires β > 0; got β = $beta."
        ),
    )
    if J <= 0
        return beta * sech(beta * h)^2
    end
    mag = _curie_weiss_solve_m(beta * J, beta * h)
    one_m_m2 = 1 - mag^2
    denom = 1 - beta * J * one_m_m2
    return beta * one_m_m2 / denom
end
