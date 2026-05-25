# ─────────────────────────────────────────────────────────────────────────────
# IsingChain1D — classical Ising chain (nearest-neighbour, 1-D).
#
# Hamiltonian:
#
#     H = -J Σ_i σ_i σ_{i+1} - h Σ_i σ_i,   σ_i ∈ {-1, +1}.
#
# Solved exactly by the 2×2 transfer matrix
#
#     T = [ e^{β(J + h)}   e^{-β J}      ;
#           e^{-β J}        e^{β(J - h)} ],
#
# with eigenvalues
#
#     λ_± = e^{β J} cosh(β h) ± √( e^{2 β J} sinh²(β h) + e^{-2 β J} ).
#
# The free energy per site in the thermodynamic limit is
#
#     f(β, h) = -β⁻¹ log λ_+,
#
# and the second-eigenvalue gap fixes the spin-spin correlation length:
#
#     ξ(β, h) = 1 / log(λ_+ / λ_-).
#
# At h = 0 these reduce to the celebrated Ising-1925 forms
#
#     f(β, 0) = -β⁻¹ log[2 cosh(β J)],
#     ξ(β, 0) = 1 / log(coth(β J)).
#
# There is no finite-temperature phase transition: the susceptibility,
# correlation length and free-energy derivatives are smooth in T for
# all T > 0; only at T = 0 do the two transfer-matrix eigenvalues
# become degenerate.
#
# Reference: E. Ising, Z. Phys. 31, 253 (1925).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: see file-header description above
#   Observable:  per src/core/quantities.jl (matches the dispatch tag)
#   Reference:   docs/src/conventions.md (project-wide convention policy)
#   STATUS:      backfilled by PR (audit gate); per-field domain content
#                left to a follow-up - see issue tracker for the model-specific
#                Hamiltonian sign / observable normalisation.

"""
    IsingChain1D(; J::Real = 1.0, h::Real = 0.0) <: AbstractQAtlasModel

Classical 1-D Ising chain with nearest-neighbour exchange `J` and
uniform longitudinal field `h`:

    H = -J Σ_i σ_i σ_{i+1} - h Σ_i σ_i,   σ_i ∈ {-1, +1}.

Solved exactly by the 2×2 transfer matrix (Ising 1925) with eigenvalues

    λ_± = e^{β J} cosh(β h) ± √( e^{2 β J} sinh²(β h) + e^{-2 β J} ).

There is no finite-temperature phase transition.

Quantities registered:

| Quantity                       | BC         | Method   |
| ------------------------------ | ---------- | -------- |
| [`CriticalTemperature`](@ref)  | `Infinite` | analytic |
| [`FreeEnergy`](@ref)           | `Infinite` | analytic |
| [`CorrelationLength`](@ref)    | `Infinite` | analytic |

# References

- E. Ising, *Z. Phys.* **31**, 253 (1925).
"""
struct IsingChain1D <: AbstractQAtlasModel
    J::Float64
    h::Float64
end
IsingChain1D(; J::Real=1.0, h::Real=0.0) = IsingChain1D(Float64(J), Float64(h))

# ═══════════════════════════════════════════════════════════════════════════════
# Transfer-matrix eigenvalues
# ═══════════════════════════════════════════════════════════════════════════════

# λ_± of the 1-D Ising transfer matrix.  Both are strictly positive
# for any (β > 0, J real, h real).
@inline function _ising_chain_1d_lambdas(β::Real, J::Real, h::Real)
    eJ = exp(β * J)
    em = exp(-β * J)
    c = eJ * cosh(β * h)
    r = sqrt(eJ^2 * sinh(β * h)^2 + em^2)
    return (c + r, c - r)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Critical temperature
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::IsingChain1D, ::CriticalTemperature, ::Infinite; kwargs...) -> Float64

Critical temperature of the 1-D Ising chain.  Ising (1925) proved
that no finite-temperature phase transition occurs in 1-D; the only
singular point is `T = 0`.  Returns `0.0`.

# References

- E. Ising, *Z. Phys.* **31**, 253 (1925).
"""
function fetch(::IsingChain1D, ::CriticalTemperature, ::Infinite; kwargs...)
    return 0.0
end

# ═══════════════════════════════════════════════════════════════════════════════
# Free energy per site
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::IsingChain1D, ::FreeEnergy, ::Infinite;
          beta::Real, J=m.J, h=m.h) -> Float64

Helmholtz free energy per site `f(β, h) = -β⁻¹ log λ_+(β, J, h)` of the
1-D Ising chain in the thermodynamic limit, with `λ_+` the larger
transfer-matrix eigenvalue.  At `h = 0` this reduces to the textbook
form `f = -β⁻¹ log[2 cosh(β J)]`.

# References

- E. Ising, *Z. Phys.* **31**, 253 (1925).
"""
function fetch(
    m::IsingChain1D,
    ::FreeEnergy,
    ::Infinite;
    beta::Real,
    J::Real=m.J,
    h::Real=m.h,
    kwargs...,
)
    beta > 0 ||
        throw(DomainError(beta, "IsingChain1D FreeEnergy requires β > 0; got β = $beta."))
    λp, _ = _ising_chain_1d_lambdas(beta, J, h)
    return -log(λp) / beta
end

# ═══════════════════════════════════════════════════════════════════════════════
# Correlation length
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::IsingChain1D, ::CorrelationLength, ::Infinite;
          beta::Real, J=m.J, h=m.h) -> Float64

Spin-spin correlation length `ξ(β, h) = 1 / log(λ_+ / λ_-)` of the 1-D
Ising chain, with `λ_±` the two transfer-matrix eigenvalues.  At
`h = 0` this reduces to the Ising 1925 form

    ξ(β, 0) = 1 / log(coth(β J)).

The correlation length is finite for every `T > 0` (no
finite-temperature phase transition) and diverges only in the
zero-temperature ferromagnetic limit `T → 0⁺` at `J > 0`.

# References

- E. Ising, *Z. Phys.* **31**, 253 (1925).
"""
function fetch(
    m::IsingChain1D,
    ::CorrelationLength,
    ::Infinite;
    beta::Real,
    J::Real=m.J,
    h::Real=m.h,
    kwargs...,
)
    beta > 0 || throw(
        DomainError(beta, "IsingChain1D CorrelationLength requires β > 0; got β = $beta."),
    )
    λp, λm = _ising_chain_1d_lambdas(beta, J, h)
    λm > 0 || return Inf       # degenerate eigenvalues ⇒ ξ → ∞ (e.g. T = 0 FM).
    ratio = λp / λm
    ratio > 1 || return Inf
    return 1 / log(ratio)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Thermodynamics at h = 0 — closed-form 1D Ising textbook expressions
#
# All five quantities below are valid at h = 0 only (textbook scope, Ising 1925).
# The h ≠ 0 generalisation requires symbolic differentiation of `λ_+(β, J, h)`
# and is deferred. A `DomainError` is raised on `h ≠ 0` for clarity.
# ═══════════════════════════════════════════════════════════════════════════════

@inline function _isingchain1d_require_zero_field(h, qname)
    iszero(h) || throw(
        DomainError(
            h,
            "IsingChain1D $(qname) currently requires h = 0; got h = $h. " *
            "The h ≠ 0 closed forms require differentiating the general λ_+ formula " *
            "and have not been implemented (closed forms exist; this is scope, not feasibility).",
        ),
    )
end

"""
    fetch(::IsingChain1D, ::Energy{:per_site}, ::Infinite;
          beta::Real, J=m.J, h=m.h) -> Float64

Internal energy per site at zero field (Ising 1925):

    u(β, h=0) = ⟨H⟩/N = -J tanh(β J).

# References
- E. Ising, *Z. Phys.* **31**, 253 (1925).
"""
function fetch(
    m::IsingChain1D,
    ::Energy{:per_site},
    ::Infinite;
    beta::Real,
    J::Real=m.J,
    h::Real=m.h,
    kwargs...,
)
    beta > 0 || throw(
        DomainError(beta, "IsingChain1D Energy{:per_site} requires β > 0; got β = $beta."),
    )
    _isingchain1d_require_zero_field(h, "Energy{:per_site}")
    return -J * tanh(beta * J)
end

"""
    fetch(::IsingChain1D, ::SpecificHeat, ::Infinite;
          beta::Real, J=m.J, h=m.h) -> Float64

Specific heat per site at zero field (Ising 1925):

    c_v(β, h=0) = (β J)² sech²(β J).

# References
- E. Ising, *Z. Phys.* **31**, 253 (1925).
"""
function fetch(
    m::IsingChain1D,
    ::SpecificHeat,
    ::Infinite;
    beta::Real,
    J::Real=m.J,
    h::Real=m.h,
    kwargs...,
)
    beta > 0 ||
        throw(DomainError(beta, "IsingChain1D SpecificHeat requires β > 0; got β = $beta."))
    _isingchain1d_require_zero_field(h, "SpecificHeat")
    bJ = beta * J
    return (bJ * sech(bJ))^2
end

"""
    fetch(::IsingChain1D, ::ThermalEntropy, ::Infinite;
          beta::Real, J=m.J, h=m.h) -> Float64

Entropy per site at zero field (Ising 1925), via the Gibbs identity
`s = β(u − f)` applied to the 1D Ising closed forms:

    s(β, h=0) = log(2 cosh(β J)) − β J tanh(β J).

Bounded between 0 (T → 0) and log 2 (T → ∞).

# References
- E. Ising, *Z. Phys.* **31**, 253 (1925).
"""
function fetch(
    m::IsingChain1D,
    ::ThermalEntropy,
    ::Infinite;
    beta::Real,
    J::Real=m.J,
    h::Real=m.h,
    kwargs...,
)
    beta > 0 || throw(
        DomainError(beta, "IsingChain1D ThermalEntropy requires β > 0; got β = $beta.")
    )
    _isingchain1d_require_zero_field(h, "ThermalEntropy")
    bJ = beta * J
    return log(2 * cosh(bJ)) - bJ * tanh(bJ)
end

"""
    fetch(::IsingChain1D, ::SusceptibilityZZ, ::Infinite;
          beta::Real, J=m.J, h=m.h) -> Float64

Zero-field longitudinal susceptibility per site (Brush 1967):

    χ(β, h=0) = β e^{2 β J}.

# References
- E. Ising, *Z. Phys.* **31**, 253 (1925).
- S. G. Brush, "History of the Lenz-Ising model", *Rev. Mod. Phys.* **39**,
  883 (1967), Eq. (4.18).
"""
function fetch(
    m::IsingChain1D,
    ::SusceptibilityZZ,
    ::Infinite;
    beta::Real,
    J::Real=m.J,
    h::Real=m.h,
    kwargs...,
)
    beta > 0 || throw(
        DomainError(beta, "IsingChain1D SusceptibilityZZ requires β > 0; got β = $beta."),
    )
    _isingchain1d_require_zero_field(h, "SusceptibilityZZ")
    return beta * exp(2 * beta * J)
end

"""
    fetch(::IsingChain1D, ::SpontaneousMagnetization, ::Infinite;
          beta::Real, J=m.J, h=m.h) -> Float64

Spontaneous magnetization per site of the 1D Ising chain — `0` for all
`T > 0` and any `J ≠ 0` (Ising 1925: no spontaneous symmetry breaking
in 1D at finite temperature).

# References
- E. Ising, *Z. Phys.* **31**, 253 (1925).
"""
function fetch(
    m::IsingChain1D,
    ::SpontaneousMagnetization,
    ::Infinite;
    beta::Real,
    J::Real=m.J,
    h::Real=m.h,
    kwargs...,
)
    beta > 0 || throw(
        DomainError(
            beta, "IsingChain1D SpontaneousMagnetization requires β > 0; got β = $beta."
        ),
    )
    return 0.0
end
