# ─────────────────────────────────────────────────────────────────────────────
# Transverse Field Ising Model — exact solutions
#
# Hamiltonian:
#   H = -J Σᵢ σᶻᵢσᶻᵢ₊₁  -  h Σᵢ σˣᵢ
#
# Solved exactly via Jordan-Wigner + Bogoliubov-de Gennes (BdG) transformation.
# The quadratic fermion Hamiltonian has quasiparticle energies Λₙ > 0, giving:
#
#   ⟨H⟩(β) = -Σₙ (Λₙ/2) tanh(β Λₙ / 2)
#
# The canonical API uses the concrete `TFIM` struct and concrete `Quantity`
# types from `src/core/quantities.jl`.  Legacy symbol-dispatch
# (`fetch(:TFIM, :energy, OBC(); …)`) routes through
# `src/deprecate/legacy_tfim.jl`.
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Pauli σ (this file)
#   Observable:  Spin S = σ/2  (QAtlas-wide spin convention; see docs/src/conventions.md)

using LinearAlgebra: eigvals, Symmetric
using QuadGK: quadgk

"""
    TFIM(; J = 1.0, h = 1.0) <: AbstractQAtlasModel

The 1D transverse field Ising model with Hamiltonian

    H = -J Σ_i σᶻ_i σᶻ_{i+1} - h Σ_i σˣ_i

`J > 0` is ferromagnetic, `h` is the transverse field.  The critical
point sits at `h = J`.

Currently registered fetches:

| Quantity                   | BC                 | Coverage                                                              |
| -------------------------- | ------------------ | --------------------------------------------------------------------- |
| [`Energy`](@ref)           | `OBC` / `Infinite` | Exact energy computed via BdG transformation                          |
| [`SpecificHeat`](@ref)     | `Infinite`         | Specific heat at finite temperature                                   |
| [`FreeEnergy`](@ref)       | `Infinite`         | Free energy density at finite temperature                             |
| [`ThermalEntropy`](@ref)   | `Infinite`         | Thermal entropy density at finite temperature                         |
| [`UniversalityClass`](@ref) | `Infinite`         | `:Ising` universality class at the critical point `h = J` (flows to `:IsingSDRG` under strong disorder) |
"""
struct TFIM <: AbstractQAtlasModel
    J::Float64
    h::Float64
end
TFIM(; J::Real=1.0, h::Real=1.0) = TFIM(Float64(J), Float64(h))

# ═══════════════════════════════════════════════════════════════════════════════
# Internal: BdG quasiparticle spectrum (OBC, finite N)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _tfim_bdg_spectrum(N, J, h) -> Vector{Float64}

Return the N positive BdG quasiparticle energies Λₙ > 0 for the OBC TFIM
with N sites, Ising coupling J, and transverse field h.

The 2N×2N BdG matrix is:
    H_BdG = [[A, B]; [-B, -A]]
where A (tridiagonal, symmetric) encodes hopping + onsite energy,
and B (antisymmetric) encodes the pairing terms from JW transformation.

    A_{ii}   = 2h
    A_{i,i±1} = -J
    B_{i,i+1} = +J,  B_{i+1,i} = -J
"""
function _tfim_bdg_spectrum(N::Int, J::Float64, h::Float64)::Vector{Float64}
    A = zeros(N, N)
    for i in 1:N
        A[i, i] = 2h
    end
    for i in 1:(N - 1)
        A[i, i + 1] = -J
        A[i + 1, i] = -J
    end

    B = zeros(N, N)
    for i in 1:(N - 1)
        B[i, i + 1] = J
        B[i + 1, i] = -J
    end

    H_bdg = [A B; -B -A]
    vals = eigvals(Symmetric(H_bdg))
    return sort!(filter(v -> v > 1e-10, vals))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy granularity convention (see src/core/quantities.jl)
# ═══════════════════════════════════════════════════════════════════════════════

native_energy_granularity(::TFIM, ::OBC) = :total
native_energy_granularity(::TFIM, ::Infinite) = :per_site

# ═══════════════════════════════════════════════════════════════════════════════
# Energy: OBC finite-N
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::TFIM, ::Energy{:total}, bc::OBC; beta, betas) -> Float64 or Vector{Float64}

Total energy ⟨H⟩(β) for the OBC TFIM with N sites.  Native granularity
for finite-N TFIM (per-site is provided by the generic conversion
fallback in `src/core/quantities.jl`).

- `N` is read from `bc.N` (`OBC(N)` / `OBC(; N)`) or from `kwargs[:N]`
  as a legacy fallback.
- `beta::Float64`: return scalar ⟨H⟩(β)
- `betas::AbstractVector{Float64}`: return vector, reusing spectrum (O(N³) once)
- no keyword: return ground-state energy E₀ = -Σₙ Λₙ/2  (β → ∞)

Uses the exact BdG formula:  ⟨H⟩ = -Σₙ (Λₙ/2) tanh(β Λₙ / 2)
"""
function fetch(
    model::TFIM,
    ::Energy{:total},
    bc::OBC;
    beta::Union{Real,Nothing}=nothing,
    betas::Union{AbstractVector{<:Real},Nothing}=nothing,
    kwargs...,
)
    N = _bc_size(bc, kwargs)
    Λ = _tfim_bdg_spectrum(N, model.J, model.h)
    if betas !== nothing
        return [-sum(λ -> (λ / 2) * tanh(β * λ / 2), Λ) for β in betas]
    elseif beta !== nothing
        return -sum(λ -> (λ / 2) * tanh(beta * λ / 2), Λ)
    else
        # Ground state: β → ∞, tanh(β Λ/2) → 1
        return -sum(Λ) / 2
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy: thermodynamic limit (PBC / Infinite)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::TFIM, ::Energy{:per_site}, ::Infinite; beta, betas) -> Float64 or Vector{Float64}

Energy *per site* ⟨H⟩/N in the thermodynamic limit (PBC, N → ∞).
Native granularity at `Infinite()` (total energy diverges and has no
defined value here).

    ε(β) = -(1/π) ∫₀^π dk  Λ(k)/2 · tanh(β Λ(k) / 2)

where the PBC dispersion is  Λ(k) = 2√(J² + h² - 2Jh cos k).

- `beta::Float64`: return scalar ε(β)
- `betas::AbstractVector{Float64}`: return vector
- no keyword: return ground-state energy per site (β → ∞)

Uses adaptive Gauss-Kronrod quadrature (QuadGK).
"""
function fetch(
    model::TFIM,
    ::Energy{:per_site},
    ::Infinite;
    beta::Union{Real,Nothing}=nothing,
    betas::Union{AbstractVector{<:Real},Nothing}=nothing,
    kwargs...,
)
    J = model.J
    h = model.h
    _energy_at_beta =
        β -> begin
            result, _ = quadgk(
                k -> begin
                    Λk = 2sqrt(J^2 + h^2 - 2J * h * cos(k))
                    (Λk / 2) * tanh(β * Λk / 2)
                end, 0.0, π; rtol=1e-10
            )
            -(1 / π) * result
        end
    if betas !== nothing
        return [_energy_at_beta(β) for β in betas]
    elseif beta !== nothing
        return _energy_at_beta(beta)
    else
        # Ground state: β → ∞
        return _energy_at_beta(1e6)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Mass gap (lowest quasi-particle energy)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::TFIM, ::MassGap, ::Infinite) -> Float64

Mass gap of the infinite-chain TFIM: the lowest single-quasiparticle
excitation energy

    Δ = min_k Λ(k),     Λ(k) = 2 √( J² + h² − 2 J h cos k ).

Closed form:

    Δ = 2 |h − J|.

Canonical values:

- ordered   (h < J): `Δ = 2(J − h)`
- disordered (h > J): `Δ = 2(h − J)`
- critical  (h = J): `Δ = 0` (Ising CFT, Δ ~ π v_F / N on finite chains)
"""
function fetch(model::TFIM, ::MassGap, ::Infinite; kwargs...)
    return 2 * abs(model.h - model.J)
end

"""
    fetch(model::TFIM, ::MassGap, bc::OBC) -> Float64

Single-quasiparticle gap of the N-site OBC TFIM read off the BdG
spectrum as `Λ_min`, the smallest positive eigenvalue of the 2N×2N
Bogoliubov-de Gennes Hamiltonian.

This is the one-particle excitation energy.  Away from the critical
point (`|h − J| > O(1/N)`) it converges to `2|h − J|` exponentially in
N.  At the critical point `h = J` the OBC gap scales as
`Δ(N) ~ π J / N` (Ising CFT).

Size is taken from `bc.N` (or `kwargs[:N]` as a legacy fallback).
"""
function fetch(model::TFIM, ::MassGap, bc::OBC; kwargs...)
    N = _bc_size(bc, kwargs)
    Λ = _tfim_bdg_spectrum(N, model.J, model.h)
    return Λ[1]
end

# ═══════════════════════════════════════════════════════════════════════════════
# Central charge (critical point h = J)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::TFIM, ::CentralCharge, ::Infinite) -> Float64

Central charge of the TFIM:

- `c = 1/2` at the critical point `h = J` (Ising CFT)
- `c = 0`   in either gapped phase (`h ≠ J`) — no low-energy CFT description

Criticality is detected by `|h/J - 1| ≤ 1e-6`.

# Example

```jldoctest
julia> QAtlas.fetch(TFIM(), CentralCharge(), Infinite())
0.5
```
"""
function fetch(model::TFIM, ::CentralCharge, ::Infinite; kwargs...)
    return abs(model.h / model.J - 1.0) ≤ 1e-6 ? 0.5 : 0.0
end

# ═══════════════════════════════════════════════════════════════════════════════
# Lieb-Robinson velocity bound (status=:bound) — status-axis worked example (v0.24)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::TFIM, ::LiebRobinsonBound, ::Infinite) -> Float64

Lieb-Robinson velocity of the infinite-chain TFIM,

    v_LR = 2 min(|J|, |h|),

the slope of the causal cone bounding commutator spread.  For the
free-fermion TFIM this tight bound is saturated by the maximum group
velocity `max_k |dΛ/dk|` of the Bogoliubov dispersion
`Λ(k) = 2√(J² + h² − 2 J h cos k)`.  Registered with `status=:bound`.

(Lieb & Robinson 1972; Hastings & Koma 2006.)
"""
function fetch(model::TFIM, ::LiebRobinsonBound, ::Infinite; kwargs...)
    return 2 * min(abs(model.J), abs(model.h))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Critical exponents at the quantum critical point — delegate to 2D Ising universality
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::TFIM, ::CriticalExponents, ::Infinite; kwargs...) -> NamedTuple

Onsager 2D-Ising critical exponents at the TFIM quantum critical
point `h = J`, delegated to the existing `Universality(:Ising)`
infrastructure at `d = 2`:

    β = 1/8,  γ = 7/4,  δ = 15,  ν = 1,  α = 0,  η = 1/4.

The 1D TFIM is exactly equivalent to the 2D classical Ising model via
the quantum-classical mapping (Pfeuty 1970), so the universal critical
exponents are identical to Onsager's 1944 result.

# References

- L. Onsager, *Phys. Rev.* **65**, 117 (1944) — 2D classical Ising exact solution.
- P. Pfeuty, *Ann. Phys.* **57**, 79 (1970) — TFIM ↔ 2D Ising equivalence.
- S. Sachdev, *Quantum Phase Transitions* (2nd ed., Cambridge 2011) — TFIM as canonical QPT example.
"""
function fetch(::TFIM, ::CriticalExponents, ::Infinite; kwargs...)
    return QAtlas.fetch(QAtlas.Universality(:Ising), CriticalExponents(); d=2)
end

# ─────────────────────────────────────────────────────────────────────────────
# Lieb-Robinson velocity (#579 inequality framework Phase 1)
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(model::TFIM, ::LiebRobinsonVelocity, ::Infinite;
          J=m.J, h=m.h) -> Float64

Lieb-Robinson velocity of the transverse-field Ising chain. Via the
Jordan-Wigner mapping the TFIM is a free Bogoliubov-fermion system
with dispersion `Λ(k) = 2 sqrt(J^2 + h^2 - 2 J h cos k)`. The tight
Lieb-Robinson velocity is the maximum single-particle group velocity
saturating the bound: differentiating `Λ(k)` and locating the
interior stationary point at `cos k = min(|J|, |h|) / max(|J|, |h|)`
gives

    v_LR = max_k |dΛ/dk| = 2 min(|J|, |h|).

At criticality `h = J` this is `2J = 2h` (Calabrese-Cardy 2006).
The Hastings-Koma upper bound `2 max(|J|, |h|)` is loose; the value
returned here is the *tight* free-fermion saturated velocity that
governs e.g. the linear-growth slope of post-quench entanglement
(see PR #588 EntanglementGrowthSlope).

At `h = 0` (classical Ising) or `J = 0` (decoupled site spins) the
chain has no quantum dynamics and `v_LR = 0`.

Reference: Lieb-Robinson 1972; Hastings-Koma 2006 (general bound);
Calabrese-Cardy 2006 (free-fermion saturation in quench dynamics).
"""
function fetch(
    model::TFIM,
    ::LiebRobinsonVelocity,
    ::Infinite;
    J::Real=model.J,
    h::Real=model.h,
    kwargs...,
)
    return 2 * min(abs(J), abs(h))
end

# ─────────────────────────────────────────────────────────────────────────────
# EntanglementGrowthSlope wrapper at h = J critical (#580 / #579 cross)
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(model::TFIM, ::EntanglementGrowthSlope, ::Infinite;
          beta_eff::Real, kwargs...) -> Float64

Linear-growth slope of post-quench half-system entanglement entropy
for the TFIM at the Ising critical point `h = J`. Wires together two
universality-layer pieces

    c = 1/2          (Universality(:Ising) CentralCharge)
    v_LR = 2 |J|     (TFIM LiebRobinsonVelocity at h = J critical)

into the Calabrese-Cardy 2005 result

    dS_A/dt = π c v_LR / (3 beta_eff) = π J / (3 beta_eff).

For non-critical TFIM (`h ≠ J`, gapped) the CC linear-growth picture
does not apply and `DomainError` is thrown.

Reference: Calabrese-Cardy *J. Stat. Mech.* P04010 (2005);
combines universality-layer dispatches from PR #588 and the TFIM
LiebRobinsonVelocity from PR #586 / fix #592.
"""
function fetch(
    model::TFIM, ::EntanglementGrowthSlope, ::Infinite; beta_eff::Real, kwargs...
)
    isapprox(model.h, model.J; atol=1e-10) || throw(
        DomainError(
            (model.J, model.h),
            "TFIM EntanglementGrowthSlope at Infinite is defined only at the Ising " *
            "critical point h = J (gapless c = 1/2 CFT); off-critical TFIM is gapped " *
            "and CC linear-growth does not apply. Got (J, h) = (\$(model.J), \$(model.h)).",
        ),
    )
    return fetch(
        Universality(:Ising),
        EntanglementGrowthSlope(),
        Infinite();
        v=fetch(model, LiebRobinsonVelocity(), Infinite()),
        beta_eff=beta_eff,
        kwargs...,
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# NMR relaxation exponent (quantum critical point h = J)
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(model::TFIM, ::NMRRelaxationExponent, ::Infinite; kwargs...) -> Float64

NMR spin-lattice relaxation rate temperature scaling exponent `θ_{NMR} = -3/4`
at the quantum critical point `h = J`. For non-critical `h ≠ J`, returns `NaN` with a warning.
"""
function fetch(model::TFIM, ::NMRRelaxationExponent, ::Infinite; kwargs...)
    if abs(model.h - model.J) ≤ 1e-6
        return -0.75
    else
        @warn "TFIM NMRRelaxationExponent is only defined at the quantum critical point h = J." h =
            model.h J = model.J
        return NaN
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Conformal tower of states (quantum critical point h = J)
# ─────────────────────────────────────────────────────────────────────────────

raw"""
    fetch(model::TFIM, q::ConformalTower, bc::Union{PBC, OBC}; kwargs...) -> Vector{NamedTuple}

Conformal tower of states excitation energies of the TFIM chain at the quantum critical
point `h = J` (Ising CFT, M(4,3), c=1/2). Delegates to `Universality(:Ising)` at
boundary condition `bc` with the exact free-fermion Fermi velocity `v = 2J` and system
size `L` extracted from `bc`.

The three Ising primary operators and their scaling dimensions are:
- Identity:  (h, h̄) = (0,    0   ),  Δ = 0
- Spin σ:    (h, h̄) = (1/16, 1/16),  Δ = 1/8
- Energy ε:  (h, h̄) = (1/2,  1/2 ),  Δ = 1

Throws a `DomainError` if the model is off-critical (`|h - J| > 1e-6`).

# Arguments
- `bc::Union{PBC, OBC}`: boundary condition; system size `L` is read from `bc.N`.
- Keyword `N` accepted as a legacy fallback if `bc` does not carry a size.

# Returns
`Vector{NamedTuple{(:energy, :dimension, :degeneracy), Tuple{Float64, Float64, Int}}}` —
see `fetch(::Universality{:Ising}, ::ConformalTower, ...)` for full field documentation.

# References
- J. Cardy, *Nucl. Phys. B* **270**, 186 (1986). — operator content of 1+1D CFTs.
- H. W. J. Blöte, J. L. Cardy, M. P. Nightingale, *Phys. Rev. Lett.* **56**, 742 (1986).
- P. Pfeuty, *Ann. Phys.* **57**, 79 (1970). — TFIM Fermi velocity `v = 2J` at h = J.
"""
function fetch(model::TFIM, q::ConformalTower, bc::Union{PBC,OBC}; kwargs...)
    isapprox(model.h, model.J; atol=1e-6) || throw(
        DomainError(
            (model.J, model.h),
            "TFIM ConformalTower is defined only at the Ising critical point h = J. Got (J, h) = ($(model.J), $(model.h)).",
        ),
    )
    L = _bc_size(bc, kwargs)
    v = 2.0 * model.J
    return fetch(Universality(:Ising), q, bc; L=L, v=v, kwargs...)
end
