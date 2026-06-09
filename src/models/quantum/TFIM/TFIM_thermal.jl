# ─────────────────────────────────────────────────────────────────────────────
# Transverse Field Ising Model — exact finite-temperature thermodynamics
#
# Hamiltonian:
#   H = -J Σᵢ σᶻᵢσᶻᵢ₊₁  -  h Σᵢ σˣᵢ
#
# All observables in this file are derived from the free-fermion (BdG)
# diagonalisation of the JW-mapped quadratic Hamiltonian.  The single-particle
# dispersion in the thermodynamic limit is
#
#   Λ(k) = 2 √(J² + h² - 2 J h cos k),       k ∈ [0, π]
#
# The grand-canonical free-fermion partition function then gives
#
#   log Z / N = (1/π) ∫₀^π dk · log( 2 cosh(β Λ(k)/2) )
#
# from which all thermodynamic potentials follow as standard derivatives.
#
# For OBC finite N the same expressions hold with the integral replaced by
# a sum over the N positive BdG quasiparticle energies returned by
# `_tfim_bdg_spectrum(N, J, h)` in `TFIM.jl`.
#
# Quantities exposed via `fetch`:
#
#   :free_energy              f(β)        = -T log Z / N
#   :entropy                  s(β)        = β (ε - f)
#   :specific_heat            c_v(β)      = ∂ε/∂T
#   :transverse_magnetization m_x(β)      = ⟨σˣ_i⟩
#   :transverse_susceptibility χ_xx(β)    = β · Var(Σᵢ σˣᵢ) / N
#
# All five are implemented for both `OBC` (per-site, exact at finite N) and
# `Infinite` (per-site, exact in the thermodynamic limit).
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra: eigvals, Symmetric
using QuadGK: quadgk

# ═══════════════════════════════════════════════════════════════════════════════
# Internal: Bogoliubov building blocks
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _tfim_dispersion(k, J, h) -> Float64

Single-particle BdG quasiparticle energy at momentum `k` for the TFIM with
couplings `J` (Z-Z coupling) and `h` (transverse field):

    Λ(k) = 2 √(J² + h² - 2 J h cos k).
"""
@inline _tfim_dispersion(k::Real, J::Real, h::Real) =
    2 * sqrt(J^2 + h^2 - 2 * J * h * cos(k))

# Kernel `g(βλ/2)` style helpers — the integrands are written in terms of `λ` and `β` only.
# A small helper avoids overflow in `log(2 cosh(x))` for large `|x|`.
@inline function _logcosh2(x::Real)
    # log(2 cosh(x)) = |x| + log(1 + exp(-2|x|))
    a = abs(x)
    return a + log1p(exp(-2 * a))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Thermodynamic potentials — Infinite (per-site)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _tfim_thermo_infinite(quantity, J, h, β) -> Float64

Compute one of the per-site thermodynamic potentials of the infinite TFIM at
inverse temperature `β`.  `quantity` is a `Symbol` from
`(:free_energy, :entropy, :specific_heat, :transverse_magnetization,
 :transverse_susceptibility)`.

The integrals are evaluated via adaptive Gauss-Kronrod quadrature.
"""
function _tfim_thermo_infinite(quantity::Symbol, J::Real, h::Real, β::Real; kwargs...)
    integrand = if quantity === :free_energy
        # f = -(1/πβ) ∫ log(2 cosh(βΛ/2)) dk
        k -> begin
            Λk = _tfim_dispersion(k, J, h)
            _logcosh2(β * Λk / 2)
        end
    elseif quantity === :entropy
        # s = (1/π) ∫ [log(2 cosh(βΛ/2)) - (βΛ/2) tanh(βΛ/2)] dk
        k -> begin
            Λk = _tfim_dispersion(k, J, h)
            x = β * Λk / 2
            _logcosh2(x) - x * tanh(x)
        end
    elseif quantity === :specific_heat
        # c_v = (1/π) ∫ (βΛ/2)² sech²(βΛ/2) dk
        k -> begin
            Λk = _tfim_dispersion(k, J, h)
            x = β * Λk / 2
            x^2 * sech(x)^2
        end
    elseif quantity === :transverse_magnetization
        # m_x = (2/π) ∫ ((h - J cos k)/Λ) tanh(βΛ/2) dk
        k -> begin
            A = h - J * cos(k)
            Λk = _tfim_dispersion(k, J, h)
            (A / Λk) * tanh(β * Λk / 2)
        end
    elseif quantity === :transverse_susceptibility
        # χ_xx = (2/π) ∫ [ (1/Λ - 4A²/Λ³) tanh(βΛ/2) + (2β A²/Λ²) sech²(βΛ/2) ] dk
        k -> begin
            A = h - J * cos(k)
            Λk = _tfim_dispersion(k, J, h)
            (1 / Λk - 4 * A^2 / Λk^3) * tanh(β * Λk / 2) +
            (2 * β * A^2 / Λk^2) * sech(β * Λk / 2)^2
        end
    elseif quantity === :nmr_relaxation
        eta = Float64(get(kwargs, :eta, 0.1))
        eta > 0 || throw(DomainError(eta, "TFIM NMRSpinRelaxationRate requires η > 0; got η = $eta."))
        β > 0 || throw(DomainError(β, "TFIM NMRSpinRelaxationRate requires β > 0; got β = $β."))
        return _tfim_nmr_relaxation_infinite(J, h, β, eta)
    else
        error("Unknown thermal quantity: $quantity")
    end

    if quantity === :nmr_relaxation
        return integrand # already computed relaxation rate
    end

    val, _ = quadgk(integrand, 0.0, π; rtol=1e-10)

    if quantity === :free_energy
        return -val / (π * β)
    elseif quantity === :transverse_magnetization || quantity === :transverse_susceptibility
        return (2 / π) * val
    else  # entropy, specific_heat
        return val / π
    end
end

function _tfim_nmr_relaxation_infinite(J::Real, h::Real, β::Real, η::Real)
    integrand_outer = k1 -> begin
        λ1 = _tfim_dispersion(k1, J, h)
        f1 = λ1 > 0 ? exp(-β * λ1) / (1.0 + exp(-β * λ1)) : 0.5
        integrand_inner = k2 -> begin
            λ2 = _tfim_dispersion(k2, J, h)
            f2 = λ2 > 0 ? exp(-β * λ2) / (1.0 + exp(-β * λ2)) : 0.5
            lorentz = η / (π * ((λ1 - λ2)^2 + η^2))
            return f1 * (1.0 - f2) * lorentz
        end
        val_inner, _ = quadgk(integrand_inner, 0.0, π; rtol=1e-6)
        return val_inner
    end
    val_outer, _ = quadgk(integrand_outer, 0.0, π; rtol=1e-6)
    return val_outer / π^2
end

# ═══════════════════════════════════════════════════════════════════════════════
# Thermodynamic potentials — OBC finite N (per-site)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _tfim_thermo_obc(quantity, N, J, h, β) -> Float64

Per-site thermodynamic potential for the OBC finite-N TFIM, computed by summing
the contribution of each BdG quasiparticle mode.

The transverse magnetisation and its susceptibility require the full
single-particle Bogoliubov coefficients, not just the spectrum, so this routine
diagonalises the BdG matrix internally to obtain them.
"""
function _tfim_thermo_obc(quantity::Symbol, N::Int, J::Float64, h::Float64, β::Real)
    if quantity === :free_energy
        Λ = _tfim_bdg_spectrum(N, J, h)
        # f/N = -(1/Nβ) Σ log(2 cosh(βΛ/2))
        return -sum(λ -> _logcosh2(β * λ / 2), Λ) / (N * β)
    elseif quantity === :entropy
        Λ = _tfim_bdg_spectrum(N, J, h)
        return sum(Λ) do λ
            x = β * λ / 2
            return _logcosh2(x) - x * tanh(x)
        end / N
    elseif quantity === :specific_heat
        Λ = _tfim_bdg_spectrum(N, J, h)
        return sum(λ -> begin
            x = β * λ / 2
            x^2 * sech(x)^2
        end, Λ) / N
    elseif quantity === :transverse_magnetization || quantity === :transverse_susceptibility
        return _tfim_transverse_obc(quantity, N, J, h, β)
    elseif quantity === :nmr_relaxation
        error("NMRSpinRelaxationRate is not implemented for OBC in TFIM.")
    else
        error("Unknown thermal quantity: $quantity")
    end
end

"""
    _xx_uniform_susceptibility(N, J, h, β) -> Float64

Exact transverse susceptibility per site for the OBC TFIM,

    χ_xx(β) = (β/N) Var(Σᵢ σˣᵢ)
            = (β/N) Σᵢⱼ [ ⟨σˣᵢ σˣⱼ⟩_β − ⟨σˣᵢ⟩_β ⟨σˣⱼ⟩_β ]

Uses the Majorana covariance matrix `Σ[a,b] = ⟨γₐγᵦ⟩ − δₐᵦ`.
With `σˣᵢ = -i γ_{2i-1} γ_{2i}` the connected correlators follow from
Wick's theorem:

  Diagonal (i = j):   ⟨(σˣᵢ)²⟩_c = 1 − Σ[2i-1, 2i]²
  Off-diagonal (i ≠ j): ⟨σˣᵢ σˣⱼ⟩_c = −Σ[2i-1,2j-1]·Σ[2i,2j]
                                        + Σ[2i-1,2j]·Σ[2i,2j-1]

No numerical differentiation; no Pfaffian library calls.
"""
function _xx_uniform_susceptibility(N::Int, J::Float64, h::Float64, β::Real)
    hmat = _majorana_ham(N, J, h)
    Σ = _majorana_thermal_covariance(hmat, β)
    # ⟨σˣᵢ⟩ = Σ[2i-1, 2i]  (from _sx_expect)
    mx = [Σ[2i - 1, 2i] for i in 1:N]
    # diagonal: ⟨(σˣᵢ)²⟩_c = 1 − ⟨σˣᵢ⟩²
    s = sum(1.0 - mx[i]^2 for i in 1:N)
    # off-diagonal: Wick contraction of -iγ_{2i-1}γ_{2i} · -iγ_{2j-1}γ_{2j}
    for i in 1:N, j in (i + 1):N
        cij = -Σ[2i - 1, 2j - 1] * Σ[2i, 2j] + Σ[2i - 1, 2j] * Σ[2i, 2j - 1]
        s += 2 * cij
    end
    return β * s / N
end

"""
    _tfim_transverse_obc(quantity, N, J, h, β) -> Float64

Compute `m_x` or `χ_xx` per site for OBC finite N by direct site-resolved BdG
expectation.  Uses the Majorana covariance formula

    Σ(β) = -i tanh(β/2 · i h_BdG)

(see `TFIM_dynamics.jl`) and identifies `⟨σˣ_i⟩ = Σ[2i-1, 2i]`.

The transverse susceptibility is computed via `_xx_uniform_susceptibility`
(exact Wick contraction, no numerical differentiation).
"""
function _tfim_transverse_obc(quantity::Symbol, N::Int, J::Float64, h::Float64, β::Real)
    hmat = _majorana_ham(N, J, h)
    Σ = _majorana_thermal_covariance(hmat, β)
    if quantity === :transverse_magnetization
        return sum(_sx_expect(Σ, i) for i in 1:N) / N
    else  # :transverse_susceptibility
        return _xx_uniform_susceptibility(N, J, h, β)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch dispatch
# ═══════════════════════════════════════════════════════════════════════════════

# (quantity struct, internal symbol) pairs — the internal helpers
# `_tfim_thermo_infinite` / `_tfim_thermo_obc` still key off the symbol.
# The concrete-struct dispatch below is meta-programmed over this list
# so adding a thermal quantity only requires adding one row + wiring
# the inner Symbol branch.
const _TFIM_THERMAL_METHODS = (
    (FreeEnergy, :free_energy),
    (ThermalEntropy, :entropy),
    (SpecificHeat, :specific_heat),
    (MagnetizationX, :transverse_magnetization),
    (SusceptibilityXX, :transverse_susceptibility),
    (NMRSpinRelaxationRate, :nmr_relaxation),
)

for (QTy, qsym) in _TFIM_THERMAL_METHODS
    @eval begin
        """
            fetch(model::TFIM, ::$($QTy), ::Infinite; beta::Real, kwargs...)

        Per-site $($(string(qsym))) of the TFIM in the thermodynamic limit at
        inverse temperature `beta`.  Uses adaptive Gauss-Kronrod quadrature
        over the BdG dispersion `Λ(k) = 2√(J² + h² − 2Jh cos k)`.
        """
        function fetch(
            model::TFIM,
            ::$QTy,
            ::Infinite;
            scheme::Symbol=:canonical,
            beta::Real,
            kwargs...,
        )
            scheme === :canonical &&
                return _tfim_thermo_infinite($(QuoteNode(qsym)), model.J, model.h, beta; kwargs...)
            return _tfim_thermo_infinite_scheme(model, $QTy(), Val(scheme); beta=beta, kwargs...)
        end

        """
            fetch(model::TFIM, ::$($QTy), bc::OBC; beta::Real, kwargs...)

        Per-site $($(string(qsym))) of the OBC TFIM with `N = bc.N` sites at
        inverse temperature `beta`.  Computed exactly via the BdG
        diagonalisation.
        """
        function fetch(model::TFIM, ::$QTy, bc::OBC; beta::Real, kwargs...)
            N = _bc_size(bc, kwargs)
            return _tfim_thermo_obc($(QuoteNode(qsym)), N, model.J, model.h, beta)
        end
    end
end

# ── Non-canonical (approximation) schemes of Infinite thermal quantities ──
# The bare / `scheme=:canonical` fetch stays the exact closed form above; a
# `scheme=:high_T` etc. routes here.  Decision C: the multi-definition selector
# lives on the native method, so `fetch(m, q, Infinite())` reproduces the exact
# row while `definitions(m, q, Infinite())` lists the registered approximations.
function _tfim_thermo_infinite_scheme(m::TFIM, q, ::Val{S}; beta) where {S}
    return error(
        "TFIM $(typeof(q)) Infinite: no scheme :$(S) " *
        "(only :canonical + registered approximations)",
    )
end

# High-temperature expansion of the free energy density:
#   f/N = -ln2/β - (β/2)(J² + h²) + O(β³),  valid for βJ ≪ 1, βh ≪ 1.
# This is the small-β limit of the exact free energy (the :canonical row).
function _tfim_thermo_infinite_scheme(m::TFIM, ::FreeEnergy, ::Val{:high_T}; beta)
    return -log(2) / beta - (beta / 2) * (m.J^2 + m.h^2)
end
