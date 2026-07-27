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
# Quantities exposed via `fetch` (the kernels dispatch on the quantity type):
#
#   FreeEnergy                f(β)        = -T log Z / N
#   ThermalEntropy            s(β)        = β (ε - f)
#   SpecificHeat              c_v(β)      = ∂ε/∂T
#   Magnetization{:x}         m_x(β)      = ⟨σˣ_i⟩
#   Susceptibility{(:x,:x)}   χ_xx(β)     = β · Var(Σᵢ σˣᵢ) / N
#
# All five are implemented for both `OBC` (per-site, exact at finite N) and
# `Infinite` (per-site, exact in the thermodynamic limit).
# `NMRSpinRelaxationRate` is available for both as well, with a Lorentzian
# broadening `eta`.
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

# ── Type-dispatched integrands ───────────────────────────────────────────────
#
# The quantity is already a TYPE in the AbstractQAtlas vocabulary (`FreeEnergy`,
# `Magnetization{:x}`, `Susceptibility{(:x,:x)}`, …).  Selecting the integrand
# from a `Symbol` threw that static information away: `integrand` was inferred
# as a `Union` of five anonymous closure types, so the whole adaptive-quadrature
# machinery was instantiated once per branch of that union — and again for every
# distinct forwarded-`kwargs` NamedTuple type and every element type.  Carrying
# the quantity type through to the kernel restores the dispatch the vocabulary
# already provides and makes each kernel type-stable.
#
# A named struct is used instead of a closure so the integrand type is shared by
# every call site rather than minted afresh at each syntactic location.
#
# The integrand expressions themselves are unchanged, and `quadgk` is still
# called with the same limits and tolerance, so the returned values are
# bit-identical to the Symbol-dispatched implementation.

# Each field keeps its OWN type, as the closures these replaced captured them.
# Promoting them to a common type drags the model parameters up to whatever `β`
# is, and `β` arrives as a `ForwardDiff.Dual` whenever the derived-input
# suppliers differentiate through `fetch` — so the whole quadrature would run in
# Dual arithmetic for parameters carrying no derivative information, which is
# the opposite of the point of the type-dispatch sweep.  (In SSH, whose
# dispersion is declared with `Float64` parameters, the same promotion was an
# outright MethodError — QAtlas #770, `shard s15`.)
struct _TFIMIntegrand{Q,TJ<:Real,TH<:Real,TB<:Real}
    J::TJ
    h::TH
    β::TB
end

function _TFIMIntegrand{Q}(J::Real, h::Real, β::Real) where {Q}
    return _TFIMIntegrand{Q,typeof(J),typeof(h),typeof(β)}(J, h, β)
end

# f = -(1/πβ) ∫ log(2 cosh(βΛ/2)) dk
@inline (g::_TFIMIntegrand{FreeEnergy})(k) =
    _logcosh2(g.β * _tfim_dispersion(k, g.J, g.h) / 2)

# s = (1/π) ∫ [log(2 cosh(βΛ/2)) - (βΛ/2) tanh(βΛ/2)] dk
@inline function (g::_TFIMIntegrand{ThermalEntropy})(k)
    x = g.β * _tfim_dispersion(k, g.J, g.h) / 2
    return _logcosh2(x) - x * tanh(x)
end

# c_v = (1/π) ∫ (βΛ/2)² sech²(βΛ/2) dk
@inline function (g::_TFIMIntegrand{SpecificHeat})(k)
    x = g.β * _tfim_dispersion(k, g.J, g.h) / 2
    return x^2 * sech(x)^2
end

# m_x = (2/π) ∫ ((h - J cos k)/Λ) tanh(βΛ/2) dk
@inline function (g::_TFIMIntegrand{MagnetizationX})(k)
    A = g.h - g.J * cos(k)
    Λk = _tfim_dispersion(k, g.J, g.h)
    return (A / Λk) * tanh(g.β * Λk / 2)
end

# χ_xx = (2/π) ∫ [ (1/Λ - 4A²/Λ³) tanh(βΛ/2) + (2β A²/Λ²) sech²(βΛ/2) ] dk
@inline function (g::_TFIMIntegrand{SusceptibilityXX})(k)
    A = g.h - g.J * cos(k)
    Λk = _tfim_dispersion(k, g.J, g.h)
    return (1 / Λk - 4 * A^2 / Λk^3) * tanh(g.β * Λk / 2) +
           (2 * g.β * A^2 / Λk^2) * sech(g.β * Λk / 2)^2
end

# The single place the Gauss-Kronrod machinery is instantiated for these paths.
_tfim_quad(g::_TFIMIntegrand) = first(quadgk(g, 0.0, π; rtol=1e-10))

"""
    _tfim_thermo_infinite(quantity, J, h, β) -> Real

Per-site thermodynamic potential of the infinite TFIM at inverse temperature
`β`, dispatched on the concrete quantity type.  The integrals are evaluated by
adaptive Gauss-Kronrod quadrature over the BdG dispersion.
"""
function _tfim_thermo_infinite end

function _tfim_thermo_infinite(::FreeEnergy, J::Real, h::Real, β::Real)
    return -_tfim_quad(_TFIMIntegrand{FreeEnergy}(J, h, β)) / (π * β)
end

function _tfim_thermo_infinite(::ThermalEntropy, J::Real, h::Real, β::Real)
    return _tfim_quad(_TFIMIntegrand{ThermalEntropy}(J, h, β)) / π
end

function _tfim_thermo_infinite(::SpecificHeat, J::Real, h::Real, β::Real)
    return _tfim_quad(_TFIMIntegrand{SpecificHeat}(J, h, β)) / π
end

function _tfim_thermo_infinite(::MagnetizationX, J::Real, h::Real, β::Real)
    return (2 / π) * _tfim_quad(_TFIMIntegrand{MagnetizationX}(J, h, β))
end

function _tfim_thermo_infinite(::SusceptibilityXX, J::Real, h::Real, β::Real)
    return (2 / π) * _tfim_quad(_TFIMIntegrand{SusceptibilityXX}(J, h, β))
end

function _tfim_thermo_infinite(::NMRSpinRelaxationRate, J::Real, h::Real, β::Real, η::Real)
    η > 0 || throw(DomainError(η, "TFIM NMRSpinRelaxationRate requires η > 0; got η = $η."))
    β > 0 || throw(DomainError(β, "TFIM NMRSpinRelaxationRate requires β > 0; got β = $β."))
    return _tfim_nmr_relaxation_infinite(J, h, β, η)
end

# Fermi occupation of a BdG mode, with the λ = 0 zero mode at half filling.
# `oftype` keeps the branch type-stable (for `Float64` it is the literal `0.5`
# the previous inline expression returned).
@inline function _fermi(β::Real, λ::Real)
    e = exp(-β * λ)
    return λ > 0 ? e / (1.0 + e) : oftype(e, 0.5)
end

# The NMR rate is a nested quadrature.  As closures the inner integrand type
# was minted inside the outer closure's own specialization, so the inner
# Gauss-Kronrod machinery was instantiated once per specialization of the outer
# one.  Naming both integrands breaks that product.
struct _TFIMNMRInner{TJ<:Real,TH<:Real,TB<:Real,TE<:Real,TL<:Real,TF<:Real}
    J::TJ
    h::TH
    β::TB
    η::TE
    λ1::TL
    f1::TF
end

@inline function (g::_TFIMNMRInner)(k2)
    λ2 = _tfim_dispersion(k2, g.J, g.h)
    f2 = _fermi(g.β, λ2)
    lorentz = g.η / (π * ((g.λ1 - λ2)^2 + g.η^2))
    return g.f1 * (1.0 - f2) * lorentz
end

struct _TFIMNMROuter{TJ<:Real,TH<:Real,TB<:Real,TE<:Real}
    J::TJ
    h::TH
    β::TB
    η::TE
end

@inline function (g::_TFIMNMROuter)(k1)
    λ1 = _tfim_dispersion(k1, g.J, g.h)
    f1 = _fermi(g.β, λ1)
    inner = _TFIMNMRInner(g.J, g.h, g.β, g.η, λ1, f1)
    return first(quadgk(inner, 0.0, π; rtol=1e-6))
end

function _tfim_nmr_relaxation_infinite(J::Real, h::Real, β::Real, η::Real)
    val_outer, _ = quadgk(_TFIMNMROuter(J, h, β, η), 0.0, π; rtol=1e-6)
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
function _tfim_thermo_obc end

function _tfim_thermo_obc(::FreeEnergy, N::Int, J::Float64, h::Float64, β::Real)
    Λ = _tfim_bdg_spectrum(N, J, h)
    # f/N = -(1/Nβ) Σ log(2 cosh(βΛ/2))
    return -sum(λ -> _logcosh2(β * λ / 2), Λ) / (N * β)
end

function _tfim_thermo_obc(::ThermalEntropy, N::Int, J::Float64, h::Float64, β::Real)
    Λ = _tfim_bdg_spectrum(N, J, h)
    return sum(Λ) do λ
        x = β * λ / 2
        return _logcosh2(x) - x * tanh(x)
    end / N
end

function _tfim_thermo_obc(::SpecificHeat, N::Int, J::Float64, h::Float64, β::Real)
    Λ = _tfim_bdg_spectrum(N, J, h)
    return sum(Λ) do λ
        x = β * λ / 2
        return x^2 * sech(x)^2
    end / N
end

function _tfim_thermo_obc(q::MagnetizationX, N::Int, J::Float64, h::Float64, β::Real)
    return _tfim_transverse_obc(q, N, J, h, β)
end

function _tfim_thermo_obc(q::SusceptibilityXX, N::Int, J::Float64, h::Float64, β::Real)
    return _tfim_transverse_obc(q, N, J, h, β)
end

function _tfim_thermo_obc(
    ::NMRSpinRelaxationRate, N::Int, J::Float64, h::Float64, β::Real, η::Real
)
    η > 0 || throw(DomainError(η, "TFIM NMRSpinRelaxationRate requires η > 0; got η = $η."))
    β > 0 || throw(DomainError(β, "TFIM NMRSpinRelaxationRate requires β > 0; got β = $β."))
    Λ = _tfim_bdg_spectrum(N, J, h)
    s = 0.0
    for λ1 in Λ
        f1 = _fermi(β, λ1)
        for λ2 in Λ
            f2 = _fermi(β, λ2)
            lorentz = η / (π * ((λ1 - λ2)^2 + η^2))
            s += f1 * (1.0 - f2) * lorentz
        end
    end
    return s / N^2
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
function _tfim_transverse_obc(::MagnetizationX, N::Int, J::Float64, h::Float64, β::Real)
    hmat = _majorana_ham(N, J, h)
    Σ = _majorana_thermal_covariance(hmat, β)
    return sum(_sx_expect(Σ, i) for i in 1:N) / N
end

function _tfim_transverse_obc(::SusceptibilityXX, N::Int, J::Float64, h::Float64, β::Real)
    return _xx_uniform_susceptibility(N, J, h, β)
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch dispatch
# ═══════════════════════════════════════════════════════════════════════════════

# (quantity type, human-readable name) pairs.  The kernels
# `_tfim_thermo_infinite` / `_tfim_thermo_obc` dispatch on the quantity type
# itself, so this list only drives the generated `fetch` methods and their
# docstrings; adding a thermal quantity means adding one row plus the
# corresponding kernel method.
#
# `NMRSpinRelaxationRate` is not in the list: it takes an extra broadening
# parameter and so gets hand-written `fetch` methods below.
const _TFIM_THERMAL_METHODS = (
    (FreeEnergy, :free_energy),
    (ThermalEntropy, :entropy),
    (SpecificHeat, :specific_heat),
    (MagnetizationX, :transverse_magnetization),
    (SusceptibilityXX, :transverse_susceptibility),
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
            q::$QTy,
            ::Infinite;
            scheme::Symbol=:canonical,
            beta::Real,
            kwargs...,
        )
            scheme === :canonical || return _tfim_thermo_infinite_scheme(
                model, q, Val(scheme); beta=beta, kwargs...
            )
            return _tfim_thermo_infinite(q, model.J, model.h, beta)
        end

        """
            fetch(model::TFIM, ::$($QTy), bc::OBC; beta::Real, kwargs...)

        Per-site $($(string(qsym))) of the OBC TFIM with `N = bc.N` sites at
        inverse temperature `beta`.  Computed exactly via the BdG
        diagonalisation.
        """
        function fetch(model::TFIM, q::$QTy, bc::OBC; beta::Real, kwargs...)
            N = _bc_size(bc, kwargs)
            return _tfim_thermo_obc(q, N, model.J, model.h, beta)
        end
    end
end

# ── NMR spin relaxation rate ────────────────────────────────────────────────
#
# The Lorentzian broadening is taken as a named keyword rather than pulled out
# of the forwarded `kwargs`, so the numeric kernels are not specialized on the
# caller's keyword NamedTuple type.  `Float64(eta)` reproduces the conversion
# the previous `Float64(get(kwargs, :eta, 0.1))` performed.
const _TFIM_DEFAULT_ETA = 0.1

"""
    fetch(model::TFIM, ::NMRSpinRelaxationRate, ::Infinite; beta::Real, eta::Real=0.1, kwargs...)

Per-site NMR spin relaxation rate `1/T₁` of the TFIM in the thermodynamic limit,
from the Lorentzian-broadened two-quasiparticle scattering integral.
"""
function fetch(
    model::TFIM,
    q::NMRSpinRelaxationRate,
    ::Infinite;
    scheme::Symbol=:canonical,
    beta::Real,
    eta::Real=_TFIM_DEFAULT_ETA,
    kwargs...,
)
    scheme === :canonical || return _tfim_thermo_infinite_scheme(
        model, q, Val(scheme); beta=beta, eta=eta, kwargs...
    )
    return _tfim_thermo_infinite(q, model.J, model.h, beta, Float64(eta))
end

"""
    fetch(model::TFIM, ::NMRSpinRelaxationRate, bc::OBC; beta::Real, eta::Real=0.1, kwargs...)

Per-site NMR spin relaxation rate `1/T₁` of the OBC TFIM with `N = bc.N` sites,
summed over the exact BdG quasiparticle spectrum.
"""
function fetch(
    model::TFIM,
    q::NMRSpinRelaxationRate,
    bc::OBC;
    beta::Real,
    eta::Real=_TFIM_DEFAULT_ETA,
    kwargs...,
)
    N = _bc_size(bc, kwargs)
    return _tfim_thermo_obc(q, N, model.J, model.h, beta, Float64(eta))
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
