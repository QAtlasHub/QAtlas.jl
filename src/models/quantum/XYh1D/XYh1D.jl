# XYh1D — anisotropic XY chain in a transverse field (Lieb-Schultz-Mattis 1961).
#
# Hamiltonian (LSM 1961 convention, ferromagnetic sign for the XY exchanges):
#
#     H = -Σ_i ( J_x σ^x_i σ^x_{i+1} + J_y σ^y_i σ^y_{i+1} )
#         - h Σ_i σ^z_i,        J_x, J_y > 0,  h ∈ ℝ.
#
# Solved exactly via Jordan-Wigner + Bogoliubov-de Gennes (BdG) transformation.
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Pauli σ (this file)
#   Observable:  Spin S = σ/2  (QAtlas-wide spin convention; see docs/src/conventions.md)

using LinearAlgebra: eigvals, Symmetric, Hermitian, eigen, Diagonal, I
using QuadGK: quadgk

"""
    XYh1D(; Jx::Real = 1.0, Jy::Real = 1.0, h::Real = 0.0) <: AbstractQAtlasModel

Anisotropic XY chain in a transverse field (Lieb-Schultz-Mattis 1961):

    H = -Σ_i ( Jx σ^x_i σ^x_{i+1} + Jy σ^y_i σ^y_{i+1} ) - h Σ_i σ^z_i.

Requires `Jx > 0` and `Jy > 0`.
"""
struct XYh1D <: AbstractQAtlasModel
    Jx::Float64
    Jy::Float64
    h::Float64
    function XYh1D(Jx::Real, Jy::Real, h::Real)
        Jx > 0 || throw(DomainError(Jx, "XYh1D requires Jx > 0; got Jx = $Jx."))
        Jy > 0 || throw(DomainError(Jy, "XYh1D requires Jy > 0; got Jy = $Jy."))
        return new(Float64(Jx), Float64(Jy), Float64(h))
    end
end
XYh1D(; Jx::Real=1.0, Jy::Real=1.0, h::Real=0.0) = XYh1D(Jx, Jy, h)

# ═══════════════════════════════════════════════════════════════════════════════
# Internal: dispersion and BdG spectrum
# ═══════════════════════════════════════════════════════════════════════════════

@inline _xyh1d_dispersion(k::Real, Jx::Real, Jy::Real, h::Real) =
    2.0 * sqrt((h - (Jx + Jy) * cos(k))^2 + (Jx - Jy)^2 * sin(k)^2)

"""
    _xyh1d_bdg_spectrum(N, Jx, Jy, h) -> Vector{Float64}

Return the N positive BdG quasiparticle energies Λₙ > 0 for the OBC XYh1D
with N sites.
"""
function _xyh1d_bdg_spectrum(N::Int, Jx::Real, Jy::Real, h::Real)::Vector{Float64}
    A = zeros(N, N)
    for i in 1:N
        A[i, i] = 2.0 * h
    end
    for i in 1:(N - 1)
        A[i, i + 1] = -(Jx + Jy)
        A[i + 1, i] = -(Jx + Jy)
    end

    B = zeros(N, N)
    for i in 1:(N - 1)
        B[i, i + 1] = Jx - Jy
        B[i + 1, i] = -(Jx - Jy)
    end

    H_bdg = [A B; -B -A]
    vals = eigvals(Symmetric(H_bdg))
    return sort!(filter(v -> v > 1e-10, vals))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Internal: Majorana Hamiltonian and Covariance
# ═══════════════════════════════════════════════════════════════════════════════

function _xyh1d_majorana_ham(N::Int, Jx::Real, Jy::Real, h::Real)
    M = zeros(Float64, 2N, 2N)
    Jx_f, Jy_f, h_f = Float64(Jx), Float64(Jy), Float64(h)
    @inbounds for i in 1:N
        M[2i - 1, 2i] = 2h_f
        M[2i, 2i - 1] = -2h_f
    end
    @inbounds for i in 1:(N - 1)
        M[2i, 2i + 1] = 2Jx_f
        M[2i + 1, 2i] = -2Jx_f
        M[2i - 1, 2i + 2] = -2Jy_f
        M[2i + 2, 2i - 1] = 2Jy_f
    end
    return M
end

function _xyh1d_majorana_covariance_gs(h::AbstractMatrix{<:Real})
    M = im .* h
    F = eigen(Hermitian((M + M') / 2))
    s = [λ > 0 ? 1.0 : (λ < 0 ? -1.0 : 0.0) for λ in F.values]
    sM = F.vectors * Diagonal(s) * F.vectors'
    Σ = real(-im .* sM)
    return (Σ - Σ') / 2
end

function _xyh1d_majorana_thermal_covariance(h::AbstractMatrix{<:Real}, β::Real)
    isinf(β) && return _xyh1d_majorana_covariance_gs(h)
    M = im .* h
    F = eigen(Hermitian((M + M') / 2))
    s = tanh.((β / 2) .* F.values)
    sM = F.vectors * Diagonal(s) * F.vectors'
    Σ = real(-im .* sM)
    return (Σ - Σ') / 2
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy granularity convention
# ═══════════════════════════════════════════════════════════════════════════════

native_energy_granularity(::XYh1D, ::OBC) = :total
native_energy_granularity(::XYh1D, ::Infinite) = :per_site

# ═══════════════════════════════════════════════════════════════════════════════
# Mass Gap — Infinite and OBC
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::XYh1D, ::MassGap, ::Infinite) -> Float64

Single-particle Bogoliubov gap of the LSM/Pfeuty XY chain in a transverse field.
"""
function fetch(
    m::XYh1D, ::MassGap, ::Infinite; Jx::Real=m.Jx, Jy::Real=m.Jy, h::Real=m.h, kwargs...
)
    x0 = (Jx + Jy) * h / (4.0 * Jx * Jy)
    min_val = if -1.0 <= x0 <= 1.0
        (Jx - Jy)^2 * (1.0 - h^2 / (4.0 * Jx * Jy))
    elseif x0 > 1.0
        (h - (Jx + Jy))^2
    else
        (h + Jx + Jy)^2
    end
    return 2.0 * sqrt(max(0.0, min_val))
end

"""
    fetch(model::XYh1D, ::MassGap, bc::OBC) -> Float64
"""
function fetch(m::XYh1D, ::MassGap, bc::OBC; kwargs...)
    N = _bc_size(bc, kwargs)
    Λ = _xyh1d_bdg_spectrum(N, m.Jx, m.Jy, m.h)
    return Λ[1]
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy: Infinite (thermodynamic limit)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::XYh1D, ::Energy{:per_site}, ::Infinite; beta, betas) -> Float64 or Vector{Float64}
"""
function fetch(
    m::XYh1D,
    ::Energy{:per_site},
    ::Infinite;
    Jx::Real=m.Jx,
    Jy::Real=m.Jy,
    h::Real=m.h,
    beta::Union{Real,Nothing}=nothing,
    betas::Union{AbstractVector{<:Real},Nothing}=nothing,
    kwargs...,
)
    _energy_at_beta =
        β -> begin
            result, _ = quadgk(
                k -> begin
                    Λk = _xyh1d_dispersion(k, Jx, Jy, h)
                    (Λk / 2.0) * tanh(β * Λk / 2.0)
                end, 0.0, π; rtol=1e-10
            )
            -(1.0 / π) * result
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
# Energy: OBC finite-N
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::XYh1D, ::Energy{:total}, bc::OBC; beta, betas) -> Float64 or Vector{Float64}
"""
function fetch(
    m::XYh1D,
    ::Energy{:total},
    bc::OBC;
    beta::Union{Real,Nothing}=nothing,
    betas::Union{AbstractVector{<:Real},Nothing}=nothing,
    kwargs...,
)
    N = _bc_size(bc, kwargs)
    Λ = _xyh1d_bdg_spectrum(N, m.Jx, m.Jy, m.h)
    if betas !== nothing
        return [-sum(λ -> (λ / 2.0) * tanh(β * λ / 2.0), Λ) for β in betas]
    elseif beta !== nothing
        return -sum(λ -> (λ / 2.0) * tanh(beta * λ / 2.0), Λ)
    else
        # Ground state
        return -sum(Λ) / 2.0
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Infinite-chain Thermodynamics via integration over dispersion
# ═══════════════════════════════════════════════════════════════════════════════

@inline function _xyh1d_logcosh2(x::Real)
    a = abs(x)
    return a + log1p(exp(-2.0 * a))
end

# ── Type-dispatched integrands ───────────────────────────────────────────────
#
# The quantity reaches these kernels as a type in the AbstractQAtlas vocabulary;
# selecting the integrand from a `Symbol` left it a union of five anonymous
# closure types at the `quadgk` call site.  Named integrand structs dispatched
# on the quantity keep the static information the vocabulary already carries.
# The integrand expressions and the `quadgk` call are unchanged, so the values
# are unchanged.

struct _XYh1DIntegrand{Q,T<:Real}
    Jx::T
    Jy::T
    h::T
    β::T
end

function _XYh1DIntegrand{Q}(Jx::Real, Jy::Real, h::Real, β::Real) where {Q}
    a, b, c, d = promote(Jx, Jy, h, β)
    return _XYh1DIntegrand{Q,typeof(a)}(a, b, c, d)
end

@inline (g::_XYh1DIntegrand{FreeEnergy})(k) =
    _xyh1d_logcosh2(g.β * _xyh1d_dispersion(k, g.Jx, g.Jy, g.h) / 2.0)

@inline function (g::_XYh1DIntegrand{ThermalEntropy})(k)
    x = g.β * _xyh1d_dispersion(k, g.Jx, g.Jy, g.h) / 2.0
    return _xyh1d_logcosh2(x) - x * tanh(x)
end

@inline function (g::_XYh1DIntegrand{SpecificHeat})(k)
    x = g.β * _xyh1d_dispersion(k, g.Jx, g.Jy, g.h) / 2.0
    return x^2 * sech(x)^2
end

@inline function (g::_XYh1DIntegrand{MagnetizationZ})(k)
    Λk = _xyh1d_dispersion(k, g.Jx, g.Jy, g.h)
    A = g.h - (g.Jx + g.Jy) * cos(k)
    return (2.0 * A / Λk) * tanh(g.β * Λk / 2.0)
end

# (2/Λ - 8A²/Λ³) tanh(βΛ/2) + (4β A²/Λ²) sech²(βΛ/2)
@inline function (g::_XYh1DIntegrand{SusceptibilityZZ})(k)
    A = g.h - (g.Jx + g.Jy) * cos(k)
    Λk = _xyh1d_dispersion(k, g.Jx, g.Jy, g.h)
    return (2.0 / Λk - 8.0 * A^2 / Λk^3) * tanh(g.β * Λk / 2.0) +
           (4.0 * g.β * A^2 / Λk^2) * sech(g.β * Λk / 2.0)^2
end

_xyh1d_quad(g::_XYh1DIntegrand) = first(quadgk(g, 0.0, π; rtol=1e-10))

"""
    _xyh1d_thermo_infinite(quantity, Jx, Jy, h, β) -> Real

Per-site thermodynamic potential of the infinite XY chain in a transverse
field, dispatched on the concrete quantity type.
"""
function _xyh1d_thermo_infinite end

function _xyh1d_thermo_infinite(::FreeEnergy, Jx::Real, Jy::Real, h::Real, β::Real)
    return -_xyh1d_quad(_XYh1DIntegrand{FreeEnergy}(Jx, Jy, h, β)) / (π * β)
end

function _xyh1d_thermo_infinite(::ThermalEntropy, Jx::Real, Jy::Real, h::Real, β::Real)
    return _xyh1d_quad(_XYh1DIntegrand{ThermalEntropy}(Jx, Jy, h, β)) / π
end

function _xyh1d_thermo_infinite(::SpecificHeat, Jx::Real, Jy::Real, h::Real, β::Real)
    return _xyh1d_quad(_XYh1DIntegrand{SpecificHeat}(Jx, Jy, h, β)) / π
end

function _xyh1d_thermo_infinite(::MagnetizationZ, Jx::Real, Jy::Real, h::Real, β::Real)
    return (1.0 / π) * _xyh1d_quad(_XYh1DIntegrand{MagnetizationZ}(Jx, Jy, h, β))
end

function _xyh1d_thermo_infinite(::SusceptibilityZZ, Jx::Real, Jy::Real, h::Real, β::Real)
    return (1.0 / π) * _xyh1d_quad(_XYh1DIntegrand{SusceptibilityZZ}(Jx, Jy, h, β))
end

# ═══════════════════════════════════════════════════════════════════════════════
# OBC Finite-size Thermodynamics via BdG
# ═══════════════════════════════════════════════════════════════════════════════

function _xyh1d_zz_uniform_susceptibility(N::Int, Jx::Real, Jy::Real, h::Real, β::Real)
    hmat = _xyh1d_majorana_ham(N, Jx, Jy, h)
    Σ = _xyh1d_majorana_thermal_covariance(hmat, β)
    mx = [Σ[2i - 1, 2i] for i in 1:N]
    s = sum(1.0 - mx[i]^2 for i in 1:N)
    for i in 1:N, j in (i + 1):N
        cij = -Σ[2i - 1, 2j - 1] * Σ[2i, 2j] + Σ[2i - 1, 2j] * Σ[2i, 2j - 1]
        s += 2.0 * cij
    end
    return β * s / N
end

"""
    _xyh1d_thermo_obc(quantity, N, Jx, Jy, h, β) -> Real

Per-site thermodynamic potential of the OBC finite-N XY chain, dispatched on
the concrete quantity type.
"""
function _xyh1d_thermo_obc end

function _xyh1d_thermo_obc(::FreeEnergy, N::Int, Jx::Real, Jy::Real, h::Real, β::Real)
    Λ = _xyh1d_bdg_spectrum(N, Jx, Jy, h)
    return -sum(λ -> _xyh1d_logcosh2(β * λ / 2.0), Λ) / (N * β)
end

function _xyh1d_thermo_obc(::ThermalEntropy, N::Int, Jx::Real, Jy::Real, h::Real, β::Real)
    Λ = _xyh1d_bdg_spectrum(N, Jx, Jy, h)
    return sum(Λ) do λ
        x = β * λ / 2.0
        return _xyh1d_logcosh2(x) - x * tanh(x)
    end / N
end

function _xyh1d_thermo_obc(::SpecificHeat, N::Int, Jx::Real, Jy::Real, h::Real, β::Real)
    Λ = _xyh1d_bdg_spectrum(N, Jx, Jy, h)
    return sum(Λ) do λ
        x = β * λ / 2.0
        return x^2 * sech(x)^2
    end / N
end

function _xyh1d_thermo_obc(::MagnetizationZ, N::Int, Jx::Real, Jy::Real, h::Real, β::Real)
    hmat = _xyh1d_majorana_ham(N, Jx, Jy, h)
    Σ = _xyh1d_majorana_thermal_covariance(hmat, β)
    return sum(Σ[2i - 1, 2i] for i in 1:N) / N
end

function _xyh1d_thermo_obc(::SusceptibilityZZ, N::Int, Jx::Real, Jy::Real, h::Real, β::Real)
    return _xyh1d_zz_uniform_susceptibility(N, Jx, Jy, h, β)
end

# ═══════════════════════════════════════════════════════════════════════════════
# site-local equilibrium observables
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::XYh1D, ::LocalMagnetization{:z}, bc::OBC; beta, kwargs...)
"""
function fetch(model::XYh1D, ::LocalMagnetization{:z}, bc::OBC; beta::Float64, kwargs...)
    N = _bc_size(bc, kwargs)
    hmat = _xyh1d_majorana_ham(N, model.Jx, model.Jy, model.h)
    Σ = _xyh1d_majorana_thermal_covariance(hmat, beta)
    return Float64[Σ[2i - 1, 2i] for i in 1:N]
end

"""
    fetch(model::XYh1D, ::LocalMagnetization{:x}, bc::OBC; beta, kwargs...)
"""
function fetch(::XYh1D, ::LocalMagnetization{:x}, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    return zeros(Float64, N)
end

"""
    fetch(model::XYh1D, ::LocalMagnetization{:y}, bc::OBC; beta, kwargs...)
"""
function fetch(::XYh1D, ::LocalMagnetization{:y}, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    return zeros(Float64, N)
end

"""
    fetch(model::XYh1D, ::EnergyLocal, bc::OBC; beta, kwargs...)
"""
function fetch(model::XYh1D, ::EnergyLocal, bc::OBC; beta::Float64, kwargs...)
    N = _bc_size(bc, kwargs)
    hmat = _xyh1d_majorana_ham(N, model.Jx, model.Jy, model.h)
    Σ = _xyh1d_majorana_thermal_covariance(hmat, beta)

    bonds = Float64[
        -model.Jx * Σ[2i, 2i + 1] + model.Jy * Σ[2i - 1, 2i + 2] for i in 1:(N - 1)
    ]

    ε = Vector{Float64}(undef, N)
    @inbounds for i in 1:N
        left = i > 1 ? bonds[i - 1] : 0.0
        right = i < N ? bonds[i] : 0.0
        ε[i] = 0.5 * (left + right) - model.h * Σ[2i - 1, 2i]
    end
    return ε
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch dispatch for thermal potentials
# ═══════════════════════════════════════════════════════════════════════════════

const _XYH1D_THERMAL_METHODS = (
    (FreeEnergy, :free_energy),
    (ThermalEntropy, :entropy),
    (SpecificHeat, :specific_heat),
    (MagnetizationZ, :transverse_magnetization),
    (SusceptibilityZZ, :transverse_susceptibility),
)

for (QTy, _) in _XYH1D_THERMAL_METHODS
    @eval begin
        function fetch(model::XYh1D, q::$QTy, ::Infinite; beta::Real, kwargs...)
            return _xyh1d_thermo_infinite(q, model.Jx, model.Jy, model.h, beta)
        end

        function fetch(model::XYh1D, q::$QTy, bc::OBC; beta::Real, kwargs...)
            N = _bc_size(bc, kwargs)
            return _xyh1d_thermo_obc(q, N, model.Jx, model.Jy, model.h, beta)
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Internal: PBC helpers — two-sector exact solution
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _xyh1d_pbc_momenta(N, sector) -> Vector{Float64}

Return the N allowed momenta for the PBC XYh1D chain.

- `:AP` (antiperiodic BC for JW fermions, **even** fermion parity, "Ramond"):
  kₙ = (2n−1)π/N,  n = 1,…,N
- `:P`  (periodic BC for JW fermions, **odd** fermion parity, "Neveu-Schwarz"):
  kₙ = 2π(n−1)/N,  n = 1,…,N
"""
@inline function _xyh1d_pbc_momenta(N::Int, sector::Symbol)
    sector === :AP && return [(2n - 1) * π / N for n in 1:N]
    return [2π * (n - 1) / N for n in 1:N]
end

"""
    _xyh1d_pbc_spectrum(N, Jx, Jy, h) -> (Λ_AP, Λ_P)

Return the N quasiparticle energies Λ(kₙ) for the AP and P sectors.
"""
function _xyh1d_pbc_spectrum(N::Int, Jx::Real, Jy::Real, h::Real)
    ks_AP = _xyh1d_pbc_momenta(N, :AP)
    ks_P = _xyh1d_pbc_momenta(N, :P)
    Λ_AP = [_xyh1d_dispersion(k, Jx, Jy, h) for k in ks_AP]
    Λ_P = [_xyh1d_dispersion(k, Jx, Jy, h) for k in ks_P]
    return Λ_AP, Λ_P
end

"""
    _xyh1d_pbc_sector_logZ(Λ, β) -> Float64

log Z for a single free-fermion sector:
  log Z_s = Σ_k log 2cosh(βΛ_k/2)
"""
@inline function _xyh1d_pbc_sector_logZ(Λ::AbstractVector, β::Real)
    return sum(λ -> _xyh1d_logcosh2(β * λ / 2.0) + log(2.0), Λ)
end

"""
    _xyh1d_pbc_logZ(N, Jx, Jy, h, β) -> Float64

log of the total PBC partition function using log-sum-exp:
  log Z_PBC = log(Z_AP + Z_P)
"""
function _xyh1d_pbc_logZ(N::Int, Jx::Real, Jy::Real, h::Real, β::Real)
    Λ_AP, Λ_P = _xyh1d_pbc_spectrum(N, Jx, Jy, h)
    a = _xyh1d_pbc_sector_logZ(Λ_AP, β)
    b = _xyh1d_pbc_sector_logZ(Λ_P, β)
    return max(a, b) + log1p(exp(-abs(a - b)))
end

"""
    _xyh1d_pbc_majorana_ham(N, Jx, Jy, h; sector=:AP) -> Matrix{Float64}

Majorana Hamiltonian for PBC XYh1D.

Extends `_xyh1d_majorana_ham` (OBC tridiagonal) with the corner bond
connecting site N back to site 1:

    AP sector (even parity): corner sign = +1
    P  sector (odd parity):  corner sign = −1

Corner bond in the OBC Majorana structure (bond i→i+1 uses indices 2i,2i+1 for
Jx-coupling and 2i-1,2i+2 for Jy-coupling).  For the N→1 bond set i=N:

    M[2N-1, 2]   ±= +2Jx   (Jx hopping, sign = ±s)
    M[2,   2N-1] ±= −2Jx
    M[2N,    1]  ±= −2Jy   (Jy hopping, sign = ±s)
    M[1,    2N]  ±= +2Jy
"""
function _xyh1d_pbc_majorana_ham(N::Int, Jx::Real, Jy::Real, h::Real; sector::Symbol=:AP)
    M = copy(_xyh1d_majorana_ham(N, Jx, Jy, h))
    s = (sector === :AP) ? 1.0 : -1.0
    @inbounds begin
        M[2N - 1, 2] += s * 2.0 * Jx
        M[2, 2N - 1] -= s * 2.0 * Jx
        M[2N, 1] -= s * 2.0 * Jy
        M[1, 2N] += s * 2.0 * Jy
    end
    return M
end

"""
    _xyh1d_thermo_pbc(quantity, N, Jx, Jy, h, β) -> Real

Per-site thermal quantity for the PBC XYh1D chain using the exact two-sector
free-fermion partition function.

The full canonical trace factorises as:
  Z_PBC = Z_AP + Z_P,    Z_s = ∏_k 2cosh(βΛₛ(k)/2)

Sector weights: w_s = Z_s / Z_total (computed via log-sum-exp).

All derived quantities follow from ∂/∂β of log Z_PBC.
"""
# Shared two-sector setup: spectra, momenta, the log-sum-exp combined log Z and
# the two macro-sector Boltzmann weights.
function _xyh1d_pbc_sectors(N::Int, Jx::Real, Jy::Real, h::Real, β::Real)
    Λ_AP, Λ_P = _xyh1d_pbc_spectrum(N, Jx, Jy, h)
    ks_AP = _xyh1d_pbc_momenta(N, :AP)
    ks_P = _xyh1d_pbc_momenta(N, :P)

    lZ_AP = _xyh1d_pbc_sector_logZ(Λ_AP, β)
    lZ_P = _xyh1d_pbc_sector_logZ(Λ_P, β)
    lZ = max(lZ_AP, lZ_P) + log1p(exp(-abs(lZ_AP - lZ_P)))
    w_AP = exp(lZ_AP - lZ)
    w_P = exp(lZ_P - lZ)
    return (Λ_AP, Λ_P, ks_AP, ks_P, lZ, w_AP, w_P)
end

"""
    _xyh1d_thermo_pbc(quantity, N, Jx, Jy, h, β) -> Real

Per-site thermodynamic potential of the PBC finite-N XY chain, dispatched on
the concrete quantity type.
"""
function _xyh1d_thermo_pbc end

function _xyh1d_thermo_pbc(::FreeEnergy, N::Int, Jx::Real, Jy::Real, h::Real, β::Real)
    _, _, _, _, lZ, _, _ = _xyh1d_pbc_sectors(N, Jx, Jy, h, β)
    return -lZ / (N * β)
end

# S_s/N = (1/N) Σ_k [ log(2cosh) − (βΛ/2) tanh(βΛ/2) ], plus the mixing entropy
# of the two macro-sectors.
function _xyh1d_thermo_pbc(::ThermalEntropy, N::Int, Jx::Real, Jy::Real, h::Real, β::Real)
    Λ_AP, Λ_P, _, _, _, w_AP, w_P = _xyh1d_pbc_sectors(N, Jx, Jy, h, β)
    _s(Λ) = sum(Λ) do λ
        x = β * λ / 2
        return _xyh1d_logcosh2(x) + log(2) - x * tanh(x)
    end / N
    _xlogx(w) = w > 1e-300 ? w * log(w) : 0.0
    S_mix = -(_xlogx(w_AP) + _xlogx(w_P)) / N
    return w_AP * _s(Λ_AP) + w_P * _s(Λ_P) + S_mix
end

# Cv/N = β² [ w_AP * (Cv_AP/N + β²(E_AP/N - Ē)²)
#           + w_P  * (Cv_P/N  + β²(E_P/N  - Ē)²) ]
# where Ē = w_AP E_AP/N + w_P E_P/N
function _xyh1d_thermo_pbc(::SpecificHeat, N::Int, Jx::Real, Jy::Real, h::Real, β::Real)
    Λ_AP, Λ_P, _, _, _, w_AP, w_P = _xyh1d_pbc_sectors(N, Jx, Jy, h, β)
    _e(Λ) = -sum(λ -> (λ / 2) * tanh(β * λ / 2), Λ) / N
    _cv(Λ) = sum(λ -> (β * λ / 2)^2 * sech(β * λ / 2)^2, Λ) / N
    E_AP, E_P = _e(Λ_AP), _e(Λ_P)
    C_AP, C_P = _cv(Λ_AP), _cv(Λ_P)
    Ē = w_AP * E_AP + w_P * E_P
    return w_AP * (C_AP + β^2 * N * (E_AP - Ē)^2) + w_P * (C_P + β^2 * N * (E_P - Ē)^2)
end

# ⟨σᶻ⟩/N = (1/(Nβ)) ∂logZ/∂h = (1/N) Σ_k (∂Λ_k/∂h) tanh(βΛ_k/2),
# with Λ = 2√(A² + C²), A = h − (Jx+Jy)cos k, C = (Jx−Jy)sin k, so ∂Λ/∂h = 4A/Λ.
# This matches the `MagnetizationZ` integrand at `Infinite`, which uses 2A/Λ over
# the half BZ [0, π]; the finite sum here runs over the full k range.
function _xyh1d_thermo_pbc(::MagnetizationZ, N::Int, Jx::Real, Jy::Real, h::Real, β::Real)
    Λ_AP, Λ_P, ks_AP, ks_P, _, w_AP, w_P = _xyh1d_pbc_sectors(N, Jx, Jy, h, β)
    _mz_sector(Λ, ks) = sum(zip(Λ, ks)) do (λ, k)
        A = h - (Jx + Jy) * cos(k)
        return (4.0 * A / λ) * tanh(β * λ / 2.0)
    end / N
    # ⟨σᶻ⟩/N = w_AP * ⟨σᶻ⟩_AP/N + w_P * ⟨σᶻ⟩_P/N
    return w_AP * _mz_sector(Λ_AP, ks_AP) + w_P * _mz_sector(Λ_P, ks_P)
end

# χ = β (⟨M²⟩ − ⟨M⟩²)/N.  Taken as a central difference of the magnetisation in
# `h`, which is the robust route here (the two-sector weights depend on `h` too,
# so the per-sector closed form is not the whole answer).
function _xyh1d_thermo_pbc(
    q::SusceptibilityZZ, N::Int, Jx::Real, Jy::Real, h::Real, β::Real
)
    δh = h * 1e-5 + 1e-8
    mz_plus = _xyh1d_thermo_pbc(MagnetizationZ(), N, Jx, Jy, h + δh, β)
    mz_minus = _xyh1d_thermo_pbc(MagnetizationZ(), N, Jx, Jy, h - δh, β)
    return (mz_plus - mz_minus) / (2δh)
end

"""
    _xyh1d_pbc_local_mz(N, Jx, Jy, h, β) -> Vector{Float64}

Site-local ⟨σᶻᵢ⟩ for PBC.  By translational invariance, all sites are equivalent;
returns a uniform vector of length N equal to the bulk ⟨σᶻ⟩/site.
"""
function _xyh1d_pbc_local_mz(N::Int, Jx::Real, Jy::Real, h::Real, β::Real)
    mz_bulk = _xyh1d_thermo_pbc(MagnetizationZ(), N, Jx, Jy, h, β)
    return fill(mz_bulk, N)
end

"""
    _xyh1d_pbc_local_energy(N, Jx, Jy, h, β) -> Vector{Float64}

Site-local energy density for PBC.  By translational invariance,
all sites share the same value E_total / N.
"""
function _xyh1d_pbc_local_energy(N::Int, Jx::Real, Jy::Real, h::Real, β::Real)
    Λ_AP, Λ_P = _xyh1d_pbc_spectrum(N, Jx, Jy, h)
    lZ_AP = _xyh1d_pbc_sector_logZ(Λ_AP, β)
    lZ_P = _xyh1d_pbc_sector_logZ(Λ_P, β)
    lZ = max(lZ_AP, lZ_P) + log1p(exp(-abs(lZ_AP - lZ_P)))
    w_AP = exp(lZ_AP - lZ)
    w_P = exp(lZ_P - lZ)
    e_AP = -sum(λ -> (λ / 2) * tanh(β * λ / 2), Λ_AP) / N
    e_P = -sum(λ -> (λ / 2) * tanh(β * λ / 2), Λ_P) / N
    e_per_site = w_AP * e_AP + w_P * e_P
    return fill(e_per_site, N)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Mass Gap — PBC
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::XYh1D, ::MassGap, bc::PBC; kwargs...) -> Float64

Lowest quasiparticle energy over both PBC sectors.

The ground-state gap is min(min(Λ_AP), min(Λ_P)), excluding any k=0, π
degeneracies at the transition point.
"""
function fetch(model::XYh1D, ::MassGap, bc::PBC; kwargs...)
    N = _bc_size(bc, kwargs)
    Λ_AP, Λ_P = _xyh1d_pbc_spectrum(N, model.Jx, model.Jy, model.h)
    return min(minimum(Λ_AP), minimum(Λ_P))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy — PBC (Phase 2, #292)
# ═══════════════════════════════════════════════════════════════════════════════

native_energy_granularity(::XYh1D, ::PBC) = :total

"""
    fetch(model::XYh1D, ::Energy{:total}, bc::PBC; beta, kwargs...) -> Float64

Ground-state (β → ∞) or thermal total energy for PBC XYh1D, via the exact
two-sector free-fermion partition function (Lieb-Schultz-Mattis 1961).
"""
function fetch(
    model::XYh1D,
    ::Energy{:total},
    bc::PBC;
    beta::Union{Real,Nothing}=nothing,
    betas::Union{AbstractVector{<:Real},Nothing}=nothing,
    kwargs...,
)
    N = _bc_size(bc, kwargs)
    Λ_AP, Λ_P = _xyh1d_pbc_spectrum(N, model.Jx, model.Jy, model.h)
    function _pbc_energy(β::Real)
        lZ_AP = _xyh1d_pbc_sector_logZ(Λ_AP, β)
        lZ_P = _xyh1d_pbc_sector_logZ(Λ_P, β)
        lZ = max(lZ_AP, lZ_P) + log1p(exp(-abs(lZ_AP - lZ_P)))
        w_AP = exp(lZ_AP - lZ)
        w_P = exp(lZ_P - lZ)
        e_AP = -sum(λ -> (λ / 2) * tanh(β * λ / 2), Λ_AP)
        e_P = -sum(λ -> (λ / 2) * tanh(β * λ / 2), Λ_P)
        return w_AP * e_AP + w_P * e_P
    end
    if betas !== nothing
        return [_pbc_energy(β) for β in betas]
    elseif beta !== nothing
        return _pbc_energy(beta)
    else
        # Ground state: β → ∞, lowest sector dominates
        e_AP_gs = -sum(Λ_AP) / 2
        e_P_gs = -sum(Λ_P) / 2
        return min(e_AP_gs, e_P_gs)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Thermodynamic potentials and magnetization — PBC (Phase 2, #292)
# ═══════════════════════════════════════════════════════════════════════════════

const _XYH1D_PBC_THERMAL_METHODS = (
    (FreeEnergy, :free_energy),
    (ThermalEntropy, :entropy),
    (SpecificHeat, :specific_heat),
    (MagnetizationZ, :transverse_magnetization),
    (SusceptibilityZZ, :transverse_susceptibility),
)

for (QTy, _) in _XYH1D_PBC_THERMAL_METHODS
    @eval begin
        function fetch(model::XYh1D, q::$QTy, bc::PBC; beta::Real, kwargs...)
            N = _bc_size(bc, kwargs)
            return _xyh1d_thermo_pbc(q, N, model.Jx, model.Jy, model.h, beta)
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Site-local observables at PBC (Phase 2, #292)
# (Translational invariance → uniform vector = bulk per-site scalar)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::XYh1D, ::LocalMagnetization{:z}, bc::PBC; beta) -> Vector{Float64}

Site-resolved ⟨σᶻ_i⟩ on the PBC chain. Translational invariance gives a
uniform vector of length N filled with the bulk MagnetizationZ value.
"""
function fetch(model::XYh1D, ::LocalMagnetization{:z}, bc::PBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    return _xyh1d_pbc_local_mz(N, model.Jx, model.Jy, model.h, beta)
end

"""
    fetch(model::XYh1D, ::LocalMagnetization{:x}, bc::PBC; beta) -> Vector{Float64}

Vanishes by Z₂ symmetry σˣ → −σˣ; returns zeros of length N.
"""
function fetch(::XYh1D, ::LocalMagnetization{:x}, bc::PBC; beta::Real, kwargs...)
    return zeros(Float64, _bc_size(bc, kwargs))
end

"""
    fetch(model::XYh1D, ::LocalMagnetization{:y}, bc::PBC; beta) -> Vector{Float64}

Vanishes by Z₂ symmetry σʸ → −σʸ; returns zeros of length N.
"""
function fetch(::XYh1D, ::LocalMagnetization{:y}, bc::PBC; beta::Real, kwargs...)
    return zeros(Float64, _bc_size(bc, kwargs))
end

"""
    fetch(model::XYh1D, ::EnergyLocal, bc::PBC; beta) -> Vector{Float64}

Site-resolved energy density on the PBC chain. Translational invariance
gives a uniform vector ε_i = E_total / N.
"""
function fetch(model::XYh1D, ::EnergyLocal, bc::PBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    return _xyh1d_pbc_local_energy(N, model.Jx, model.Jy, model.h, beta)
end
