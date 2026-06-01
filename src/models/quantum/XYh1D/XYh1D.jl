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

function _xyh1d_thermo_infinite(quantity::Symbol, Jx::Real, Jy::Real, h::Real, β::Real)
    integrand = if quantity === :free_energy
        k -> begin
            Λk = _xyh1d_dispersion(k, Jx, Jy, h)
            _xyh1d_logcosh2(β * Λk / 2.0)
        end
    elseif quantity === :entropy
        k -> begin
            Λk = _xyh1d_dispersion(k, Jx, Jy, h)
            x = β * Λk / 2.0
            _xyh1d_logcosh2(x) - x * tanh(x)
        end
    elseif quantity === :specific_heat
        k -> begin
            Λk = _xyh1d_dispersion(k, Jx, Jy, h)
            x = β * Λk / 2.0
            x^2 * sech(x)^2
        end
    elseif quantity === :transverse_magnetization
        k -> begin
            Λk = _xyh1d_dispersion(k, Jx, Jy, h)
            A = h - (Jx + Jy) * cos(k)
            (2.0 * A / Λk) * tanh(β * Λk / 2.0)
        end
    elseif quantity === :transverse_susceptibility
        k -> begin
            A = h - (Jx + Jy) * cos(k)
            Λk = _xyh1d_dispersion(k, Jx, Jy, h)
            # (2/Λ - 8A²/Λ³) tanh(βΛ/2) + (4β A²/Λ²) sech²(βΛ/2)
            (2.0 / Λk - 8.0 * A^2 / Λk^3) * tanh(β * Λk / 2.0) +
            (4.0 * β * A^2 / Λk^2) * sech(β * Λk / 2.0)^2
        end
    else
        error("Unknown thermal quantity: $quantity")
    end

    val, _ = quadgk(integrand, 0.0, π; rtol=1e-10)

    if quantity === :free_energy
        return -val / (π * β)
    elseif quantity === :transverse_magnetization || quantity === :transverse_susceptibility
        return (1.0 / π) * val
    else  # entropy, specific_heat
        return val / π
    end
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

function _xyh1d_thermo_obc(quantity::Symbol, N::Int, Jx::Real, Jy::Real, h::Real, β::Real)
    if quantity === :free_energy
        Λ = _xyh1d_bdg_spectrum(N, Jx, Jy, h)
        return -sum(λ -> _xyh1d_logcosh2(β * λ / 2.0), Λ) / (N * β)
    elseif quantity === :entropy
        Λ = _xyh1d_bdg_spectrum(N, Jx, Jy, h)
        return sum(Λ) do λ
            x = β * λ / 2.0
            _xyh1d_logcosh2(x) - x * tanh(x)
        end / N
    elseif quantity === :specific_heat
        Λ = _xyh1d_bdg_spectrum(N, Jx, Jy, h)
        return sum(λ -> begin
            x = β * λ / 2.0
            x^2 * sech(x)^2
        end, Λ) / N
    elseif quantity === :transverse_magnetization
        hmat = _xyh1d_majorana_ham(N, Jx, Jy, h)
        Σ = _xyh1d_majorana_thermal_covariance(hmat, β)
        return sum(Σ[2i - 1, 2i] for i in 1:N) / N
    elseif quantity === :transverse_susceptibility
        return _xyh1d_zz_uniform_susceptibility(N, Jx, Jy, h, β)
    else
        error("Unknown thermal quantity: $quantity")
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# site-local equilibrium observables
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::XYh1D, ::MagnetizationZLocal, bc::OBC; beta, kwargs...)
"""
function fetch(model::XYh1D, ::MagnetizationZLocal, bc::OBC; beta::Float64, kwargs...)
    N = _bc_size(bc, kwargs)
    hmat = _xyh1d_majorana_ham(N, model.Jx, model.Jy, model.h)
    Σ = _xyh1d_majorana_thermal_covariance(hmat, beta)
    return Float64[Σ[2i - 1, 2i] for i in 1:N]
end

"""
    fetch(model::XYh1D, ::MagnetizationXLocal{:equilibrium}, bc::OBC; beta, kwargs...)
"""
function fetch(::XYh1D, ::MagnetizationXLocal{:equilibrium}, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    return zeros(Float64, N)
end

"""
    fetch(model::XYh1D, ::MagnetizationYLocal, bc::OBC; beta, kwargs...)
"""
function fetch(::XYh1D, ::MagnetizationYLocal, bc::OBC; beta::Real, kwargs...)
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

for (QTy, qsym) in _XYH1D_THERMAL_METHODS
    @eval begin
        function fetch(model::XYh1D, ::$QTy, ::Infinite; beta::Real, kwargs...)
            return _xyh1d_thermo_infinite(
                $(QuoteNode(qsym)), model.Jx, model.Jy, model.h, beta
            )
        end

        function fetch(model::XYh1D, ::$QTy, bc::OBC; beta::Real, kwargs...)
            N = _bc_size(bc, kwargs)
            return _xyh1d_thermo_obc(
                $(QuoteNode(qsym)), N, model.Jx, model.Jy, model.h, beta
            )
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
function _xyh1d_thermo_pbc(quantity::Symbol, N::Int, Jx::Real, Jy::Real, h::Real, β::Real)
    Λ_AP, Λ_P = _xyh1d_pbc_spectrum(N, Jx, Jy, h)
    ks_AP = _xyh1d_pbc_momenta(N, :AP)
    ks_P = _xyh1d_pbc_momenta(N, :P)

    lZ_AP = _xyh1d_pbc_sector_logZ(Λ_AP, β)
    lZ_P = _xyh1d_pbc_sector_logZ(Λ_P, β)
    lZ = max(lZ_AP, lZ_P) + log1p(exp(-abs(lZ_AP - lZ_P)))
    w_AP = exp(lZ_AP - lZ)
    w_P = exp(lZ_P - lZ)

    if quantity === :free_energy
        return -lZ / (N * β)

    elseif quantity === :entropy
        # S_s/N = (1/N) Σ_k [ log(2cosh) − (βΛ/2) tanh(βΛ/2) ]
        _s(Λ) = sum(λ -> begin
            x = β * λ / 2;
            _xyh1d_logcosh2(x) + log(2) - x * tanh(x)
        end, Λ) / N
        # Mixing entropy of the two macro-sectors
        _xlogx(w) = w > 1e-300 ? w * log(w) : 0.0
        S_mix = -(_xlogx(w_AP) + _xlogx(w_P)) / N
        return w_AP * _s(Λ_AP) + w_P * _s(Λ_P) + S_mix

    elseif quantity === :specific_heat
        # Cv/N = β² [ w_AP * (Cv_AP/N + β²(E_AP/N - Ē)²)
        #           + w_P  * (Cv_P/N  + β²(E_P/N  - Ē)²) ]
        # where Ē = w_AP E_AP/N + w_P E_P/N
        _e(Λ) = -sum(λ -> (λ / 2) * tanh(β * λ / 2), Λ) / N
        _cv(Λ) = sum(λ -> (β * λ / 2)^2 * sech(β * λ / 2)^2, Λ) / N
        E_AP, E_P = _e(Λ_AP), _e(Λ_P)
        C_AP, C_P = _cv(Λ_AP), _cv(Λ_P)
        Ē = w_AP * E_AP + w_P * E_P
        return w_AP * (C_AP + β^2 * N * (E_AP - Ē)^2) + w_P * (C_P + β^2 * N * (E_P - Ē)^2)

    elseif quantity === :transverse_magnetization
        # ⟨σᶻ⟩/N = -(1/(Nβ)) ∂logZ/∂h
        #         = -(1/N) Σ_k (∂Λ/∂h) tanh(βΛ/2)
        # ∂Λ/∂h = (h − (Jx+Jy)cos k) / (Λ/2)  [using our dispersion with factor 2]
        #       Actually: Λ = 2√(A² + C²), A = h-(Jx+Jy)cos k, C = (Jx-Jy)sin k
        #       ∂Λ/∂h = 2A / Λ * (A/|A|) ... let me be careful:
        #       ∂Λ/∂h = 2 * (h - (Jx+Jy)cos k) / (Λ/2) ... no:
        #       Λ = 2√(A²+C²), ∂Λ/∂h = 2 * A / √(A²+C²) = 2A / (Λ/2) = 4A/Λ
        _mz_sector(Λ, ks) = sum(zip(Λ, ks)) do (λ, k)
            A = h - (Jx + Jy) * cos(k)
            (4.0 * A / λ) * tanh(β * λ / 2.0)
        end / N
        # The sign: ⟨σᶻ⟩ = (1/(Nβ)) ∂logZ/∂h = (1/N) Σ_k (∂logZ_k/∂(βh)) / β * β
        #   = -(1/N) Σ_k (∂Λ_k/∂h) * tanh(βΛ_k/2) * (-1) ... need to be careful.
        # Free energy f = -logZ/(Nβ),  ∂f/∂h = -1/(Nβ) ∂logZ/∂h
        # ⟨σᶻ⟩ = -∂f/∂h ... depends on sign convention.
        # Using our convention: H = - h Σ σᶻ, so ⟨σᶻ⟩ = -∂F/∂(Nh) = (1/Nβ)∂logZ/∂h
        # ∂logZ_s/∂h = Σ_k tanh(βΛ_k/2) * ∂(βΛ_k)/∂h / 2 * 2
        #            = Σ_k tanh(βΛ_k/2) * β * ∂Λ_k/∂h
        # ∂Λ_k/∂h = 4A/Λ where A = h-(Jx+Jy)cos k
        # Therefore ⟨σᶻ⟩_s = (1/N) Σ_k tanh(βΛ/2) * ∂Λ/∂h
        #                    = (1/N) Σ_k tanh(βΛ/2) * 4A/Λ
        # Note: this matches what's used in _xyh1d_thermo_infinite for :transverse_magnetization
        #   integrand = (2A/Λ) * tanh(βΛ/2) which integrates over [0,π] (half BZ), so factor of 2.
        # For the finite sum, we sum over the full k-range [0,2π), which gives the same.
        mz_AP = _mz_sector(Λ_AP, ks_AP)
        mz_P = _mz_sector(Λ_P, ks_P)
        # Weighted average from both sectors:
        # ∂logZ/∂h = w_AP * ∂logZ_AP/∂h + w_P * ∂logZ_P/∂h  (in terms of N*β factor)
        # ⟨σᶻ⟩/N = w_AP * ⟨σᶻ⟩_AP/N + w_P * ⟨σᶻ⟩_P/N
        return w_AP * mz_AP + w_P * mz_P

    elseif quantity === :transverse_susceptibility
        # χ = β * (⟨(σᶻ_total)²⟩ - ⟨σᶻ_total⟩²) / N
        # For each sector, the susceptibility within the sector:
        # χ_s = β/N Σ_{i,j} ⟨σᶻᵢ σᶻⱼ⟩_s,c
        # For the free-fermion sector, this is:
        # χ_s/N = (1/N) Σ_k [β sech²(βΛ_k/2) * (∂Λ_k/∂h)²/4
        #                   + (1/N) cross terms]
        # For translational invariance in PBC, the connected correlator is diagonal in k:
        # χ_s/N = (1/N) Σ_k sech²(βΛ_k/2) * β * (4A/Λ)²/4
        #       = (1/N) Σ_k sech²(βΛ_k/2) * β * (2A/Λ)²  ... let me re-derive.
        # ∂²logZ_s/∂(βh)² = Σ_k sech²(βΛ_k/2) * (β ∂Λ_k/∂h)²/4 + ... 
        # Actually: ∂²logZ_s/∂h² = β² Var_s(H_field) where H_field = -h Σ σᶻ
        # χ = -(1/N) ∂²F/∂h² = (1/(Nβ)) ∂²logZ/∂h²
        # ∂²logZ_s/∂h² = β Σ_k sech²(βΛ_k/2) * (∂Λ_k/∂h)² 
        #               + β Σ_k tanh(βΛ_k/2) * ∂²Λ_k/∂h²
        # ∂²Λ_k/∂h² = 4/Λ_k - 16A²/Λ_k³
        # Giving: χ_s = (1/N) Σ_k [ sech²(βΛ/2)*(4A/Λ)²*β/4 + tanh(βΛ/2)*(4/Λ - 16A²/Λ³) ]
        # Wait, let's be more careful. χ = (β/N) (⟨M²⟩ - ⟨M⟩²)
        # For PBC with translational invariance:
        # Use finite-difference approximation numerically (most robust):
        δh = h * 1e-5 + 1e-8
        mz_plus = _xyh1d_thermo_pbc(:transverse_magnetization, N, Jx, Jy, h + δh, β)
        mz_minus = _xyh1d_thermo_pbc(:transverse_magnetization, N, Jx, Jy, h - δh, β)
        return (mz_plus - mz_minus) / (2δh)

    else
        error("Unknown PBC thermal quantity: $quantity")
    end
end

"""
    _xyh1d_pbc_local_mz(N, Jx, Jy, h, β) -> Vector{Float64}

Site-local ⟨σᶻᵢ⟩ for PBC.  By translational invariance, all sites are equivalent;
returns a uniform vector of length N equal to the bulk ⟨σᶻ⟩/site.
"""
function _xyh1d_pbc_local_mz(N::Int, Jx::Real, Jy::Real, h::Real, β::Real)
    mz_bulk = _xyh1d_thermo_pbc(:transverse_magnetization, N, Jx, Jy, h, β)
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
