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
