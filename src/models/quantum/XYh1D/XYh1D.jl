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

"""
    fetch(model::XYh1D, ::FreeEnergy, ::Infinite; beta) -> Float64

Per-site free-energy density of the infinite XYh1D chain at inverse
temperature β, via Gauss-Kronrod integration over the Bogoliubov
dispersion (Lieb-Schultz-Mattis 1961).
"""
function fetch(m::XYh1D, ::FreeEnergy, ::Infinite; beta::Real, kwargs...)
    return _xyh1d_thermo_infinite(:free_energy, m.Jx, m.Jy, m.h, beta)
end

"""
    fetch(model::XYh1D, ::ThermalEntropy, ::Infinite; beta) -> Float64
"""
function fetch(m::XYh1D, ::ThermalEntropy, ::Infinite; beta::Real, kwargs...)
    return _xyh1d_thermo_infinite(:entropy, m.Jx, m.Jy, m.h, beta)
end

"""
    fetch(model::XYh1D, ::SpecificHeat, ::Infinite; beta) -> Float64
"""
function fetch(m::XYh1D, ::SpecificHeat, ::Infinite; beta::Real, kwargs...)
    return _xyh1d_thermo_infinite(:specific_heat, m.Jx, m.Jy, m.h, beta)
end

"""
    fetch(model::XYh1D, ::FreeEnergy, bc::OBC; beta) -> Float64
"""
function fetch(m::XYh1D, ::FreeEnergy, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    return _xyh1d_thermo_obc(:free_energy, N, m.Jx, m.Jy, m.h, beta)
end

"""
    fetch(model::XYh1D, ::ThermalEntropy, bc::OBC; beta) -> Float64
"""
function fetch(m::XYh1D, ::ThermalEntropy, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    return _xyh1d_thermo_obc(:entropy, N, m.Jx, m.Jy, m.h, beta)
end

"""
    fetch(model::XYh1D, ::SpecificHeat, bc::OBC; beta) -> Float64
"""
function fetch(m::XYh1D, ::SpecificHeat, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    return _xyh1d_thermo_obc(:specific_heat, N, m.Jx, m.Jy, m.h, beta)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Magnetization & Susceptibility dispatchers (Phase 2, #292)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::XYh1D, ::MagnetizationZ, ::Infinite; beta) -> Float64
"""
function fetch(m::XYh1D, ::MagnetizationZ, ::Infinite; beta::Real, kwargs...)
    return _xyh1d_thermo_infinite(:transverse_magnetization, m.Jx, m.Jy, m.h, beta)
end

"""
    fetch(model::XYh1D, ::MagnetizationZ, bc::OBC; beta) -> Float64
"""
function fetch(m::XYh1D, ::MagnetizationZ, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    return _xyh1d_thermo_obc(:transverse_magnetization, N, m.Jx, m.Jy, m.h, beta)
end

"""
    fetch(model::XYh1D, ::SusceptibilityZZ, ::Infinite; beta) -> Float64
"""
function fetch(m::XYh1D, ::SusceptibilityZZ, ::Infinite; beta::Real, kwargs...)
    return _xyh1d_thermo_infinite(:transverse_susceptibility, m.Jx, m.Jy, m.h, beta)
end

"""
    fetch(model::XYh1D, ::SusceptibilityZZ, bc::OBC; beta) -> Float64
"""
function fetch(m::XYh1D, ::SusceptibilityZZ, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    return _xyh1d_thermo_obc(:transverse_susceptibility, N, m.Jx, m.Jy, m.h, beta)
end

"""
    fetch(model::XYh1D, ::MagnetizationZLocal, bc::OBC; beta) -> Vector{Float64}
"""
function fetch(m::XYh1D, ::MagnetizationZLocal, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    hmat = _xyh1d_majorana_ham(N, m.Jx, m.Jy, m.h)
    Σ = _xyh1d_majorana_thermal_covariance(hmat, beta)
    return Float64[Σ[2i - 1, 2i] for i in 1:N]
end

# ═══════════════════════════════════════════════════════════════════════════════
# Site-local observables (Phase 2, #292)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::XYh1D, ::MagnetizationXLocal{:equilibrium}, bc::OBC; beta) -> Vector{Float64}

Site-resolved ⟨σˣ_i⟩ at equilibrium. By the residual Z₂ symmetry σˣ → −σˣ of
the XYh1D Hamiltonian at zero longitudinal field, the equilibrium expectation
vanishes identically: returns the zero vector of length N.
"""
function fetch(::XYh1D, ::MagnetizationXLocal{:equilibrium}, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    return zeros(Float64, N)
end

"""
    fetch(model::XYh1D, ::MagnetizationYLocal, bc::OBC; beta) -> Vector{Float64}

Site-resolved ⟨σʸ_i⟩. Vanishes by Z₂ symmetry σʸ → −σʸ; returns zeros.
"""
function fetch(::XYh1D, ::MagnetizationYLocal, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    return zeros(Float64, N)
end

"""
    fetch(model::XYh1D, ::EnergyLocal, bc::OBC; beta) -> Vector{Float64}

Site-resolved energy density ε_i = ⟨H_i⟩ on the OBC chain, where H_i is the
symmetric bond split. Computed from the Majorana thermal covariance.
"""
function fetch(model::XYh1D, ::EnergyLocal, bc::OBC; beta::Real, kwargs...)
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
