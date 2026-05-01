# ─────────────────────────────────────────────────────────────────────────────
# TFIM — static σˣσˣ and σʸσʸ structure factors at OBC and Infinite.
#
# Defined as the Fourier transform of the static (equal-time) correlator
#
#     S_αα(q, β) = (1/N) Σ_{i,j} e^{-iq(i-j)} ⟨σᵅ_i σᵅ_j⟩_β.
#
# For OBC the lattice lacks translation invariance so this is the
# boundary-aware "naïve" double sum; in a long enough chain the bulk
# contribution dominates and the value converges to the translation-
# invariant result.  An `Infinite()` method routes through an OBC
# `N_proxy` calculation, mirroring the existing `ZZStructureFactor,
# Infinite` and `SusceptibilityZZ, Infinite` proxies in `TFIM_zaxis.jl`.
#
# σˣ_i / σʸ_i correlators come from the Pfaffian / Wick-contraction
# helpers in `TFIM_dynamics.jl` and `TFIM_yy.jl`.  Diagonal entries
# (`i = j`) use the Pauli identity `(σ^α)² = I` so `⟨(σ^α_i)²⟩ = 1`
# without a Pfaffian call.
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Static N × N correlator matrices (XX, YY) — builds the full C[i,j] from
# one Majorana diagonalisation, then a single Fourier sweep.
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _sx_sx_static_thermal(N, J, h, β) -> Matrix{Float64}

`C[i, j] = ⟨σˣ_i σˣ_j⟩_β` at OBC, evaluated as the t = 0 slice of the
free-fermion Pfaffian / Wick formula.  Symmetric, diagonal `= 1` exactly
(`σˣ² = I`).  One Majorana diagonalisation amortised over N(N-1)/2
small (4×4) Pfaffians.
"""
function _sx_sx_static_thermal(N::Int, J::Float64, h::Float64, β::Real)
    hmat = _majorana_ham(N, J, h)
    Σ = _majorana_thermal_covariance(hmat, β)
    R = _majorana_evolution(hmat, 0.0)
    C = zeros(Float64, N, N)
    @inbounds for i in 1:N
        C[i, i] = 1.0
        for j in (i + 1):N
            v = real(_sx_sx_corr_from_cached(Σ, R, i, j))
            C[i, j] = v
            C[j, i] = v
        end
    end
    return C
end

"""
    _sy_sy_static_thermal(N, J, h, β) -> Matrix{Float64}

`C[i, j] = ⟨σʸ_i σʸ_j⟩_β` at OBC.  Same structure as
[`_sx_sx_static_thermal`](@ref) but uses the σʸ Majorana index list
defined in `TFIM_yy.jl`.
"""
function _sy_sy_static_thermal(N::Int, J::Float64, h::Float64, β::Real)
    hmat = _majorana_ham(N, J, h)
    Σ = _majorana_thermal_covariance(hmat, β)
    R = _majorana_evolution(hmat, 0.0)
    C = zeros(Float64, N, N)
    @inbounds for i in 1:N
        C[i, i] = 1.0
        for j in (i + 1):N
            v = real(_sy_sy_corr_from_cached(Σ, R, i, j))
            C[i, j] = v
            C[j, i] = v
        end
    end
    return C
end

# ═══════════════════════════════════════════════════════════════════════════════
# Static structure factor (Fourier sum over the full C[i, j] matrix)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _xx_static_structure_factor(N, J, h, β, q) -> Float64

`S_xx(q, β) = (1/N) Σ_{i,j} e^{-iq(i-j)} ⟨σˣ_i σˣ_j⟩_β` for OBC.
Direct double sum over the N×N matrix returned by
[`_sx_sx_static_thermal`](@ref).
"""
function _xx_static_structure_factor(N::Int, J::Float64, h::Float64, β::Real, q::Real)
    C = _sx_sx_static_thermal(N, J, h, β)
    s = 0.0 + 0.0im
    @inbounds for i in 1:N, j in 1:N
        s += exp(-im * q * (i - j)) * C[i, j]
    end
    return real(s) / N
end

"""
    _yy_static_structure_factor(N, J, h, β, q) -> Float64

`S_yy(q, β) = (1/N) Σ_{i,j} e^{-iq(i-j)} ⟨σʸ_i σʸ_j⟩_β` for OBC.
"""
function _yy_static_structure_factor(N::Int, J::Float64, h::Float64, β::Real, q::Real)
    C = _sy_sy_static_thermal(N, J, h, β)
    s = 0.0 + 0.0im
    @inbounds for i in 1:N, j in 1:N
        s += exp(-im * q * (i - j)) * C[i, j]
    end
    return real(s) / N
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch — XXStructureFactor / YYStructureFactor at OBC and Infinite
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::TFIM, ::XXStructureFactor, bc::OBC; beta::Real, q::Real, kwargs...)
        -> Float64

Static transverse structure factor `S_xx(q, β)` for the OBC TFIM with
N sites.  Defined as `(1/N) Σ_{i,j} e^{-iq(i-j)} ⟨σˣ_i σˣ_j⟩_β` with σˣ
correlators from the t = 0 slice of the free-fermion Pfaffian formula.
"""
function fetch(model::TFIM, ::XXStructureFactor, bc::OBC; beta::Real, q::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    return _xx_static_structure_factor(N, model.J, model.h, beta, q)
end

"""
    fetch(model::TFIM, ::YYStructureFactor, bc::OBC; beta::Real, q::Real, kwargs...)
        -> Float64

Static σʸ structure factor for the OBC TFIM.  Companion of
[`XXStructureFactor`](@ref).
"""
function fetch(model::TFIM, ::YYStructureFactor, bc::OBC; beta::Real, q::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    return _yy_static_structure_factor(N, model.J, model.h, beta, q)
end

"""
    fetch(model::TFIM, ::XXStructureFactor, ::Infinite;
          beta::Real, q::Real, N_proxy::Int = 80, kwargs...) -> Float64

Static transverse structure factor `S_xx(q, β)` in the thermodynamic
limit, computed as the OBC large-N proxy at `N_proxy = 80` (default,
~3-digit accuracy at moderate β in the gapped phase; raise `N_proxy`
to tighten).  Same convention and proxy strategy as the existing
[`SusceptibilityZZ`](@ref) / [`ZZStructureFactor`](@ref) Infinite
methods in `TFIM_zaxis.jl`.
"""
function fetch(
    model::TFIM,
    ::XXStructureFactor,
    ::Infinite;
    beta::Real,
    q::Real,
    N_proxy::Int=80,
    kwargs...,
)
    return _xx_static_structure_factor(N_proxy, model.J, model.h, beta, q)
end

"""
    fetch(model::TFIM, ::YYStructureFactor, ::Infinite;
          beta::Real, q::Real, N_proxy::Int = 80, kwargs...) -> Float64

Static σʸ structure factor in the thermodynamic limit; OBC large-N
proxy at `N_proxy`.
"""
function fetch(
    model::TFIM,
    ::YYStructureFactor,
    ::Infinite;
    beta::Real,
    q::Real,
    N_proxy::Int=80,
    kwargs...,
)
    return _yy_static_structure_factor(N_proxy, model.J, model.h, beta, q)
end
