# ─────────────────────────────────────────────────────────────────────────────
# Transverse Field Ising Model — site-local thermal observables (OBC)
#
# Site-resolved thermal expectation values returned as a Vector of length N
# (not averaged over sites).  Exist so that bulk self-averaging (∑ᵢ / N) can
# be bypassed when comparing exact baselines against random-sampling / TPQ
# estimators, which themselves produce per-site data.
#
# Exposed quantities (OBC only):
#   :magnetization_x_local — ⟨σˣ_i⟩_β                   Vector{Float64} (N)
#   :magnetization_z_local — ⟨σᶻ_i⟩_β ≡ 0  by Z₂       Vector{Float64} (N)
#   :energy_local          — local energy density ε_i   Vector{Float64} (N)
#                            bonds split symmetrically so that
#                            Σᵢ ε_i = ⟨H⟩ exactly
#
# The N×N equal-time correlator ⟨σᶻ_i σᶻ_j⟩_β is provided by the existing
# `:zz_static_thermal` quantity in TFIM_dynamics.jl.
# ─────────────────────────────────────────────────────────────────────────────

# Nearest-neighbour ZZ bond expectation from the Majorana covariance.
#   σᶻ_i σᶻ_{i+1} = -i γ_{2i} γ_{2i+1}
#   ⟨σᶻ_i σᶻ_{i+1}⟩_β = -i · (i Σ[2i, 2i+1]) = Σ[2i, 2i+1]
_zz_bond(Σ::AbstractMatrix, i::Int) = Σ[2i, 2i + 1]

"""
    fetch(model::TFIM, ::MagnetizationXLocal, bc::OBC; beta::Float64, kwargs...)
        -> Vector{Float64}

Site-resolved transverse magnetisation `[⟨σˣ_i⟩_β for i = 1:N]` of the OBC
TFIM at inverse temperature `beta`, read off from the Majorana thermal
covariance as `Σ[2i-1, 2i]`.  `N` is taken from `bc.N`.

`beta = Inf` falls back to the ground state.
"""
function fetch(model::TFIM, ::MagnetizationXLocal, bc::OBC; beta::Float64, kwargs...)
    N = _bc_size(bc, kwargs)
    hmat = _majorana_ham(N, model.J, model.h)
    Σ = _majorana_thermal_covariance(hmat, beta)
    return Float64[_sx_expect(Σ, i) for i in 1:N]
end

"""
    fetch(model::TFIM, ::MagnetizationZLocal, bc::OBC; beta::Float64, kwargs...)
        -> Vector{Float64}

Site-resolved longitudinal magnetisation `[⟨σᶻ_i⟩_β for i = 1:N]`.
Identically zero in the OBC TFIM by the Z₂ symmetry σᶻ → −σᶻ of the
Hamiltonian (Gaussian state, odd product of Majoranas).  Returned as an
explicit zero vector so consumers can use it as an exact baseline against
finite random-sample estimates that fluctuate around zero.
"""
function fetch(::TFIM, ::MagnetizationZLocal, bc::OBC; beta::Float64, kwargs...)
    N = _bc_size(bc, kwargs)
    return zeros(Float64, N)
end

"""
    fetch(model::TFIM, ::EnergyLocal, bc::OBC; beta::Float64, kwargs...)
        -> Vector{Float64}

Site-local energy density `ε_i` of the OBC TFIM at inverse temperature
`beta`, defined so that `Σᵢ ε_i = ⟨H⟩_β`.  Each bond is split symmetrically
between its two endpoints:

    ε_i = -(J/2) (⟨σᶻ_{i-1} σᶻ_i⟩_β + ⟨σᶻ_i σᶻ_{i+1}⟩_β) - h ⟨σˣ_i⟩_β

with the missing bonds at the i = 1 and i = N boundaries taken to be zero.
Bond expectations are read off as `Σ(β)[2i, 2i+1]` from the Majorana thermal
covariance (exact, O(N) after the single 2N×2N diagonalisation).
"""
function fetch(model::TFIM, ::EnergyLocal, bc::OBC; beta::Float64, kwargs...)
    N = _bc_size(bc, kwargs)
    hmat = _majorana_ham(N, model.J, model.h)
    Σ = _majorana_thermal_covariance(hmat, beta)

    bonds = Float64[_zz_bond(Σ, i) for i in 1:(N - 1)]

    ε = Vector{Float64}(undef, N)
    @inbounds for i in 1:N
        left = i > 1 ? bonds[i - 1] : 0.0
        right = i < N ? bonds[i] : 0.0
        ε[i] = -(model.J / 2) * (left + right) - model.h * _sx_expect(Σ, i)
    end
    return ε
end
