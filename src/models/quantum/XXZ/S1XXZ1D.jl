# ─────────────────────────────────────────────────────────────────────────────
# Spin-1 XXZ chain (S1XXZ1D) — small-N dense ED reference.
#
# Hamiltonian (OBC):
#
#   H = J Σᵢ [ Sˣᵢ Sˣᵢ₊₁ + Sʸᵢ Sʸᵢ₊₁ + Δ Sᶻᵢ Sᶻᵢ₊₁ ],   spin = 1.
#
# Spin-1 operators with eigenvalues ±1, 0 (off-diagonal of Sˣ is 1/√2;
# see `_S1_x`, `_S1_y`, `_S1_z` in `HeisenbergS1.jl`).  At Δ = 1 this
# coincides with `S1Heisenberg1D` — the registered tests cross-check
# both paths.  For general Δ the spin-1 XXZ interpolates between the
# spin-1 XY model (Δ = 0) and the Haldane chain (Δ = 1); for large
# positive Δ it has Néel order, and Δ ≲ -1 enters a ferromagnetic
# phase (Schulz 1986).  No closed-form thermodynamic-limit formulae
# are exposed here — only finite-N exact dense ED suitable as a
# ThermalMPS reference (3^N Hilbert space, capped at N ≤ 8).
#
# Reference:
#   H. J. Schulz, Phys. Rev. B 34, 6372 (1986) — phase diagram of the
#       spin-1 XXZ chain.
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra: I, Diagonal, Hermitian, eigen, eigvals, kron, tr

"""
    S1XXZ1D(; J::Real = 1.0, Δ::Real = 0.0) <: AbstractQAtlasModel

Spin-1 XXZ chain,

    H = J Σᵢ [ Sˣ Sˣ + Sʸ Sʸ + Δ Sᶻ Sᶻ ],   spin = 1.

`Δ = 1` reproduces `S1Heisenberg1D` (Haldane chain).  Distinct from the
spin-1/2 [`XXZ1D`](@ref): different local dimension (3 vs 2), different
phase diagram (Haldane phase at Δ ≈ 1 — a gapped, topologically
non-trivial phase absent in the spin-1/2 case).  Dense-ED reference
path only (N ≤ 8 by `_MAX_ED_SITES_S1`).
"""
struct S1XXZ1D <: AbstractQAtlasModel
    J::Float64
    Δ::Float64
end
S1XXZ1D(; J::Real=1.0, Δ::Real=0.0) = S1XXZ1D(Float64(J), Float64(Δ))

"""
    _s1_xxz_hamiltonian_matrix(model::S1XXZ1D, N::Int) -> Matrix{ComplexF64}

Assemble the `3^N × 3^N` OBC Hamiltonian via spin-1 primitives.
Capped by `_MAX_ED_SITES_S1`.
"""
function _s1_xxz_hamiltonian_matrix(model::S1XXZ1D, N::Int)
    N ≥ 2 || throw(ArgumentError("S1XXZ1D OBC chain needs N ≥ 2 (got N = $N)"))
    N ≤ _MAX_ED_SITES_S1 || throw(
        ArgumentError("spin-1 dense ED is capped at N ≤ $(_MAX_ED_SITES_S1) (got N = $N)"),
    )
    J = model.J
    Δ = model.Δ
    D = 3^N
    bond = J * (kron(_S1_x, _S1_x) + kron(_S1_y, _S1_y) + Δ * kron(_S1_z, _S1_z))
    H = zeros(ComplexF64, D, D)
    for i in 1:(N - 1)
        d_left = 3^(i - 1)
        d_right = 3^(N - i - 1)
        H .+= kron(
            Matrix{ComplexF64}(I, d_left, d_left),
            bond,
            Matrix{ComplexF64}(I, d_right, d_right),
        )
    end
    return H
end

"""
    _s1_xxz_thermal_kernel(model, N, beta) -> NamedTuple

Mirror of `_s1_thermal_kernel` for `S1XXZ1D`.
"""
function _s1_xxz_thermal_kernel(model::S1XXZ1D, N::Int, beta::Real)
    H = _s1_xxz_hamiltonian_matrix(model, N)
    F = eigen(Hermitian(H))
    evals = F.values
    evecs = F.vectors
    if isinf(beta) && beta > 0
        emin = evals[1]
        gs_mask = (evals .- emin) .≤ 1e-12
        ng = count(gs_mask)
        weights = zeros(Float64, length(evals))
        weights[gs_mask] .= 1.0 / ng
    else
        emin = minimum(evals)
        ws = exp.(-beta .* (evals .- emin))
        Z = sum(ws)
        weights = ws ./ Z
    end
    ρ = evecs * Diagonal(weights) * evecs'
    return (; H=H, evals=evals, evecs=evecs, weights=weights, ρ=ρ)
end

native_energy_granularity(::S1XXZ1D, ::OBC) = :total

"""
    fetch(model::S1XXZ1D, ::Energy{:total}, ::OBC; beta) -> Float64

Total thermal energy `⟨H⟩_β` of the spin-1 OBC XXZ chain at finite
N ≤ 8 via dense ED.
"""
function fetch(model::S1XXZ1D, ::Energy{:total}, bc::OBC; beta::Real, kwargs...)
    H = _s1_xxz_hamiltonian_matrix(model, bc.N)
    return _ed_thermal_energy(H, beta)
end

"""
    fetch(model::S1XXZ1D, ::FreeEnergy, ::OBC; beta) -> Float64
"""
function fetch(model::S1XXZ1D, ::FreeEnergy, bc::OBC; beta::Real, kwargs...)
    H = _s1_xxz_hamiltonian_matrix(model, bc.N)
    return _ed_thermal_free_energy(H, beta) / bc.N
end

"""
    fetch(model::S1XXZ1D, ::ThermalEntropy, ::OBC; beta) -> Float64
"""
function fetch(model::S1XXZ1D, ::ThermalEntropy, bc::OBC; beta::Real, kwargs...)
    H = _s1_xxz_hamiltonian_matrix(model, bc.N)
    return _ed_thermal_entropy(H, beta) / bc.N
end

"""
    fetch(model::S1XXZ1D, ::SpecificHeat, ::OBC; beta) -> Float64
"""
function fetch(model::S1XXZ1D, ::SpecificHeat, bc::OBC; beta::Real, kwargs...)
    H = _s1_xxz_hamiltonian_matrix(model, bc.N)
    return _ed_thermal_specific_heat(H, beta) / bc.N
end

function _s1_xxz_total_mag(N::Int, S::Matrix{ComplexF64})
    D = 3^N
    M = zeros(ComplexF64, D, D)
    @inbounds for i in 1:N
        M .+= _spin1_string(N, i => S)
    end
    return M
end

for (Q, axis_sym) in ((:MagnetizationX, :x), (:MagnetizationY, :y), (:MagnetizationZ, :z))
    axis_str = string(axis_sym)
    @eval begin
        """
            fetch(model::S1XXZ1D, ::$($Q), ::OBC; beta) -> Float64

        Per-site bulk magnetisation `⟨Σᵢ S^$($axis_str)_i⟩_β / N` of the
        spin-1 OBC XXZ chain.
        """
        function fetch(model::S1XXZ1D, ::$Q, bc::OBC; beta::Real, kwargs...)
            N = bc.N
            kernel = _s1_xxz_thermal_kernel(model, N, beta)
            S = _S1_AXIS_MATS.$axis_sym
            M = _s1_xxz_total_mag(N, S)
            return real(tr(M * kernel.ρ)) / N
        end
    end
end

for (Q, axis_sym) in
    ((:SusceptibilityXX, :x), (:SusceptibilityYY, :y), (:SusceptibilityZZ, :z))
    axis_str = string(axis_sym)
    @eval begin
        """
            fetch(model::S1XXZ1D, ::$($Q), ::OBC; beta) -> Float64

        Per-site uniform susceptibility `χ(β) = β·Var(M)/N`.
        """
        function fetch(model::S1XXZ1D, ::$Q, bc::OBC; beta::Real, kwargs...)
            N = bc.N
            kernel = _s1_xxz_thermal_kernel(model, N, beta)
            S = _S1_AXIS_MATS.$axis_sym
            M = _s1_xxz_total_mag(N, S)
            M2 = M * M
            mean_M = real(tr(M * kernel.ρ))
            mean_M2 = real(tr(M2 * kernel.ρ))
            return beta * (mean_M2 - mean_M^2) / N
        end
    end
end

for (Q, axis_sym) in ((:XXCorrelation, :x), (:YYCorrelation, :y), (:ZZCorrelation, :z))
    axis_str = string(axis_sym)
    @eval begin
        """
            fetch(model::S1XXZ1D, ::$($Q){:static}, ::OBC; beta, i, j) -> Float64

        Static thermal correlator `⟨S^$($axis_str)_i S^$($axis_str)_j⟩_β`
        for the spin-1 OBC XXZ chain.
        """
        function fetch(
            model::S1XXZ1D, ::$Q{:static}, bc::OBC; beta::Real, i::Int, j::Int, kwargs...
        )
            N = bc.N
            (1 ≤ i ≤ N && 1 ≤ j ≤ N) || throw(
                ArgumentError(
                    "$($Q){:static}: site indices must satisfy 1 ≤ i,j ≤ N (got i=$i, j=$j, N=$N)",
                ),
            )
            kernel = _s1_xxz_thermal_kernel(model, N, beta)
            S = _S1_AXIS_MATS.$axis_sym
            O = i == j ? _spin1_string(N, i => S * S) : _spin1_string(N, i => S, j => S)
            return real(tr(O * kernel.ρ))
        end

        """
            fetch(model::S1XXZ1D, ::$($Q){:connected}, ::OBC; beta, i, j) -> Float64

        Connected correlator for the spin-1 OBC XXZ chain.
        """
        function fetch(
            model::S1XXZ1D, ::$Q{:connected}, bc::OBC; beta::Real, i::Int, j::Int, kwargs...
        )
            N = bc.N
            (1 ≤ i ≤ N && 1 ≤ j ≤ N) || throw(
                ArgumentError(
                    "$($Q){:connected}: site indices must satisfy 1 ≤ i,j ≤ N (got i=$i, j=$j, N=$N)",
                ),
            )
            kernel = _s1_xxz_thermal_kernel(model, N, beta)
            S = _S1_AXIS_MATS.$axis_sym
            O = i == j ? _spin1_string(N, i => S * S) : _spin1_string(N, i => S, j => S)
            Si = _spin1_string(N, i => S)
            Sj = i == j ? Si : _spin1_string(N, j => S)
            mean_SiSj = real(tr(O * kernel.ρ))
            mean_Si = real(tr(Si * kernel.ρ))
            mean_Sj = i == j ? mean_Si : real(tr(Sj * kernel.ρ))
            return mean_SiSj - mean_Si * mean_Sj
        end
    end
end

"""
    fetch(model::S1XXZ1D, ::MassGap, ::OBC) -> Float64

Single-particle gap `Δ = E₁ - E₀` of the spin-1 OBC XXZ chain at finite N ≤ 8.
"""
function fetch(model::S1XXZ1D, ::MassGap, bc::OBC; kwargs...)
    H = _s1_xxz_hamiltonian_matrix(model, bc.N)
    evals = eigvals(Hermitian(H))
    return evals[2] - evals[1]
end

"""
    fetch(model::S1XXZ1D, ::ExactSpectrum, ::OBC) -> Vector{Float64}

Sorted real eigenvalues of the OBC spin-1 XXZ Hamiltonian (3^N × 3^N, N ≤ 8).
"""
function fetch(model::S1XXZ1D, ::ExactSpectrum, bc::OBC; kwargs...)
    H = _s1_xxz_hamiltonian_matrix(model, bc.N)
    return sort!(real.(eigvals(Hermitian(H))))
end
