# ─────────────────────────────────────────────────────────────────────────────
# Spin-1 transverse-field Ising chain (S1TFIM) — small-N dense ED reference.
#
# Hamiltonian (OBC):
#
#   H = -J Σᵢ Sᶻᵢ Sᶻᵢ₊₁  -  h Σᵢ Sˣᵢ,    spin = 1.
#
# `J > 0` ferromagnetic, `h` transverse field.  Local Hilbert space is
# 3-dimensional with Sᵅ matrices that have eigenvalues `±1, 0` (spin-1
# convention, *not* Pauli — see `_S1_x`, `_S1_y`, `_S1_z` defined in
# `HeisenbergS1.jl`).  This intentionally matches `S1Heisenberg1D`'s
# convention so cross-model identities (e.g. TFIM + ZZ-correlator
# verification) line up without rescaling.
#
# Unlike the spin-1/2 TFIM which factors via Jordan-Wigner +
# Bogoliubov-de Gennes (free fermions; see `TFIM.jl`), the spin-1 case
# is *not* a free-fermion theory: a single site already lives in a
# 3-state space and `(Sˣ)²` is not the identity, breaking the JW path.
# We therefore expose a small-N OBC dense-ED reference (3ᴺ Hilbert space,
# capped at N ≤ 8 by `_MAX_ED_SITES_S1`) intended for ThermalMPS-style
# benchmarks.  Note the spin-1 TFIM exhibits a non-trivial symmetry-
# breaking transition at finite `h_c/J` (numerical, no closed form);
# the tag here returns finite-N exact thermal observables only.
#
# References:
#   M. Suzuki, Prog. Theor. Phys. 56, 1454 (1976) — generalised quantum
#       Ising / spin-S transverse-field models.
#   F. C. Alcaraz and A. L. Malvezzi, J. Phys. A 28, 1521 (1995) — spin-1
#       quantum Ising chain phase diagram.
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra: I, Diagonal, Hermitian, eigen, eigvals, kron, tr

"""
    S1TFIM(; J::Real = 1.0, h::Real = 1.0) <: AbstractQAtlasModel

Spin-1 transverse-field Ising chain,

    H = -J Σᵢ Sᶻᵢ Sᶻᵢ₊₁ - h Σᵢ Sˣᵢ,    spin = 1.

Distinct from the spin-1/2 [`TFIM`](@ref): the local dimension is 3,
the model is *not* solved by Jordan-Wigner, and the critical field
`h_c` differs from the spin-1/2 self-dual point `h = J`.  Dense-ED
reference path only (N ≤ 8 by `_MAX_ED_SITES_S1`).
"""
struct S1TFIM <: AbstractQAtlasModel
    J::Float64
    h::Float64
end
S1TFIM(; J::Real=1.0, h::Real=1.0) = S1TFIM(Float64(J), Float64(h))

"""
    _s1_tfim_hamiltonian_matrix(model::S1TFIM, N::Int) -> Matrix{ComplexF64}

Assemble the `3^N × 3^N` OBC Hamiltonian

    H = -J Σᵢ Sᶻᵢ Sᶻᵢ₊₁ - h Σᵢ Sˣᵢ

via explicit tensor products built from the spin-1 primitives
`_S1_x`, `_S1_z` defined in `HeisenbergS1.jl`.  Capped by
`_MAX_ED_SITES_S1`.
"""
function _s1_tfim_hamiltonian_matrix(model::S1TFIM, N::Int)
    N ≥ 2 || throw(ArgumentError("S1TFIM OBC chain needs N ≥ 2 (got N = $N)"))
    N ≤ _MAX_ED_SITES_S1 || throw(
        ArgumentError("spin-1 dense ED is capped at N ≤ $(_MAX_ED_SITES_S1) (got N = $N)"),
    )
    J = model.J
    h = model.h
    D = 3^N
    H = zeros(ComplexF64, D, D)
    bond = (-J) * kron(_S1_z, _S1_z)
    for i in 1:(N - 1)
        d_left = 3^(i - 1)
        d_right = 3^(N - i - 1)
        H .+= kron(
            Matrix{ComplexF64}(I, d_left, d_left),
            bond,
            Matrix{ComplexF64}(I, d_right, d_right),
        )
    end
    for i in 1:N
        H .+= (-h) * _spin1_string(N, i => _S1_x)
    end
    return H
end

"""
    _s1_tfim_thermal_kernel(model, N, beta) -> NamedTuple

Mirror of `_s1_thermal_kernel` for `S1TFIM`.  Returns the dense
Hamiltonian, its eigendecomposition, normalised Boltzmann weights, and
the thermal density matrix in the product basis.
"""
function _s1_tfim_thermal_kernel(model::S1TFIM, N::Int, beta::Real)
    H = _s1_tfim_hamiltonian_matrix(model, N)
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

native_energy_granularity(::S1TFIM, ::OBC) = :total

"""
    fetch(model::S1TFIM, ::Energy{:total}, ::OBC; beta) -> Float64

Total thermal energy of the spin-1 OBC TFIM at finite N ≤ 8 via dense ED.
"""
function fetch(model::S1TFIM, ::Energy{:total}, bc::OBC; beta::Real, kwargs...)
    H = _s1_tfim_hamiltonian_matrix(model, bc.N)
    return _ed_thermal_energy(H, beta)
end

"""
    fetch(model::S1TFIM, ::FreeEnergy, ::OBC; beta) -> Float64

Per-site Helmholtz free energy `f(β) = -log Z / (Nβ)` for the spin-1 OBC TFIM.
"""
function fetch(model::S1TFIM, ::FreeEnergy, bc::OBC; beta::Real, kwargs...)
    H = _s1_tfim_hamiltonian_matrix(model, bc.N)
    return _ed_thermal_free_energy(H, beta) / bc.N
end

"""
    fetch(model::S1TFIM, ::ThermalEntropy, ::OBC; beta) -> Float64

Per-site Gibbs entropy `s(β) = β·(ε - f)` for the spin-1 OBC TFIM.
"""
function fetch(model::S1TFIM, ::ThermalEntropy, bc::OBC; beta::Real, kwargs...)
    H = _s1_tfim_hamiltonian_matrix(model, bc.N)
    return _ed_thermal_entropy(H, beta) / bc.N
end

"""
    fetch(model::S1TFIM, ::SpecificHeat, ::OBC; beta) -> Float64

Per-site heat capacity `c(β) = β²·Var(H)/N` for the spin-1 OBC TFIM.
"""
function fetch(model::S1TFIM, ::SpecificHeat, bc::OBC; beta::Real, kwargs...)
    H = _s1_tfim_hamiltonian_matrix(model, bc.N)
    return _ed_thermal_specific_heat(H, beta) / bc.N
end

function _s1_tfim_total_mag(N::Int, S::Matrix{ComplexF64})
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
            fetch(model::S1TFIM, ::$($Q), ::OBC; beta) -> Float64

        Per-site bulk magnetisation `⟨Σᵢ S^$($axis_str)_i⟩_β / N` of the
        spin-1 OBC TFIM (eigenvalues ±1, 0 — spin-1 convention).
        """
        function fetch(model::S1TFIM, ::$Q, bc::OBC; beta::Real, kwargs...)
            N = bc.N
            kernel = _s1_tfim_thermal_kernel(model, N, beta)
            S = _S1_AXIS_MATS.$axis_sym
            M = _s1_tfim_total_mag(N, S)
            return real(tr(M * kernel.ρ)) / N
        end
    end
end

for (Q, axis_sym) in
    ((:SusceptibilityXX, :x), (:SusceptibilityYY, :y), (:SusceptibilityZZ, :z))
    axis_str = string(axis_sym)
    @eval begin
        """
            fetch(model::S1TFIM, ::$($Q), ::OBC; beta) -> Float64

        Per-site uniform susceptibility `χ(β) = β·Var(M)/N` of the spin-1 OBC TFIM.
        """
        function fetch(model::S1TFIM, ::$Q, bc::OBC; beta::Real, kwargs...)
            N = bc.N
            kernel = _s1_tfim_thermal_kernel(model, N, beta)
            S = _S1_AXIS_MATS.$axis_sym
            M = _s1_tfim_total_mag(N, S)
            M2 = M * M
            mean_M = real(tr(M * kernel.ρ))
            mean_M2 = real(tr(M2 * kernel.ρ))
            return beta * (mean_M2 - mean_M^2) / N
        end
    end
end

"""
    fetch(model::S1TFIM, ::ZZCorrelation{:static}, ::OBC; beta, i, j) -> Float64

Static thermal correlator `⟨Sᶻ_i Sᶻ_j⟩_β` of the spin-1 OBC TFIM.
"""
function fetch(
    model::S1TFIM, ::ZZCorrelation{:static}, bc::OBC; beta::Real, i::Int, j::Int, kwargs...
)
    N = bc.N
    (1 ≤ i ≤ N && 1 ≤ j ≤ N) || throw(
        ArgumentError(
            "ZZCorrelation{:static}: site indices must satisfy 1 ≤ i,j ≤ N (got i=$i, j=$j, N=$N)",
        ),
    )
    kernel = _s1_tfim_thermal_kernel(model, N, beta)
    O = if i == j
        _spin1_string(N, i => _S1_z * _S1_z)
    else
        _spin1_string(N, i => _S1_z, j => _S1_z)
    end
    return real(tr(O * kernel.ρ))
end

"""
    fetch(model::S1TFIM, ::ZZCorrelation{:connected}, ::OBC; beta, i, j) -> Float64

Connected correlator `⟨Sᶻ_i Sᶻ_j⟩ - ⟨Sᶻ_i⟩·⟨Sᶻ_j⟩` for the spin-1 OBC TFIM.
"""
function fetch(
    model::S1TFIM,
    ::ZZCorrelation{:connected},
    bc::OBC;
    beta::Real,
    i::Int,
    j::Int,
    kwargs...,
)
    N = bc.N
    (1 ≤ i ≤ N && 1 ≤ j ≤ N) || throw(
        ArgumentError(
            "ZZCorrelation{:connected}: site indices must satisfy 1 ≤ i,j ≤ N (got i=$i, j=$j, N=$N)",
        ),
    )
    kernel = _s1_tfim_thermal_kernel(model, N, beta)
    O = if i == j
        _spin1_string(N, i => _S1_z * _S1_z)
    else
        _spin1_string(N, i => _S1_z, j => _S1_z)
    end
    Si = _spin1_string(N, i => _S1_z)
    Sj = i == j ? Si : _spin1_string(N, j => _S1_z)
    mean_SiSj = real(tr(O * kernel.ρ))
    mean_Si = real(tr(Si * kernel.ρ))
    mean_Sj = i == j ? mean_Si : real(tr(Sj * kernel.ρ))
    return mean_SiSj - mean_Si * mean_Sj
end

"""
    fetch(model::S1TFIM, ::MassGap, ::OBC) -> Float64

Many-body gap `Δ = E₁ - E₀` (first-excitation gap) of the spin-1 OBC
TFIM at finite N ≤ 8.  Computed from dense ED of the 3^N Hamiltonian.

Not a single-particle gap: the spin-1 TFIM is *not* a free-fermion theory
(no JW factorisation, since (Sˣ)² ≠ I for spin-1); this is the genuine
many-body gap from the lowest excited eigenvalue minus the ground-state
eigenvalue of H.
"""
function fetch(model::S1TFIM, ::MassGap, bc::OBC; kwargs...)
    H = _s1_tfim_hamiltonian_matrix(model, bc.N)
    evals = eigvals(Hermitian(H))
    return evals[2] - evals[1]
end

"""
    fetch(model::S1TFIM, ::ExactSpectrum, ::OBC) -> Vector{Float64}

Sorted real eigenvalues of the OBC spin-1 TFIM Hamiltonian (3^N × 3^N, N ≤ 8).
"""
function fetch(model::S1TFIM, ::ExactSpectrum, bc::OBC; kwargs...)
    H = _s1_tfim_hamiltonian_matrix(model, bc.N)
    return sort!(real.(eigvals(Hermitian(H))))
end
