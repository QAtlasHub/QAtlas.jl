# ─────────────────────────────────────────────────────────────────────────────
# Spin-1 Heisenberg chain — extended observable coverage
#
# Companion to `HeisenbergS1.jl`.  Adds the OBC observable suite that
# matches the granularity of the TFIM coverage (per-site magnetisations,
# variance susceptibilities, two-point correlators with `:static` and
# `:connected` modes, site-resolved local quantities, the gap and
# entanglement entropies).  Every routine reuses the single
# eigendecomposition produced by `_s1_thermal_kernel` to keep multi-
# observable callers cheap.
#
# Convention for spin-1 operators (see also `_S1_x`, `_S1_y`, `_S1_z`):
# the matrices `Sᵅ` carry physical spin-1 eigenvalues `±1, 0`, i.e. the
# off-diagonals of `Sˣ` are `1/√2`.  This is the "spin operator"
# convention, **not** the Pauli convention used by the TFIM σ-matrices
# (eigenvalues ±1, off-diagonal 1).  Magnetisation / susceptibility /
# correlator outputs of this file are therefore `⟨Sᵅ⟩`, `⟨SᵅSᵅ⟩`, …
# (not `⟨σᵅ⟩`).  Downstream comparisons against TFIM-style data should
# multiply by the appropriate factor of 2 / 4.
#
# Infinite-system observables exposed here are literature constants
# (Haldane gap, GS energy density of the AFM spin-1 chain) and tagged
# `reliability=:medium` in the registry so consumers can filter them
# distinctly from finite-N exact data.
#
# References:
#   F. D. M. Haldane, PRL 50, 1153 (1983).
#   S. R. White, D. A. Huse, Phys. Rev. B 48, 3844 (1993) —
#     DMRG values  Δ ≈ 0.41048 J,  e₀ ≈ -1.401484 J.
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Magnetisations  (per-site, OBC)
# ═══════════════════════════════════════════════════════════════════════════════

# Build the bulk magnetisation operator M_α = Σᵢ Sᵅᵢ once for a given (α, N).
function _s1_total_mag(N::Int, S::Matrix{ComplexF64})
    D = 3^N
    M = zeros(ComplexF64, D, D)
    @inbounds for i in 1:N
        M .+= _spin1_string(N, i => S)
    end
    return M
end

const _S1_AXIS_MATS = (x=_S1_x, y=_S1_y, z=_S1_z)

for (Q, axis_sym) in ((:MagnetizationX, :x), (:MagnetizationY, :y), (:MagnetizationZ, :z))
    axis_str = string(axis_sym)
    @eval begin
        """
            fetch(model::S1Heisenberg1D, ::$($Q), ::OBC; beta) -> Float64

        Per-site bulk magnetisation `⟨Σᵢ S^$($axis_str)_i⟩_β / N` of the spin-1
        Heisenberg OBC chain at finite `N ≤ $(_MAX_ED_SITES_S1)`, computed
        from the dense thermal density matrix.

        Spin-1 convention: `S^$($axis_str)` carries eigenvalues `±1, 0`
        (off-diagonal `1/√2` for `S^x`); see file header.  The SU(2)-symmetric
        AFM ground state has `⟨S^α⟩ = 0` in every direction at any
        temperature, so this is a non-trivial check only for symmetry-
        breaking finite-β samples (and for benchmarking sampler bias).
        """
        function fetch(model::S1Heisenberg1D, ::$Q, bc::OBC; beta::Real, kwargs...)
            N = bc.N
            kernel = _s1_thermal_kernel(model, N, beta)
            M = _s1_total_mag(N, _S1_AXIS_MATS.$axis_sym)
            return _s1_thermal_expectation(kernel, M) / N
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Susceptibilities  (per-site, OBC)
# ═══════════════════════════════════════════════════════════════════════════════
#
# χ_αα(β) = β · Var(M_α) / N,   M_α = Σᵢ Sᵅᵢ
#         = β · (⟨M_α²⟩ − ⟨M_α⟩²) / N

for (Q, axis_sym) in
    ((:SusceptibilityXX, :x), (:SusceptibilityYY, :y), (:SusceptibilityZZ, :z))
    axis_str = string(axis_sym)
    @eval begin
        """
            fetch(model::S1Heisenberg1D, ::$($Q), ::OBC; beta) -> Float64

        Per-site uniform $(uppercase($axis_str * $axis_str)) susceptibility
        `χ_$($axis_str * $axis_str)(β) = β · Var(M_$($axis_str)) / N` of the
        spin-1 Heisenberg OBC chain via the dense thermal density matrix.

        At infinite temperature each site contributes `Tr((Sᵅ)²)/3 = 2/3` so
        `χ_αα → 2β/3` for any axis.
        """
        function fetch(model::S1Heisenberg1D, ::$Q, bc::OBC; beta::Real, kwargs...)
            N = bc.N
            kernel = _s1_thermal_kernel(model, N, beta)
            M = _s1_total_mag(N, _S1_AXIS_MATS.$axis_sym)
            M2 = M * M
            mean_M = _s1_thermal_expectation(kernel, M)
            mean_M2 = _s1_thermal_expectation(kernel, M2)
            return beta * (mean_M2 - mean_M^2) / N
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Two-point correlators  (static / connected, OBC)
# ═══════════════════════════════════════════════════════════════════════════════

for (Q, axis_sym) in ((:XXCorrelation, :x), (:YYCorrelation, :y), (:ZZCorrelation, :z))
    axis_str = string(axis_sym)
    @eval begin
        """
            fetch(model::S1Heisenberg1D, ::$($Q){:static}, ::OBC;
                  beta, i::Int, j::Int) -> Float64

        Static thermal correlator `⟨S^$($axis_str)_i S^$($axis_str)_j⟩_β`
        of the spin-1 OBC Heisenberg chain.  Caller passes both site
        indices `1 ≤ i, j ≤ N`; equal sites give `⟨(S^α)²⟩_β`.
        """
        function fetch(
            model::S1Heisenberg1D,
            ::$Q{:static},
            bc::OBC;
            beta::Real,
            i::Int,
            j::Int,
            kwargs...,
        )
            N = bc.N
            (1 ≤ i ≤ N && 1 ≤ j ≤ N) || throw(
                ArgumentError(
                    "$($Q){:static}: site indices must satisfy 1 ≤ i,j ≤ N (got i=$i, j=$j, N=$N)",
                ),
            )
            kernel = _s1_thermal_kernel(model, N, beta)
            S = _S1_AXIS_MATS.$axis_sym
            O = i == j ? _spin1_string(N, i => S * S) : _spin1_string(N, i => S, j => S)
            return _s1_thermal_expectation(kernel, O)
        end

        """
            fetch(model::S1Heisenberg1D, ::$($Q){:connected}, ::OBC;
                  beta, i::Int, j::Int) -> Float64

        Connected (cumulant) correlator
        `⟨S^$($axis_str)_i S^$($axis_str)_j⟩ - ⟨S^$($axis_str)_i⟩·⟨S^$($axis_str)_j⟩`
        for the spin-1 OBC Heisenberg chain.
        """
        function fetch(
            model::S1Heisenberg1D,
            ::$Q{:connected},
            bc::OBC;
            beta::Real,
            i::Int,
            j::Int,
            kwargs...,
        )
            N = bc.N
            (1 ≤ i ≤ N && 1 ≤ j ≤ N) || throw(
                ArgumentError(
                    "$($Q){:connected}: site indices must satisfy 1 ≤ i,j ≤ N " *
                    "(got i=$i, j=$j, N=$N)",
                ),
            )
            kernel = _s1_thermal_kernel(model, N, beta)
            S = _S1_AXIS_MATS.$axis_sym
            O = i == j ? _spin1_string(N, i => S * S) : _spin1_string(N, i => S, j => S)
            Si = _spin1_string(N, i => S)
            Sj = i == j ? Si : _spin1_string(N, j => S)
            mean_SiSj = _s1_thermal_expectation(kernel, O)
            mean_Si = _s1_thermal_expectation(kernel, Si)
            mean_Sj = i == j ? mean_Si : _s1_thermal_expectation(kernel, Sj)
            return mean_SiSj - mean_Si * mean_Sj
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Local (site-resolved) observables  (OBC)
# ═══════════════════════════════════════════════════════════════════════════════

for (Q, axis_sym) in ((:MagnetizationXLocal, :x), (:MagnetizationZLocal, :z))
    axis_str = string(axis_sym)
    @eval begin
        """
            fetch(model::S1Heisenberg1D, ::$($Q), ::OBC; beta) -> Vector{Float64}

        Site-resolved `[⟨S^$($axis_str)_i⟩_β  for i = 1:N]` of the spin-1 OBC
        Heisenberg chain.  Sums to `N · MagnetizationX/Y/Z` (per-site bulk
        average) by construction.
        """
        function fetch(model::S1Heisenberg1D, ::$Q, bc::OBC; beta::Real, kwargs...)
            N = bc.N
            kernel = _s1_thermal_kernel(model, N, beta)
            S = _S1_AXIS_MATS.$axis_sym
            return Float64[
                _s1_thermal_expectation(kernel, _spin1_string(N, i => S)) for i in 1:N
            ]
        end
    end
end

"""
    fetch(model::S1Heisenberg1D, ::EnergyLocal, ::OBC; beta) -> Vector{Float64}

Site-resolved local energy density `ε_i` of the OBC spin-1 Heisenberg
chain such that `Σᵢ ε_i = ⟨H⟩_β`.  Each bond contribution
`b_i = J ⟨Sᵢ · Sᵢ₊₁⟩_β` is split symmetrically between its two
endpoints, with the missing left/right bond at `i = 1` / `i = N` set to
zero — i.e.

    ε_i = (1/2) (b_{i-1} + b_i),    b_0 ≡ b_N ≡ 0.
"""
function fetch(model::S1Heisenberg1D, ::EnergyLocal, bc::OBC; beta::Real, kwargs...)
    N = bc.N
    kernel = _s1_thermal_kernel(model, N, beta)
    J = model.J
    bonds = Vector{Float64}(undef, N - 1)
    @inbounds for i in 1:(N - 1)
        Bxx = _spin1_string(N, i => _S1_x, i + 1 => _S1_x)
        Byy = _spin1_string(N, i => _S1_y, i + 1 => _S1_y)
        Bzz = _spin1_string(N, i => _S1_z, i + 1 => _S1_z)
        bonds[i] =
            J * (
                _s1_thermal_expectation(kernel, Bxx) +
                _s1_thermal_expectation(kernel, Byy) +
                _s1_thermal_expectation(kernel, Bzz)
            )
    end
    ε = Vector{Float64}(undef, N)
    @inbounds for i in 1:N
        left = i > 1 ? bonds[i - 1] : 0.0
        right = i < N ? bonds[i] : 0.0
        ε[i] = (left + right) / 2
    end
    return ε
end

# ═══════════════════════════════════════════════════════════════════════════════
# MassGap (OBC) — energy difference between ground and first excited state
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::S1Heisenberg1D, ::MassGap, ::OBC) -> Float64

Single-particle gap of the spin-1 Heisenberg OBC chain at finite
`N ≤ $(_MAX_ED_SITES_S1)`,

    Δ = E₁ - E₀

between the lowest two eigenvalues of the dense Hamiltonian.  At fixed
`N` this contains finite-size and edge-state corrections relative to the
bulk Haldane gap `Δ_∞ ≈ 0.41048 J`.
"""
function fetch(model::S1Heisenberg1D, ::MassGap, bc::OBC; kwargs...)
    H = _s1_heisenberg_hamiltonian_matrix(model, bc.N)
    evals = eigvals(Hermitian(H))
    return evals[2] - evals[1]
end

# ═══════════════════════════════════════════════════════════════════════════════
# Entanglement entropies (OBC, ground state or thermal)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _s1_partial_trace_A(ρ::AbstractMatrix, ℓ::Int, N::Int) -> Matrix{ComplexF64}

Trace out sites `ℓ+1 .. N` of the `3^N × 3^N` density matrix `ρ` and
return the `3^ℓ × 3^ℓ` reduced density matrix on sites `1 .. ℓ` (local
dimension `d = 3` for spin-1).

The convention matches `_spin1_string`: site 1 is the leftmost (slowest-
running) tensor factor, so the row/column index decomposes as
`(i_A, i_B)` with `i_A = 0:dA-1, i_B = 0:dB-1` where `dA = d^ℓ,
dB = d^(N-ℓ)`.  Reshaping ρ as `(dA, dB, dA, dB)` and contracting the B
indices then gives `ρ_A`.
"""
function _s1_partial_trace_A(ρ::AbstractMatrix, ℓ::Int, N::Int)
    1 ≤ ℓ ≤ N - 1 ||
        throw(ArgumentError("partial trace: ℓ must satisfy 1 ≤ ℓ ≤ N-1 (got ℓ=$ℓ, N=$N)"))
    d = 3
    dA = d^ℓ
    dB = d^(N - ℓ)
    R = reshape(Array(ρ), (dB, dA, dB, dA))
    ρA = zeros(eltype(ρ), dA, dA)
    @inbounds for a1 in 1:dA, a2 in 1:dA
        s = zero(eltype(ρ))
        for b in 1:dB
            s += R[b, a1, b, a2]
        end
        ρA[a1, a2] = s
    end
    return ρA
end

"""
    fetch(model::S1Heisenberg1D, ::VonNeumannEntropy, ::OBC;
          ℓ::Int, beta::Real = Inf) -> Float64

Von Neumann entanglement entropy `S_vN = -Tr ρ_A log ρ_A` of the first
`ℓ` sites of the OBC spin-1 Heisenberg chain.  The reduced density
matrix `ρ_A` is the partial trace of the thermal density matrix
`ρ = exp(-βH) / Z` (or the GS projector when `beta = Inf`) over sites
`ℓ+1 .. N`.

Both `ℓ` and `N = bc.N` are bounded by `_MAX_ED_SITES_S1` because the
full `3^N × 3^N` ρ is built explicitly.
"""
function fetch(
    model::S1Heisenberg1D, ::VonNeumannEntropy, bc::OBC; ℓ::Int, beta::Real=Inf, kwargs...
)
    N = bc.N
    1 ≤ ℓ ≤ N - 1 || throw(
        ArgumentError("VonNeumannEntropy: ℓ must satisfy 1 ≤ ℓ ≤ N - 1 (got ℓ=$ℓ, N=$N)"),
    )
    kernel = _s1_thermal_kernel(model, N, beta)
    ρA = _s1_partial_trace_A(kernel.ρ, ℓ, N)
    λs = real.(eigvals(Hermitian((ρA + ρA') / 2)))
    S = 0.0
    @inbounds for λ in λs
        if λ > 1e-15
            S -= λ * log(λ)
        end
    end
    return S
end

"""
    fetch(model::S1Heisenberg1D, q::RenyiEntropy, ::OBC;
          ℓ::Int, beta::Real = Inf) -> Float64

Rényi entropy `S_α = log Tr ρ_A^α / (1 - α)` for the OBC spin-1
Heisenberg chain.  See [`VonNeumannEntropy`](@ref) for the partial-trace
convention.
"""
function fetch(
    model::S1Heisenberg1D, q::RenyiEntropy, bc::OBC; ℓ::Int, beta::Real=Inf, kwargs...
)
    N = bc.N
    1 ≤ ℓ ≤ N - 1 ||
        throw(ArgumentError("RenyiEntropy: ℓ must satisfy 1 ≤ ℓ ≤ N - 1 (got ℓ=$ℓ, N=$N)"))
    α = q.α
    kernel = _s1_thermal_kernel(model, N, beta)
    ρA = _s1_partial_trace_A(kernel.ρ, ℓ, N)
    λs = real.(eigvals(Hermitian((ρA + ρA') / 2)))
    s = 0.0
    @inbounds for λ in λs
        if λ > 1e-15
            s += λ^α
        end
    end
    return log(s) / (1 - α)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Infinite-system literature values
# ═══════════════════════════════════════════════════════════════════════════════

"""
    native_energy_granularity(::S1Heisenberg1D, ::Infinite) -> :per_site

The infinite spin-1 Heisenberg chain has a well-defined per-site
ground-state energy density `e₀ ≈ -1.401484 J` (White-Huse 1993) but no
finite total energy.  Per-site is the only meaningful granularity.
"""
native_energy_granularity(::S1Heisenberg1D, ::Infinite) = :per_site

"""
    fetch(model::S1Heisenberg1D, ::Energy{:per_site}, ::Infinite) -> Float64

Ground-state energy per site of the spin-1 antiferromagnetic Heisenberg
chain (Haldane chain) in the thermodynamic limit:

    e₀ ≈ -1.40148403897 J

Numerical value from the high-precision DMRG study of
S. R. White & D. A. Huse, Phys. Rev. B **48**, 3844 (1993).  No closed-
form expression is known; `reliability=:medium` in the registry.
"""
function fetch(model::S1Heisenberg1D, ::Energy{:per_site}, ::Infinite; kwargs...)
    return -1.40148403897 * model.J
end

"""
    fetch(model::S1Heisenberg1D, ::MassGap, ::Infinite) -> Float64

Bulk Haldane gap of the spin-1 antiferromagnetic Heisenberg chain,

    Δ_∞ ≈ 0.41048 J

(literature value, S. R. White & D. A. Huse, Phys. Rev. B **48**, 3844
(1993)).  No closed form; `reliability=:medium` in the registry.
"""
function fetch(model::S1Heisenberg1D, ::MassGap, ::Infinite; kwargs...)
    return 0.41048 * model.J
end
