# ─────────────────────────────────────────────────────────────────────────────
# XXZ chain (1D, spin-1/2) — finite-N OBC thermal observables via dense ED.
#
# Hamiltonian (spin convention):
#
#   H = J Σᵢ [Sˣᵢ Sˣᵢ₊₁ + Sʸᵢ Sʸᵢ₊₁ + Δ Sᶻᵢ Sᶻᵢ₊₁]
#     = (J/4) Σᵢ [σˣ σˣ + σʸ σʸ + Δ σᶻ σᶻ]
#
# Observables here are returned in the *Pauli* (`σ`) convention to stay
# numerically compatible with the TFIM spec sheet, e.g.
#
#   ⟨σᵅ_i⟩ = 2 ⟨Sᵅ_i⟩,    χ_αα = β Var(Σᵢ σᵅᵢ) / N.
#
# Method.  Hilbert space dimension `D = 2^N` (capped at `N ≤
# _MAX_ED_SITES = 12`).  We diagonalise H once in `eigen(Hermitian(H))`
# and reuse the eigendecomposition for every thermal expectation:
#
#     ⟨A⟩_β = Σₙ wₙ ⟨n|A|n⟩,    wₙ = exp(-β eₙ) / Z.
#
# This is the same pattern as TFIM's `_majorana_thermal_covariance`,
# adapted to a non-Gaussian chain where we keep the full eigenbasis
# instead of projecting to a 2N×2N covariance.
#
# Cost: O(D³) for `eigen`, O(D²) per observable evaluation through the
# eigenbasis.  At N = 12 (D = 4096) the eigendecomposition is ~1 s on a
# laptop and every subsequent expectation is ~10 ms.
#
# Why a private helper module instead of `_ed_thermal_*` in
# `core/dense_ed.jl`?  The S1Heisenberg agent shares `dense_ed.jl`; we
# avoid editing that file and keep XXZ-specific helpers here so the two
# extension Tier-1 PRs do not collide.
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra: Hermitian, eigen, Eigen, Diagonal, diagm, diag, log, tr, I

# ═══════════════════════════════════════════════════════════════════════════════
# Internal: thermal kernel — one eigendecomposition reused across observables
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _xxz1d_thermal_kernel(model, N, beta) -> NamedTuple

Diagonalise the OBC XXZ Hamiltonian once and return everything that
downstream observables need:

- `evals`     :: Vector{Float64}      — eigenvalues
- `evecs`     :: Matrix{ComplexF64}   — eigenvectors (columns)
- `weights`   :: Vector{Float64}      — Boltzmann weights wₙ, normalized so
                                       Σ wₙ = 1  (uses an `emin` shift)
- `H`         :: Matrix{ComplexF64}   — the Hamiltonian (kept for callers
                                       that need it, e.g. EnergyLocal)

Callers that only need the spectrum can ignore `evecs`.
"""
function _xxz1d_thermal_kernel(model::XXZ1D, N::Int, β::Real)
    H = _xxz1d_hamiltonian_matrix(model, N)
    F = eigen(Hermitian(H))
    evals = F.values
    evecs = F.vectors
    emin = minimum(evals)
    ws = exp.(-β .* (evals .- emin))
    weights = ws ./ sum(ws)
    return (; evals=evals, evecs=evecs, weights=weights, H=H)
end

"""
    _xxz1d_thermal_expectation_op(F::NamedTuple, A::AbstractMatrix) -> Float64

Compute `⟨A⟩_β = Σₙ wₙ ⟨n|A|n⟩` for a Hermitian operator `A`, using a
`_xxz1d_thermal_kernel` result `F`.  Returns the real part (assumes the
imaginary part is round-off only).
"""
function _xxz1d_thermal_expectation_op(F::NamedTuple, A::AbstractMatrix)
    Adiag = real.(diag(F.evecs' * A * F.evecs))
    return sum(F.weights .* Adiag)
end

"""
    _xxz1d_thermal_density_matrix(F::NamedTuple) -> Matrix{ComplexF64}

Build the density matrix `ρ = exp(-βH)/Z = U diag(w) U†` from a kernel
result.  Allocates a `2^N × 2^N` matrix; only call for entanglement
quantities that actually need ρ explicitly.
"""
function _xxz1d_thermal_density_matrix(F::NamedTuple)
    return F.evecs * (Diagonal(ComplexF64.(F.weights)) * F.evecs')
end

# Pauli-string single-site operators σᵅ_i
@inline _xxz1d_sx(N::Int, i::Int) = _pauli_string(N, i => _σx)
@inline _xxz1d_sy(N::Int, i::Int) = _pauli_string(N, i => _σy)
@inline _xxz1d_sz(N::Int, i::Int) = _pauli_string(N, i => _σz)

# Pauli pair σᵅ_i σᵅ_j (i ≠ j)
function _xxz1d_sasa(N::Int, i::Int, j::Int, σα::AbstractMatrix)
    i == j && return _pauli_string(N, i => σα * σα)  # = I when σα is Pauli
    return _pauli_string(N, i => σα, j => σα)
end

# Total magnetisation Σᵢ σᵅ_i
function _xxz1d_total_M(N::Int, σα::AbstractMatrix)
    D = 2^N
    M = zeros(ComplexF64, D, D)
    for i in 1:N
        M .+= _pauli_string(N, i => σα)
    end
    return M
end

# ═══════════════════════════════════════════════════════════════════════════════
# Native granularity declarations for the thermodynamic potentials we expose
# ═══════════════════════════════════════════════════════════════════════════════
#
# `Energy` is already declared in XXZ.jl (`:total` at OBC).  The other
# thermal scalars (FreeEnergy / ThermalEntropy / SpecificHeat) follow
# the TFIM precedent of returning per-site by default; no
# native_energy_granularity entry needed for those.

# ═══════════════════════════════════════════════════════════════════════════════
# Thermodynamic potentials (per-site) — F, S, C
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::XXZ1D, ::FreeEnergy, ::OBC; beta) -> Float64

Per-site Helmholtz free energy `f(β) = -log Z / (Nβ)` of the spin-½ OBC
XXZ chain at finite `N ≤ $(_MAX_ED_SITES)`, computed by dense ED.
"""
function fetch(model::XXZ1D, ::FreeEnergy, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    H = _xxz1d_hamiltonian_matrix(model, N)
    return _ed_thermal_free_energy(H, beta) / N
end

"""
    fetch(model::XXZ1D, ::ThermalEntropy, ::OBC; beta) -> Float64

Per-site Gibbs entropy `s(β) = β · (ε - f)` of the OBC XXZ chain.
"""
function fetch(model::XXZ1D, ::ThermalEntropy, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    H = _xxz1d_hamiltonian_matrix(model, N)
    return _ed_thermal_entropy(H, beta) / N
end

"""
    fetch(model::XXZ1D, ::SpecificHeat, ::OBC; beta) -> Float64

Per-site heat capacity `c(β) = β² · Var(H) / N`, computed exactly from
the energy variance in the eigenbasis (no numerical differentiation).
"""
function fetch(model::XXZ1D, ::SpecificHeat, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    H = _xxz1d_hamiltonian_matrix(model, N)
    return _ed_thermal_specific_heat(H, beta) / N
end

# ═══════════════════════════════════════════════════════════════════════════════
# Magnetisation (per-site, Pauli convention)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::XXZ1D, ::MagnetizationX, ::OBC; beta) -> Float64

Per-site bulk magnetisation `⟨Σᵢ σˣᵢ⟩_β / N` of the OBC XXZ chain.
"""
function fetch(model::XXZ1D, ::MagnetizationX, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    F = _xxz1d_thermal_kernel(model, N, beta)
    Mx = _xxz1d_total_M(N, _σx)
    return _xxz1d_thermal_expectation_op(F, Mx) / N
end

"""
    fetch(model::XXZ1D, ::MagnetizationY, ::OBC; beta) -> Float64

Per-site bulk magnetisation `⟨Σᵢ σʸᵢ⟩_β / N`.  Identically zero for the
XXZ Hamiltonian (the eigenvectors of a real-symmetric H — the XXZ
matrix in the σᶻ-product basis is real symmetric — give purely
imaginary expectations of σʸ that cancel mode-by-mode).  Returned by
explicit calculation rather than hard-coded zero so the caller still
sees the dense-ED noise floor.
"""
function fetch(model::XXZ1D, ::MagnetizationY, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    F = _xxz1d_thermal_kernel(model, N, beta)
    My = _xxz1d_total_M(N, _σy)
    return _xxz1d_thermal_expectation_op(F, My) / N
end

"""
    fetch(model::XXZ1D, ::MagnetizationZ, ::OBC; beta) -> Float64

Per-site bulk magnetisation `⟨Σᵢ σᶻᵢ⟩_β / N`.  Σᵢ σᶻᵢ commutes with `H`
(U(1) symmetry), so the thermal average is `Tr(M_z exp(-βH))/Z`; for
even `N` and any real `β` the σᶻ-product basis groups symmetrically
between sectors of opposite total Sᶻ and the average is zero.
For odd `N` the unique smallest-Sᶻ-magnitude sector is half-filled
±½ (still S_z = ±½), so the average is again zero up to round-off.
"""
function fetch(model::XXZ1D, ::MagnetizationZ, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    F = _xxz1d_thermal_kernel(model, N, beta)
    Mz = _xxz1d_total_M(N, _σz)
    return _xxz1d_thermal_expectation_op(F, Mz) / N
end

# ─── Site-resolved magnetisations (Vector of length N) ────────────────────

"""
    fetch(model::XXZ1D, ::MagnetizationXLocal, ::OBC; beta) -> Vector{Float64}

Site-resolved `[⟨σˣ_i⟩_β for i = 1:N]`.  Identically zero up to
dense-ED round-off because `σˣ_i` flips a single Sᶻ and the XXZ
Hamiltonian conserves total Sᶻ.
"""
function fetch(model::XXZ1D, ::MagnetizationXLocal, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    F = _xxz1d_thermal_kernel(model, N, beta)
    return Float64[_xxz1d_thermal_expectation_op(F, _xxz1d_sx(N, i)) for i in 1:N]
end

"""
    fetch(model::XXZ1D, ::MagnetizationYLocal, ::OBC; beta) -> Vector{Float64}

Site-resolved `[⟨σʸ_i⟩_β for i = 1:N]`.  Identically zero by the same
U(1) argument as `MagnetizationXLocal` plus parity (σʸ_i is purely
imaginary in the σᶻ basis).
"""
function fetch(model::XXZ1D, ::MagnetizationYLocal, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    F = _xxz1d_thermal_kernel(model, N, beta)
    return Float64[_xxz1d_thermal_expectation_op(F, _xxz1d_sy(N, i)) for i in 1:N]
end

"""
    fetch(model::XXZ1D, ::MagnetizationZLocal, ::OBC; beta) -> Vector{Float64}

Site-resolved `[⟨σᶻ_i⟩_β for i = 1:N]`.  Each ⟨σᶻ_i⟩ is identically
zero up to round-off in the canonical Boltzmann ensemble (sectors of
opposite S_z come in equal weight pairs).
"""
function fetch(model::XXZ1D, ::MagnetizationZLocal, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    F = _xxz1d_thermal_kernel(model, N, beta)
    return Float64[_xxz1d_thermal_expectation_op(F, _xxz1d_sz(N, i)) for i in 1:N]
end

# ─── Local energy density (length N, Σ ε_i = ⟨H⟩) ─────────────────────────

"""
    fetch(model::XXZ1D, ::EnergyLocal, ::OBC; beta) -> Vector{Float64}

Site-local energy density `ε_i` of the OBC XXZ chain at inverse
temperature `beta`, defined so that `Σᵢ ε_i = ⟨H⟩_β`.  Each bond
`b_{i,i+1} = (J/4)(σˣσˣ + σʸσʸ + Δ σᶻσᶻ)` is split symmetrically
between its two endpoints:

    ε_i = ½ ⟨b_{i-1,i}⟩_β + ½ ⟨b_{i,i+1}⟩_β

with the missing bonds at `i = 1` / `i = N` taken to be zero.
"""
function fetch(model::XXZ1D, ::EnergyLocal, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    F = _xxz1d_thermal_kernel(model, N, beta)
    J, Δ = model.J, model.Δ
    prefac = J / 4
    bonds = Vector{Float64}(undef, N - 1)
    for i in 1:(N - 1)
        bxx = _pauli_string(N, i => _σx, i + 1 => _σx)
        byy = _pauli_string(N, i => _σy, i + 1 => _σy)
        bzz = _pauli_string(N, i => _σz, i + 1 => _σz)
        bxx_v = _xxz1d_thermal_expectation_op(F, bxx)
        byy_v = _xxz1d_thermal_expectation_op(F, byy)
        bzz_v = _xxz1d_thermal_expectation_op(F, bzz)
        bonds[i] = prefac * (bxx_v + byy_v + Δ * bzz_v)
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
# Susceptibilities (per-site, σ-convention variance)
# ═══════════════════════════════════════════════════════════════════════════════

# Common helper: χ_αα(β) = β · (⟨M_α²⟩ - ⟨M_α⟩²) / N
function _xxz1d_uniform_susceptibility(F::NamedTuple, N::Int, σα::AbstractMatrix, β::Real)
    Mα = _xxz1d_total_M(N, σα)
    M1 = _xxz1d_thermal_expectation_op(F, Mα)
    M2 = _xxz1d_thermal_expectation_op(F, Mα * Mα)
    return β * (M2 - M1^2) / N
end

"""
    fetch(model::XXZ1D, ::SusceptibilityXX, ::OBC; beta) -> Float64

Static transverse susceptibility per site
`χ_xx(β) = (β/N) · (⟨M_x²⟩ - ⟨M_x⟩²)`, with `M_x = Σᵢ σˣᵢ`.
"""
function fetch(model::XXZ1D, ::SusceptibilityXX, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    F = _xxz1d_thermal_kernel(model, N, beta)
    return _xxz1d_uniform_susceptibility(F, N, _σx, beta)
end

"""
    fetch(model::XXZ1D, ::SusceptibilityYY, ::OBC; beta) -> Float64

Static y-axis susceptibility per site, `χ_yy(β) = (β/N) Var(M_y)`.
By the SU(2) → U(1) reduction the XXZ model has `⟨M_y⟩ = 0`, so
this is `β ⟨M_y²⟩ / N`.
"""
function fetch(model::XXZ1D, ::SusceptibilityYY, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    F = _xxz1d_thermal_kernel(model, N, beta)
    return _xxz1d_uniform_susceptibility(F, N, _σy, beta)
end

"""
    fetch(model::XXZ1D, ::SusceptibilityZZ, ::OBC; beta) -> Float64

Static longitudinal susceptibility per site, `χ_zz(β) = (β/N) Var(M_z)`.
At `Δ = 1` (Heisenberg) this equals `χ_xx = χ_yy` by SU(2) symmetry.
"""
function fetch(model::XXZ1D, ::SusceptibilityZZ, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    F = _xxz1d_thermal_kernel(model, N, beta)
    return _xxz1d_uniform_susceptibility(F, N, _σz, beta)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Two-point correlators (static, equal-time thermal)
# ═══════════════════════════════════════════════════════════════════════════════

# Common helper: returns ⟨σᵅ_i σᵅ_j⟩_β (or the connected version when the
# `connected = true` flag is set).
function _xxz1d_static_corr(
    F::NamedTuple, N::Int, σα::AbstractMatrix, i::Int, j::Int; connected::Bool=false
)
    1 ≤ i ≤ N || throw(ArgumentError("static correlation: i = $i out of 1..$N"))
    1 ≤ j ≤ N || throw(ArgumentError("static correlation: j = $j out of 1..$N"))
    if i == j
        # σᵅ² = I (Pauli) so ⟨σᵅᵢ σᵅᵢ⟩ = 1 always.
        c2 = 1.0
    else
        Op = _pauli_string(N, i => σα, j => σα)
        c2 = _xxz1d_thermal_expectation_op(F, Op)
    end
    if !connected
        return c2
    end
    Oi = _pauli_string(N, i => σα)
    Oj = _pauli_string(N, j => σα)
    ci = _xxz1d_thermal_expectation_op(F, Oi)
    cj = _xxz1d_thermal_expectation_op(F, Oj)
    return c2 - ci * cj
end

# Generate fetch methods for the four (mode, axis) combinations we support.
const _XXZ1D_CORR_AXES = ((XXCorrelation, _σx), (YYCorrelation, _σy), (ZZCorrelation, _σz))

for (CorrTy, σα) in _XXZ1D_CORR_AXES
    @eval begin
        function fetch(
            model::XXZ1D, ::$CorrTy{:static}, bc::OBC; beta::Real, i::Int, j::Int, kwargs...
        )
            N = _bc_size(bc, kwargs)
            F = _xxz1d_thermal_kernel(model, N, beta)
            return _xxz1d_static_corr(F, N, $σα, i, j; connected=false)
        end

        function fetch(
            model::XXZ1D,
            ::$CorrTy{:connected},
            bc::OBC;
            beta::Real,
            i::Int,
            j::Int,
            kwargs...,
        )
            N = _bc_size(bc, kwargs)
            F = _xxz1d_thermal_kernel(model, N, beta)
            return _xxz1d_static_corr(F, N, $σα, i, j; connected=true)
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Mass gap (OBC)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::XXZ1D, ::MassGap, ::OBC) -> Float64

Energy gap `E₁ - E₀` between the two lowest eigenvalues of the OBC
XXZ Hamiltonian at finite `N ≤ $(_MAX_ED_SITES)`.
"""
function fetch(model::XXZ1D, ::MassGap, bc::OBC; kwargs...)
    N = _bc_size(bc, kwargs)
    H = _xxz1d_hamiltonian_matrix(model, N)
    evals = sort(real.(eigvals(Hermitian(H))))
    return evals[2] - evals[1]
end

"""
    fetch(model::XXZ1D, ::MassGap, ::Infinite) -> Float64

Mass gap of the spin-½ XXZ chain in the thermodynamic limit:

- Critical regime `-1 < Δ ≤ 1`: gapless Luttinger liquid, returns `0.0`.
- Gapped regimes (`Δ > 1` antiferromagnetic Ising-like, `Δ < -1`
  ferromagnetic Ising-like): closed-form gap is non-trivial (Bethe
  ansatz integrals), not yet implemented; returns `NaN` with a
  warning.
"""
function fetch(model::XXZ1D, ::MassGap, ::Infinite; kwargs...)
    Δ = model.Δ
    if -1.0 < Δ ≤ 1.0
        return 0.0
    end
    @warn "XXZ1D MassGap (Infinite) gapped-regime closed-form not yet implemented; " *
        "use OBC at small N for an ED reference." Δ = Δ
    return NaN
end

# ═══════════════════════════════════════════════════════════════════════════════
# Entanglement entropies (mixed-state partial trace)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _xxz1d_partial_trace_B(ρ, ℓ, N) -> Matrix{ComplexF64}

Trace out the right `N - ℓ` sites from a `2^N × 2^N` density matrix
`ρ`, leaving the reduced state `ρ_A` on sites `1..ℓ` (a `2^ℓ × 2^ℓ`
matrix).  Implementation: reshape `ρ` to a 4-tensor
`(d_A, d_B, d_A, d_B)` and contract the two `d_B` legs.
"""
function _xxz1d_partial_trace_B(ρ::AbstractMatrix, ℓ::Int, N::Int)
    dA = 2^ℓ
    dB = 2^(N - ℓ)
    R = reshape(ρ, (dA, dB, dA, dB))
    ρA = zeros(eltype(ρ), dA, dA)
    @inbounds for a in 1:dA, ap in 1:dA, b in 1:dB
        ρA[a, ap] += R[a, b, ap, b]
    end
    return ρA
end

"""
    _xxz1d_reduced_density_matrix(model, N, ℓ, β) -> Matrix{ComplexF64}

Reduced density matrix of the first `ℓ` sites at inverse temperature `β`.
For `β = Inf` we use the ground-state pure state `|0⟩` (the lowest
eigenvector of `H`); for finite `β` we build the full thermal density
matrix and partial-trace.

The full thermal path costs `O(D²)` memory (D = 2^N); at the
`_MAX_ED_SITES = 12` ceiling that's a 4096×4096 complex matrix
(~256 MB), still cheap.
"""
function _xxz1d_reduced_density_matrix(model::XXZ1D, N::Int, ℓ::Int, β::Real)
    1 ≤ ℓ ≤ N - 1 ||
        throw(ArgumentError("entanglement: ℓ must satisfy 1 ≤ ℓ ≤ N-1; got ℓ=$ℓ, N=$N"))
    if isinf(β)
        H = _xxz1d_hamiltonian_matrix(model, N)
        F = eigen(Hermitian(H))
        ψ = F.vectors[:, 1]
        # ρ_A = Tr_B |ψ⟩⟨ψ|  =  reshape(ψ, (dA, dB)) · adjoint
        dA = 2^ℓ
        dB = 2^(N - ℓ)
        Ψ = reshape(ψ, (dA, dB))
        return Ψ * Ψ'
    end
    F = _xxz1d_thermal_kernel(model, N, β)
    ρ = _xxz1d_thermal_density_matrix(F)
    return _xxz1d_partial_trace_B(ρ, ℓ, N)
end

"""
    fetch(model::XXZ1D, ::VonNeumannEntropy, ::OBC; ℓ, beta=Inf) -> Float64

Von Neumann entanglement entropy `S = -Tr ρ_A log ρ_A` of the first `ℓ`
sites of the OBC XXZ chain at inverse temperature `beta` (or the
ground state when `beta = Inf`).

Computed by exact ED + partial trace; cost `O(2^{2N})` memory and
`O(2^{3ℓ})` for the `eigen` of `ρ_A`.  Capped by `_MAX_ED_SITES`.
"""
function fetch(
    model::XXZ1D, ::VonNeumannEntropy, bc::OBC; ℓ::Int, beta::Real=Inf, kwargs...
)
    N = _bc_size(bc, kwargs)
    ρA = _xxz1d_reduced_density_matrix(model, N, ℓ, beta)
    λ = eigvals(Hermitian(ρA))
    S = 0.0
    @inbounds for p in λ
        if real(p) > 1e-15
            S -= real(p) * log(real(p))
        end
    end
    return S
end

"""
    fetch(model::XXZ1D, q::RenyiEntropy, ::OBC; ℓ, beta=Inf) -> Float64

Rényi entropy of order `α = q.α` for the first `ℓ` sites,

    S_α = log(Tr ρ_A^α) / (1 - α).

`α = 1` is rejected at the `RenyiEntropy` constructor; use
`VonNeumannEntropy()` for that limit.
"""
function fetch(model::XXZ1D, q::RenyiEntropy, bc::OBC; ℓ::Int, beta::Real=Inf, kwargs...)
    N = _bc_size(bc, kwargs)
    ρA = _xxz1d_reduced_density_matrix(model, N, ℓ, beta)
    λ = real.(eigvals(Hermitian(ρA)))
    α = q.α
    s = 0.0
    @inbounds for p in λ
        # Clamp to ≥ 0 to handle round-off; skip exactly-zero modes.
        pp = max(p, 0.0)
        pp > 1e-15 || continue
        s += pp^α
    end
    return log(s) / (1 - α)
end
