# ─────────────────────────────────────────────────────────────────────────────
# Transverse Field Ising Model — exact real-time correlation functions
#
# Hamiltonian (OBC, N sites):
#
#   H = -J Σ_i σᶻ_i σᶻ_{i+1}  -  h Σ_i σˣ_i
#
# Jordan–Wigner mapping (standard): σˣ diagonal, σᶻ nonlocal.  Introducing
# Majorana operators
#
#   γ_{2i-1} = c_i + c_i†
#   γ_{2i}   = i (c_i† - c_i)            {γ_a, γ_b} = 2 δ_{ab}
#
# one finds
#
#   σˣ_i                = -i γ_{2i-1} γ_{2i}
#   σᶻ_i σᶻ_{i+1}       = -i γ_{2i}   γ_{2i+1}
#   σᶻ_i                = (-i)^{i-1} γ_1 γ_2 … γ_{2i-2} γ_{2i-1}  (JW string,
#                                                                  consecutive
#                                                                  Majoranas)
#
# (the single-operator expression contains an odd number of Majoranas — its
# expectation value in a Gaussian state vanishes, so ⟨σᶻ⟩_GS = 0 on OBC.)
#
# The Hamiltonian takes the quadratic Majorana form
#
#   H = (i/4) Σ_{ab} h_{ab} γ_a γ_b
#
# with h a 2N × 2N real antisymmetric matrix whose only non-zero entries
# (upper triangle) are
#
#   h_{2i-1, 2i} = 2h    (i = 1 … N)
#   h_{2i  , 2i+1} = 2J  (i = 1 … N-1)
#
# The Heisenberg-picture evolution reads
#
#   γ(t) = exp(h t) γ(0),          exp(h t) ∈ SO(2N)
#
# and the ground-state Majorana 2-point function is
#
#   ⟨γ_a γ_b⟩_GS = δ_{ab} + i Σ_{ab},       Σ = -i · sign(i h)
#
# (derivation: in the canonical basis where h is block-diagonal with 2×2
# blocks [0 Λ; -Λ 0], the ground-state block is `Σ̃_block = [0 1; -1 0]`,
# while a direct computation gives `i sign(i h)|block = [0 -1; 1 0]`, hence
# the overall minus sign).
# Higher-point correlators reduce to a Pfaffian via Wick's theorem.
#
# The quantities exposed to `fetch` are
#
#   :sz_sz_correlation  — ⟨σᶻ_i(t) σᶻ_j(0)⟩_GS          (Pfaffian, size 2(i+j)-2)
#   :sx_sx_correlation  — ⟨σˣ_i(t) σˣ_j(0)⟩_GS          (4×4 Pfaffian)
#   :sz_sz_spreading    — matrix C[it, ix] = ⟨σᶻ_{ix}(t) σᶻ_center(0)⟩
#                          over all sites × supplied `times`
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra: Hermitian, eigen, Diagonal, I

# ═══════════════════════════════════════════════════════════════════════════════
# Building blocks
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _majorana_ham(N, J, h) -> Matrix{Float64}

Return the 2N × 2N real antisymmetric Majorana Hamiltonian matrix `h`
of the OBC TFIM such that `H = (i/4) Σ h_{ab} γ_a γ_b`.
"""
function _majorana_ham(N::Int, J::Float64, h::Float64)
    M = zeros(Float64, 2N, 2N)
    @inbounds for i in 1:N
        M[2i - 1, 2i] = 2h
        M[2i, 2i - 1] = -2h
    end
    @inbounds for i in 1:(N - 1)
        M[2i, 2i + 1] = 2J
        M[2i + 1, 2i] = -2J
    end
    return M
end

"""
    _majorana_covariance_gs(h::AbstractMatrix) -> Matrix{Float64}

Ground-state Majorana covariance Σ of a quadratic Hamiltonian whose
Majorana matrix is `h`.  Σ is real antisymmetric, with
`⟨γ_a γ_b⟩_GS = δ_{ab} + i Σ_{ab}`.  The formula `Σ = -i · sign(i h)` is
used, evaluated through an eigendecomposition of the Hermitian matrix
`i h`.
"""
function _majorana_covariance_gs(h::AbstractMatrix{<:Real})
    M = im .* h                                # Hermitian 2N × 2N
    F = eigen(Hermitian((M + M') / 2))
    s = [λ > 0 ? 1.0 : (λ < 0 ? -1.0 : 0.0) for λ in F.values]
    sM = F.vectors * Diagonal(s) * F.vectors'
    Σ = real(-im .* sM)
    # Symmetrise to strict antisymmetry to cancel round-off drift.
    return (Σ - Σ') / 2
end

"""
    _majorana_thermal_covariance(h::AbstractMatrix, β::Real) -> Matrix{Float64}

Thermal-equilibrium Majorana covariance at inverse temperature `β`.  The
T = 0 ground-state formula `Σ_GS = -i sign(i h)` generalises to

    Σ(β) = -i tanh((β/2) · i h)

so that `⟨γ_a γ_b⟩_β = δ_{ab} + i Σ(β)_{ab}`.  In the canonical 2×2 BdG
block this gives the off-diagonal entry `tanh(βΛ/2)`, recovering the
Fermi–Dirac population of each quasiparticle and reducing to `±1` as
`β → ∞`.

`β = Inf` falls back to `_majorana_covariance_gs`.
"""
function _majorana_thermal_covariance(h::AbstractMatrix{<:Real}, β::Real)
    isinf(β) && return _majorana_covariance_gs(h)
    M = im .* h
    F = eigen(Hermitian((M + M') / 2))
    s = tanh.((β / 2) .* F.values)
    sM = F.vectors * Diagonal(s) * F.vectors'
    Σ = real(-im .* sM)
    return (Σ - Σ') / 2
end

"""
    _majorana_evolution(h::AbstractMatrix, t::Real) -> Matrix{Float64}

`exp(h * t)` for a real antisymmetric `h`, returned as a real
orthogonal matrix.  Tiny imaginary noise from the matrix exponential
is projected away.
"""
function _majorana_evolution(h::AbstractMatrix{<:Real}, t::Real)
    if t == 0
        return Matrix{Float64}(I, size(h)...)
    end
    return real(exp(h .* t))
end

# Majorana index list for σᶻ_k = (-i)^{k-1} γ_1 γ_2 … γ_{2k-2} γ_{2k-1}.
# (Equivalent to σˣ_k in the standard Kitaev convention; QAtlas's H uses
# Z-Z coupling, X-field, which is unitarily related by a Hadamard rotation.)
_sz_majorana_indices(k::Int) = collect(1:(2k - 1))

# For σˣ_k = -i γ_{2k-1} γ_{2k}, the two Majoranas are (2k-1, 2k).
_sx_majorana_indices(k::Int) = Int[2k - 1, 2k]

"""
    _build_wick_matrix(idx_t, idx_0, Σ, R, RΣ) -> Matrix{ComplexF64}

Assemble the Wick contraction matrix `F` for a product

    γ_{idx_t[1]}(t) γ_{idx_t[2]}(t) … γ_{idx_0[1]}(0) γ_{idx_0[2]}(0) …

from the ground-state Majorana covariance `Σ`, the time evolution
matrix `R = exp(h t)`, and the precomputed product `RΣ = R * Σ`.
`F` is complex antisymmetric; its Pfaffian gives `⟨…⟩_GS`.
"""
function _build_wick_matrix(
    idx_t::Vector{Int},
    idx_0::Vector{Int},
    Σ::AbstractMatrix,
    R::AbstractMatrix,
    RΣ::AbstractMatrix,
)
    na = length(idx_t)
    nb = length(idx_0)
    M = na + nb
    F = zeros(ComplexF64, M, M)

    # (t, t) block
    @inbounds for k in 1:na, l in (k + 1):na
        a = idx_t[k]
        b = idx_t[l]
        val = (a == b ? 1.0 : 0.0) + im * Σ[a, b]
        F[k, l] = val
        F[l, k] = -val
    end

    # (t, 0) cross block
    @inbounds for k in 1:na, l in 1:nb
        a = idx_t[k]
        b = idx_0[l]
        val = R[a, b] + im * RΣ[a, b]
        F[k, na + l] = val
        F[na + l, k] = -val
    end

    # (0, 0) block
    @inbounds for k in 1:nb, l in (k + 1):nb
        a = idx_0[k]
        b = idx_0[l]
        val = (a == b ? 1.0 : 0.0) + im * Σ[a, b]
        F[na + k, na + l] = val
        F[na + l, na + k] = -val
    end

    return F
end

# ═══════════════════════════════════════════════════════════════════════════════
# σᶻ σᶻ and σˣ σˣ real-time correlation functions
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _sz_sz_corr(N, J, h, i, j, t; β = Inf) -> ComplexF64

`⟨σᶻ_i(t) σᶻ_j(0)⟩_β` for the OBC TFIM at inverse temperature `β`
(default `Inf` ⇒ ground state).

Implementation: Wick / Pfaffian over the `(2i-1) + (2j-1) = 2(i+j)-2`
Majorana operators constituting the product.  The overall phase is
`(-i)^{i+j-2}`.  The thermal generalisation enters only through the
Majorana 2-point function

    ⟨γ_a γ_b⟩_β = δ_{ab} + i Σ(β)_{ab},   Σ(β) = -i tanh((β/2) i h),

so the body of the Pfaffian assembly is unchanged.
"""
function _sz_sz_corr(N::Int, J::Float64, h::Float64, i::Int, j::Int, t::Real; β::Real=Inf)
    (1 ≤ i ≤ N && 1 ≤ j ≤ N) || throw(ArgumentError("site indices out of range"))
    hmat = _majorana_ham(N, J, h)
    Σ = _majorana_thermal_covariance(hmat, β)
    R = _majorana_evolution(hmat, t)
    return _sz_sz_corr_from_cached(Σ, R, i, j)
end

function _sz_sz_corr_from_cached(Σ::AbstractMatrix, R::AbstractMatrix, i::Int, j::Int)
    RΣ = R * Σ
    idx_t = _sz_majorana_indices(i)
    idx_0 = _sz_majorana_indices(j)
    F = _build_wick_matrix(idx_t, idx_0, Σ, R, RΣ)
    pf = pfaffian(F)
    # σᶻ_k = (-i)^{k-1} γ_1 … γ_{2k-1}, so the prefactor of the product is
    # (-i)^{i+j-2}.
    return ((-im)^(i + j - 2)) * pf
end

"""
    _sx_sx_corr(N, J, h, i, j, t; β = Inf) -> ComplexF64

`⟨σˣ_i(t) σˣ_j(0)⟩_β` for the OBC TFIM.  Reduces to a 4×4 Pfaffian since
`σˣ_k = -i γ_{2k-1} γ_{2k}`.  `β = Inf` ⇒ ground state.
"""
function _sx_sx_corr(N::Int, J::Float64, h::Float64, i::Int, j::Int, t::Real; β::Real=Inf)
    (1 ≤ i ≤ N && 1 ≤ j ≤ N) || throw(ArgumentError("site indices out of range"))
    hmat = _majorana_ham(N, J, h)
    Σ = _majorana_thermal_covariance(hmat, β)
    R = _majorana_evolution(hmat, t)
    return _sx_sx_corr_from_cached(Σ, R, i, j)
end

function _sx_sx_corr_from_cached(Σ::AbstractMatrix, R::AbstractMatrix, i::Int, j::Int)
    RΣ = R * Σ
    idx_t = _sx_majorana_indices(i)
    idx_0 = _sx_majorana_indices(j)
    F = _build_wick_matrix(idx_t, idx_0, Σ, R, RΣ)
    pf = pfaffian(F)
    # σˣ carries factor (-i)·(-i) = -1 from the two -i γγ prefactors.
    return -pf
end

# ═══════════════════════════════════════════════════════════════════════════════
# Single-operator expectations
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _sx_expect(Σ, i) -> Float64

`⟨σˣ_i⟩_GS = Σ[2i-1, 2i]`, since `σˣ_i = -i γ_{2i-1} γ_{2i}`.
"""
_sx_expect(Σ::AbstractMatrix, i::Int) = Σ[2i - 1, 2i]

# ⟨σᶻ_i⟩_GS = 0 for OBC (odd-product expectation in Gaussian state).
_sz_expect(::AbstractMatrix, ::Int) = 0.0

# ═══════════════════════════════════════════════════════════════════════════════
# Spreading-correlation convenience
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _sz_sz_spreading(N, J, h, center, times) -> Matrix{ComplexF64}

Return `C[it, ix] = ⟨σᶻ_ix(t_it) σᶻ_center(0)⟩_GS` for `ix ∈ 1:N`,
`t_it ∈ 1:length(times)`.  Reuses the Majorana Hamiltonian, covariance,
and (per time-step) evolution matrix, so the cost is
`O(length(times) · N · M³)` with `M = 2(center + N) - 2`.
"""
function _sz_sz_spreading(
    N::Int, J::Float64, h::Float64, center::Int, times::AbstractVector{<:Real}; β::Real=Inf
)
    (1 ≤ center ≤ N) || throw(ArgumentError("center site out of range"))
    hmat = _majorana_ham(N, J, h)
    Σ = _majorana_thermal_covariance(hmat, β)
    nt = length(times)
    C = zeros(ComplexF64, nt, N)
    for (it, t) in enumerate(times)
        R = _majorana_evolution(hmat, t)
        for ix in 1:N
            C[it, ix] = _sz_sz_corr_from_cached(Σ, R, ix, center)
        end
    end
    return C
end

# ═══════════════════════════════════════════════════════════════════════════════
# Static (equal-time) thermal correlators and structure factor
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _sz_sz_static_thermal(N, J, h, β; i = nothing, j = nothing) -> Matrix{Float64}

Static `⟨σᶻ_i σᶻ_j⟩_β` for the OBC TFIM at inverse temperature `β`.
If both `i` and `j` are given, returns a single value (wrapped in a 1×1
matrix); otherwise returns the full N×N matrix of equal-time correlators.
"""
function _sz_sz_static_thermal(
    N::Int,
    J::Float64,
    h::Float64,
    β::Real;
    i::Union{Int,Nothing}=nothing,
    j::Union{Int,Nothing}=nothing,
)
    hmat = _majorana_ham(N, J, h)
    Σ = _majorana_thermal_covariance(hmat, β)
    R = Matrix{Float64}(I, size(hmat)...)
    if i !== nothing && j !== nothing
        v = real(_sz_sz_corr_from_cached(Σ, R, i, j))
        return fill(v, 1, 1)
    end
    C = zeros(Float64, N, N)
    for a in 1:N, b in a:N
        v = real(_sz_sz_corr_from_cached(Σ, R, a, b))
        C[a, b] = v
        C[b, a] = v
    end
    return C
end

"""
    _zz_static_structure_factor(N, J, h, β, q) -> Float64

`S_zz(q, β) = (1/N) Σ_{i,j} e^{-i q (i-j)} ⟨σᶻ_i σᶻ_j⟩_β` evaluated by direct
double sum from the thermal Pfaffian correlator.  For OBC the lattice lacks
translation invariance so this is the boundary-aware definition; in the bulk
of a long enough chain it converges to the translation-invariant value.
"""
function _zz_static_structure_factor(N::Int, J::Float64, h::Float64, β::Real, q::Real)
    C = _sz_sz_static_thermal(N, J, h, β)
    s = 0.0 + 0.0im
    for i in 1:N, j in 1:N
        s += exp(-im * q * (i - j)) * C[i, j]
    end
    return real(s) / N
end

"""
    _zz_uniform_susceptibility(N, J, h, β) -> Float64

Uniform (q = 0) longitudinal susceptibility per site,

    χ_zz(β) = (β/N) Σ_{i,j} ⟨σᶻ_i σᶻ_j⟩_β,

obtained from a finite-temperature direct double sum (assumes ⟨σᶻ⟩_β = 0 in
the OBC ground manifold of the TFIM, which holds for any `h ≠ 0` finite N).
This is the static *isothermal* susceptibility — the fluctuation-dissipation
form `χ = β · ⟨M²⟩_c / N` for the magnetisation `M = Σᵢ σᶻᵢ`.
"""
function _zz_uniform_susceptibility(N::Int, J::Float64, h::Float64, β::Real)
    C = _sz_sz_static_thermal(N, J, h, β)
    return β * sum(C) / N
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch dispatch
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::TFIM, ::ZZCorrelation{:dynamic}, bc::OBC;
          i::Int, j::Int, t::Float64, beta::Float64 = Inf) -> ComplexF64

Exact `⟨σᶻ_i(t) σᶻ_j(0)⟩_β` for the OBC TFIM.  `beta = Inf` (the default)
gives the ground-state result.  `N` comes from `bc.N`.
"""
function fetch(
    model::TFIM,
    ::ZZCorrelation{:dynamic},
    bc::OBC;
    i::Int,
    j::Int,
    t::Float64,
    beta::Real=Inf,
    kwargs...,
)
    N = _bc_size(bc, kwargs)
    return _sz_sz_corr(N, model.J, model.h, i, j, t; β=beta)
end

"""
    fetch(model::TFIM, ::XXCorrelation{:dynamic}, bc::OBC;
          i::Int, j::Int, t::Float64, beta::Float64 = Inf) -> ComplexF64

Exact `⟨σˣ_i(t) σˣ_j(0)⟩_β` for the OBC TFIM.
"""
function fetch(
    model::TFIM,
    ::XXCorrelation{:dynamic},
    bc::OBC;
    i::Int,
    j::Int,
    t::Float64,
    beta::Real=Inf,
    kwargs...,
)
    N = _bc_size(bc, kwargs)
    return _sx_sx_corr(N, model.J, model.h, i, j, t; β=beta)
end

"""
    fetch(model::TFIM, ::ZZCorrelation{:lightcone}, bc::OBC;
          center::Int, times::AbstractVector{<:Real}, beta::Float64 = Inf) -> Matrix{ComplexF64}

Exact spreading correlation `C[it, ix] = ⟨σᶻ_ix(t_it) σᶻ_center(0)⟩_β` for all
sites `ix ∈ 1:N` and `t_it ∈ times`.
"""
function fetch(
    model::TFIM,
    ::ZZCorrelation{:lightcone},
    bc::OBC;
    center::Int,
    times::AbstractVector{<:Real},
    beta::Real=Inf,
    kwargs...,
)
    N = _bc_size(bc, kwargs)
    return _sz_sz_spreading(N, model.J, model.h, center, times; β=beta)
end

"""
    fetch(model::TFIM, ::ZZCorrelation{:static}, bc::OBC;
          beta::Float64, [i::Int, j::Int]) -> Matrix{Float64} or Float64

Static (equal-time) thermal correlator `⟨σᶻ_i σᶻ_j⟩_β` for the OBC TFIM.
With both `i` and `j` given returns a scalar; otherwise returns the full
N×N matrix.
"""
function fetch(
    model::TFIM,
    ::ZZCorrelation{:static},
    bc::OBC;
    beta::Float64,
    i::Union{Int,Nothing}=nothing,
    j::Union{Int,Nothing}=nothing,
    kwargs...,
)
    N = _bc_size(bc, kwargs)
    if i !== nothing && j !== nothing
        return _sz_sz_static_thermal(N, model.J, model.h, beta; i=i, j=j)[1, 1]
    end
    return _sz_sz_static_thermal(N, model.J, model.h, beta)
end

"""
    fetch(model::TFIM, ::ZZCorrelation{:connected}, bc::OBC;
          beta::Float64, [i::Int, j::Int]) -> Matrix{Float64} or Float64

Connected static thermal correlator
`C^c_{ij} = ⟨σᶻ_i σᶻ_j⟩_β − ⟨σᶻ_i⟩_β ⟨σᶻ_j⟩_β` for the OBC TFIM.

In the OBC TFIM the Z₂ symmetry σᶻ → −σᶻ is unbroken at any finite N
(Gaussian state of the JW fermions; odd-product expectation
vanishes), so `⟨σᶻ_i⟩_β = 0` and the connected correlator coincides
with the bare static one.  This method therefore re-uses the
`:static` routine and is provided as a separate dispatch for caller
clarity / API completeness.

If at some point a TFIM variant breaks Z₂ explicitly (e.g. by adding
a longitudinal field), the implementation will still be correct
provided the per-site `⟨σᶻ_i⟩_β` is taken from
[`MagnetizationZLocal`](@ref) and subtracted off — see the comment
in the source.
"""
function fetch(
    model::TFIM,
    ::ZZCorrelation{:connected},
    bc::OBC;
    beta::Float64,
    i::Union{Int,Nothing}=nothing,
    j::Union{Int,Nothing}=nothing,
    kwargs...,
)
    # In TFIM at OBC ⟨σᶻ_i⟩ = 0 by Z₂ symmetry.  The static routine already
    # returns the connected correlator under that assumption.  Wire to it
    # directly; if a future TFIM variant adds explicit symmetry breaking,
    # subtract `mz[i] * mz[j]` here using `MagnetizationZLocal`.
    return fetch(model, ZZCorrelation{:static}(), bc; beta=beta, i=i, j=j, kwargs...)
end

"""
    fetch(model::TFIM, ::ZZStructureFactor, bc::OBC;
          beta::Float64, q::Real) -> Float64

Static structure factor `S_zz(q, β)` for the OBC TFIM at wave vector `q`.
"""
function fetch(model::TFIM, ::ZZStructureFactor, bc::OBC; beta::Float64, q::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    return _zz_static_structure_factor(N, model.J, model.h, beta, q)
end

"""
    fetch(model::TFIM, ::SusceptibilityZZ, bc::OBC;
          beta::Float64) -> Float64

Static uniform longitudinal (`q = 0`) susceptibility per site,

    χ_zz(β) = (β/N) Σ_{i,j} ⟨σᶻ_i σᶻ_j⟩_β.
"""
function fetch(model::TFIM, ::SusceptibilityZZ, bc::OBC; beta::Float64, kwargs...)
    N = _bc_size(bc, kwargs)
    J = model.J
    h = model.h
    return _zz_uniform_susceptibility(N, J, h, beta)
end
