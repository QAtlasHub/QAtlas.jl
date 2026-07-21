# ─────────────────────────────────────────────────────────────────────────────
# TFIM — y-axis observables on OBC (σʸ_i correlators + MagnetizationY +
# SusceptibilityYY).
#
# Conventions: see the JW header in `TFIM_dynamics.jl`.  In QAtlas's
# X-string convention the Majorana decomposition reads
#
#   σˣ_i = -i γ_{2i-1} γ_{2i}                        (local pair)
#   σᶻ_i = (-i)^{i-1} γ_1 γ_2 … γ_{2i-2} γ_{2i-1}    (string ending at 2i-1)
#   σʸ_i = i σˣ_i σᶻ_i = -(-i)^{i-1} γ_1 … γ_{2i-2} γ_{2i}    (string ending at 2i)
#
# σʸ_i carries the same Majorana count `2i-1` as σᶻ_i but its last
# Majorana is shifted from 2i-1 → 2i.  The product `σʸ_i σʸ_j` has the
# same (-i)^{i+j-2} prefactor as σᶻ_i σᶻ_j: each σʸ contributes
# `-(-i)^{i-1}` and `(-1)·(-1) = 1`.
#
# Single-operator expectations:
#   ⟨σʸ_i⟩_β = 0   (odd Majorana count → Wick vanishes in a Gaussian state)
#
# Static / dynamic 2-point correlators reduce to the same Pfaffian
# machinery as σᶻσᶻ — only the index lists change — so we re-use
# `_build_wick_matrix` from `TFIM_dynamics.jl`.
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Majorana index list for σʸ_k (= [1, 2, …, 2k-2, 2k]; length 2k-1).
# ═══════════════════════════════════════════════════════════════════════════════

function _sy_majorana_indices(k::Int)
    out = Vector{Int}(undef, 2k - 1)
    @inbounds for n in 1:(2k - 2)
        out[n] = n
    end
    @inbounds out[2k - 1] = 2k
    return out
end

# ═══════════════════════════════════════════════════════════════════════════════
# σʸ σʸ correlators — Pfaffian / Wick contraction
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _sy_sy_corr_from_cached(Σ, R, i, j) -> ComplexF64

`⟨σʸ_i(t) σʸ_j(0)⟩_β` from a precomputed thermal Majorana covariance `Σ`
and time-evolution matrix `R = exp(h·t)`.  Same structure as
`_sz_sz_corr_from_cached` (in `TFIM_dynamics.jl`) — only the index lists
differ.

The overall phase is `(-i)^{i+j-2}`: each σʸ_k contributes
`-(-i)^{k-1}`; the two minus signs cancel, leaving `(-i)^{i+j-2}`.
"""
function _sy_sy_corr_from_cached(Σ::AbstractMatrix, R::AbstractMatrix, i::Int, j::Int)
    RΣ = R * Σ
    idx_t = _sy_majorana_indices(i)
    idx_0 = _sy_majorana_indices(j)
    F = _build_wick_matrix(idx_t, idx_0, Σ, R, RΣ)
    pf = pfaffian(F)
    return ((-im)^(i + j - 2)) * pf
end

"""
    _sy_sy_corr(N, J, h, i, j, t; β = Inf) -> ComplexF64

`⟨σʸ_i(t) σʸ_j(0)⟩_β` for the OBC TFIM at inverse temperature `β`.
Wraps [`_sy_sy_corr_from_cached`](@ref) with a fresh Majorana
diagonalisation and time-evolution matrix.
"""
function _sy_sy_corr(N::Int, J::Float64, h::Float64, i::Int, j::Int, t::Real; β::Real=Inf)
    (1 ≤ i ≤ N && 1 ≤ j ≤ N) ||
        throw(ArgumentError("YY correlator: site indices out of range (i=$i, j=$j, N=$N)"))
    hmat = _majorana_ham(N, J, h)
    Σ = _majorana_thermal_covariance(hmat, β)
    R = _majorana_evolution(hmat, t)
    return _sy_sy_corr_from_cached(Σ, R, i, j)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Diagonal value: ⟨(σʸ_i)²⟩ = 1  (Pauli identity)
# ═══════════════════════════════════════════════════════════════════════════════
#
# Wick / Pfaffian on the i = j case would compute the squared product of
# 2(2i-1) Majoranas, which evaluates to `(-i)^{2(i-1)} · (-1) = -1 · -1 · …`
# but we know the answer analytically: `σʸ_i² = I`.  Returning `1.0`
# directly avoids a degenerate Pfaffian from the perfectly-correlated
# t=0 contraction.

# ═══════════════════════════════════════════════════════════════════════════════
# fetch — YY static / connected / dynamic at OBC
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::TFIM, ::SpinCorrelation{:y,:y}, bc::OBC;
          beta=Inf, i, j) -> Float64

Static thermal correlator `⟨σʸ_i σʸ_j⟩_β` on the OBC TFIM.  Equivalent
to the `t = 0` slice of [`DynamicalCorrelation`](@ref)`(:y, :y)`,
returned as `Float64` (the imaginary part is round-off only at t = 0).
"""
function fetch(
    model::TFIM,
    ::SpinCorrelation{:y,:y},
    bc::OBC;
    beta::Real=Inf,
    i::Int,
    j::Int,
    kwargs...,
)
    N = _bc_size(bc, kwargs)
    1 ≤ i ≤ N || throw(ArgumentError("YY static: site i = $i out of 1..$N"))
    1 ≤ j ≤ N || throw(ArgumentError("YY static: site j = $j out of 1..$N"))
    if i == j
        return 1.0  # σʸ² = I exactly
    end
    return real(_sy_sy_corr(N, model.J, model.h, i, j, 0.0; β=beta))
end

"""
    fetch(model::TFIM, ::ConnectedSpinCorrelation{:y,:y}, bc::OBC;
          beta=Inf, i, j) -> Float64

Connected thermal correlator `⟨σʸ_i σʸ_j⟩_β − ⟨σʸ_i⟩ ⟨σʸ_j⟩`.  Since
`⟨σʸ⟩ = 0` in any Gaussian state of the TFIM (odd-Majorana product),
the connected and static values coincide off-diagonal; the diagonal
returns `1 - 0² = 1`.
"""
function fetch(
    model::TFIM,
    ::ConnectedSpinCorrelation{:y,:y},
    bc::OBC;
    beta::Real=Inf,
    i::Int,
    j::Int,
    kwargs...,
)
    N = _bc_size(bc, kwargs)
    return fetch(model, SpinCorrelation(:y, :y), bc; beta=beta, i=i, j=j, N=N)
end

"""
    fetch(model::TFIM, ::DynamicalCorrelation{(:y, :y)}, bc::OBC;
          beta=Inf, i, j, t) -> ComplexF64

Real-time correlator `⟨σʸ_i(t) σʸ_j(0)⟩_β`.  Returns `ComplexF64`; the
imaginary part is non-zero in general (Re part is even in t, Im part is
odd in t — see the time-domain identity tests in
`test/identities/test_TFIM_dynamic_symmetries.jl`).
"""
function fetch(
    model::TFIM,
    ::DynamicalCorrelation{(:y, :y)},
    bc::OBC;
    i::Int,
    j::Int,
    t::Float64,
    beta::Real=Inf,
    kwargs...,
)
    N = _bc_size(bc, kwargs)
    return _sy_sy_corr(N, model.J, model.h, i, j, t; β=beta)
end

# ═══════════════════════════════════════════════════════════════════════════════
# MagnetizationY OBC — identically zero by Wick / parity
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::TFIM, ::MagnetizationY, bc::OBC; beta) -> Float64

Per-site bulk magnetisation `⟨Σᵢ σʸᵢ⟩_β / N` of the OBC TFIM.
Identically zero in any Gaussian state because `σʸ_i` reduces to an
odd product of Majoranas.  Returned as exact `0.0` so callers can use
it as a deterministic baseline against random-sample estimators that
fluctuate around zero.
"""
function fetch(::TFIM, ::MagnetizationY, ::OBC; beta::Real, kwargs...)
    return 0.0
end

# ═══════════════════════════════════════════════════════════════════════════════
# SusceptibilityYY OBC — β·Var(M_y)/N via Wick / Pfaffian
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::TFIM, ::SusceptibilityYY, bc::OBC; beta) -> Float64

Per-site equal-time fluctuation `χ_yy(β) = (β / N) · Σ_{i,j} ⟨σʸ_i σʸ_j⟩_β`
of the OBC TFIM.  ⟨σʸ⟩ = 0 in this Gaussian state, so the variance
form simplifies to `(β/N) · ⟨M_y²⟩`.

Implementation: per pair `(i, j)` evaluate the Pfaffian of the static
Majorana Wick matrix.  Diagonal contribution `⟨(σʸᵢ)²⟩ = 1` (Pauli
identity) is added directly, off-diagonal twice (symmetric).  Cost is
`O(N² · M³)` with `M = 2 max(i, j) − 1`; same scaling as
`SusceptibilityXX OBC` and `_xx_uniform_susceptibility` in
`TFIM_thermal.jl`.
"""
function fetch(model::TFIM, ::SusceptibilityYY, bc::OBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    hmat = _majorana_ham(N, model.J, model.h)
    Σ = _majorana_thermal_covariance(hmat, beta)
    R = _majorana_evolution(hmat, 0.0)  # = identity, but kept for cache shape
    s = N * 1.0  # diagonal: Σᵢ ⟨(σʸᵢ)²⟩ = N
    @inbounds for i in 1:N, j in (i + 1):N
        cij = real(_sy_sy_corr_from_cached(Σ, R, i, j))
        s += 2 * cij
    end
    return beta * s / N
end
