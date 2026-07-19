# ─────────────────────────────────────────────────────────────────────────────
# Transverse Field Ising Model — Fidelity susceptibility χ_F(h)
#
# Closed-form analytical implementation via Bogoliubov-de Gennes diagonalisation.
# The fidelity susceptibility is the curvature of the ground-state overlap as a
# function of the driving parameter, here taken to be the transverse field h:
#
#   χ_F(h) = lim_{δ → 0} -2 log|⟨ψ_0(h) | ψ_0(h + δ)⟩| / δ²
#          = Σ_{n ≠ 0}  |⟨n | ∂_h H | 0⟩|² / (E_n − E_0)²    (2nd-order PT form)
#
# For the TFIM   H = -J Σᵢ σᶻᵢ σᶻᵢ₊₁ - h Σᵢ σˣᵢ
# we have ∂_h H = -Σᵢ σˣᵢ.  After Jordan-Wigner + Bogoliubov, only 2-quasiparticle
# excitations couple to the ground state.  In momentum space (Infinite, PBC) each
# (k, -k) pair contributes a single 2-qp channel of energy 2 Λ_k, giving
#
#   χ_F(h) per pair = (∂_h θ_k)²
#
# where θ_k is the Bogoliubov rotation angle that diagonalises the
# 2×2 BdG block at momentum k.  The standard convention (Sachdev,
# "Quantum Phase Transitions" 2e §5.1; Damski PRB 87, 165101) is
#
#   tan(2 θ_k) = J sin k / (h - J cos k)                        (★)
#
# with the factor of 2 inside `tan` reflecting that θ_k mixes the
# pair (k, -k) modes (the rotation is by 2θ_k in the SU(2) Bogoliubov
# subspace).  Differentiating (★) with respect to h:
#
#   ∂_h(2 θ_k) = -J sin k / ε_k²,    ε_k = √(J² + h² - 2Jh cos k) = Λ_k / 2,
#   ⇒ ∂_h θ_k = -J sin k / (2 ε_k²).
#
# Per (k, -k) pair:
#
#   χ_F(h)/pair = (∂_h θ_k)² = (J sin k)² / (4 ε_k⁴),
#
# and integrating over k > 0 in the thermodynamic limit gives
#
#   χ_F / L = (1 / 8π) ∫₀^π (J sin k)² / ε_k⁴ dk
#           = (2 / π) ∫₀^π (J sin k)² / Λ_k(h)⁴ dk          (Λ_k = 2 ε_k).
#
# The closed forms follow from the residue identity
#   ∫₀^π sin² k / (J² + h² − 2Jh cos k)² dk = (π / (4 J² h²))[(J² + h²)/|J² − h²| − 1]
#
# which simplifies in each phase:
#   - ordered  (h < J):  χ_F / L = 1 / (16 (J² − h²))
#   - disordered (h > J): χ_F / L = J² / (16 h² (h² − J²))
#
# Both diverge linearly at criticality (`χ_F / L ∝ 1/|h − J|` as `h → J`),
# consistent with the leading singular scaling χ_F / L ~ |h − h_c|^{−1} for
# the 1D Ising universality class (`ν = 1, d = 1, z = 1`; singular exponent
# ν(2 + 2/ν − d) − 1 = 1).
#
# OBC:  build the BdG matrix H_BdG (2N × 2N) directly,
#
#         H_BdG = [ A   B ; -B  -A ],   A_ii = 2h, A_{i,i±1} = -J,
#                                       B_{i,i+1} = +J, B_{i+1,i} = -J,
#
#       and diagonalise it as one Hermitian eigenproblem.  Particle-hole
#       symmetry guarantees eigenvalues come in ±Λ_n pairs.  Take the N
#       eigenvectors with positive eigenvalues, written as (u_n, v_n) ∈ ℝ^N
#       × ℝ^N — these are the standard Bogoliubov amplitudes
#
#         η_n = Σⱼ (u_{n,j} c_j + v_{n,j} c_j†).
#
#       The fidelity susceptibility is then
#
#         χ_F = Σ_{p < q} 4 X_{pq}² / (Λ_p + Λ_q)²,
#         X_{pq} = Σⱼ (u_{q,j} v_{p,j} − u_{p,j} v_{q,j}),
#
#       i.e. the antisymmetric 2-qp matrix element of (1/2) ∂_h H = Σᵢ c†ᵢ cᵢ
#       + const, the only piece coupling to (η_p†η_q†|0⟩) channels.
#
#       Using the full 2N × 2N diagonalisation (rather than the squared
#       N × N Lieb–Schultz–Mattis form (A−B)(A+B)) sidesteps the zero-mode
#       division `ψ_n = (A+B)φ_n / Λ_n` which becomes ill-conditioned in
#       the ordered phase where the symmetry-breaking edge mode has
#       Λ ≈ 10⁻¹⁵.
#
# References:
#   - Gu, "Fidelity approach to quantum phase transitions",
#     [Gu2010](@cite) — review (Eq. 5.55 for TFIM).
#   - Damski, "Fidelity approach to quantum phase transitions", PRB 87, 165101
#     (2013) — closed form for TFIM (his Eq. 23–24, in the convention with
#     J = 1 and h ≡ g; differs from ours by an overall factor of 2 due to a
#     different bookkeeping of (k, -k) pairs and the factor in ∂_h H).
#
# Convention chosen here: χ_F is the *quantum information geometric* fidelity
# susceptibility per the Zanardi-Paunkovic / De Grandi definition,
#   χ_F = -∂²|⟨ψ(λ) | ψ(λ + δλ)⟩|² / ∂(δλ)² / 2  evaluated at δλ = 0.
# For TFIM at Infinite this gives χ_F / L = 1/(16(J²−h²)) for h<J and
# χ_F / L = J²/(16 h²(h²−J²)) for h>J.
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra: Symmetric, eigen
using QuadGK: quadgk

# ═══════════════════════════════════════════════════════════════════════════════
# Internal helpers
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _tfim_bdg_full(N, J, h) -> (Λ::Vector{Float64}, U::Matrix{Float64}, V::Matrix{Float64})

Full diagonalisation of the 2N × 2N TFIM OBC BdG Hamiltonian.  Returns the
N positive quasiparticle energies `Λ[n]` and the Bogoliubov amplitudes
`U[i, n]`, `V[i, n]` (site-index `i`, mode-index `n`) defining

    η_n = Σᵢ U[i, n] cᵢ + V[i, n] cᵢ†,

normalised by `Σᵢ (U[i, n]² + V[i, n]²) = 1`.

The 2N × 2N BdG matrix is `H_BdG = [A B; -B -A]` with `A` symmetric
(transverse field + hopping) and `B` antisymmetric (pairing); its
spectrum is symmetric about zero by particle-hole symmetry.  We pick
the N eigenvectors with positive eigenvalues; for each such
eigenvector `(x, y) ∈ ℝ²ᴺ`, we identify `U[:, n] = x`, `V[:, n] = y`.

Diagonalising the full 2N × 2N matrix (rather than the squared
N × N system `(A−B)(A+B) φ = Λ² φ`) is robust to the near-zero
edge mode in the ordered phase (`h < J`), which becomes
`Λ ≈ 10⁻¹⁵` and renders the `ψ = (A+B)φ / Λ` reconstruction
ill-conditioned at certain `N`.
"""
function _tfim_bdg_full(N::Int, J::Float64, h::Float64)
    A = zeros(Float64, N, N)
    @inbounds for i in 1:N
        A[i, i] = 2h
    end
    @inbounds for i in 1:(N - 1)
        A[i, i + 1] = -J
        A[i + 1, i] = -J
    end

    B = zeros(Float64, N, N)
    @inbounds for i in 1:(N - 1)
        B[i, i + 1] = J
        B[i + 1, i] = -J
    end

    HBdG = [A B; -B -A]
    F = eigen(Symmetric((HBdG + HBdG') / 2))
    # Eigenvalues come in ±Λ pairs; take the upper N (positive ones).
    # `eigen(Symmetric, …)` returns eigenvalues sorted ascending.
    pos_idx = (N + 1):(2N)
    Λ = F.values[pos_idx]
    # Numerical noise: floor tiny negatives at zero (ordered-phase edge mode).
    Λ = max.(Λ, 0.0)
    W = F.vectors[:, pos_idx]  # 2N × N
    U = W[1:N, :]               # particle component
    V = W[(N + 1):(2N), :]      # hole component
    return Λ, U, V
end

"""
    _tfim_chi_F_obc(N, J, h) -> Float64

Total fidelity susceptibility (NOT per site) of the OBC TFIM with `N`
sites, with respect to the transverse field `h`.

Formula (2nd-order perturbation theory on the BdG diagonalisation):

    χ_F = Σ_{p < q} 4 X_{pq}² / (Λ_p + Λ_q)²,
    X_{pq} = Σⱼ [U_{q,j} V_{p,j} − U_{p,j} V_{q,j}],

where `U[j, n]`, `V[j, n]` are the Bogoliubov particle/hole amplitudes
returned by `_tfim_bdg_full`.

Cost: O(N³) eigendecomposition + O(N³) for the X matrix + O(N²)
summation.
"""
function _tfim_chi_F_obc(N::Int, J::Float64, h::Float64)
    Λ, U, V = _tfim_bdg_full(N, J, h)
    # X[p, q] = Σⱼ (U[j,q] V[j,p] - U[j,p] V[j,q]) = (V' U)[p,q] - (U' V)[p,q]
    VU = V' * U
    X = VU .- VU'  # antisymmetric

    χ = 0.0
    @inbounds for p in 1:N, q in (p + 1):N
        denom = Λ[p] + Λ[q]
        if denom > 1e-12
            χ += 4 * X[p, q]^2 / denom^2
        end
    end
    return χ
end

"""
    _tfim_chi_F_infinite(J, h; rtol=1e-10) -> Float64

Per-site fidelity susceptibility χ_F / L for the infinite TFIM with
respect to the transverse field h, computed by Gauss–Kronrod quadrature
of the closed-form Bogoliubov-vacuum overlap integral

    χ_F / L = (1 / 8π) ∫₀^π (J sin k)² / ε_k(h)⁴ dk,
    ε_k(h)  = √(J² + h² − 2 J h cos k) = Λ_k(h) / 2.

Equivalent to `(2/π) ∫₀^π (J sin k)² / Λ_k(h)⁴ dk` if expressed in
terms of the doubled BdG dispersion `Λ_k = 2 ε_k`.

The closed form `χ_F / L = J² / (8 (J² - h²)²)` (`h ≠ J`) follows by
residue integration; the routine still does numerical quadrature so
that small numerical errors propagate consistently with the OBC
counterpart and so that the implementation is robust to future
generalisations (e.g. χ_F with respect to J).

Throws `DomainError` at the critical point `h = J` (logarithmic
divergence of the integrand at `k = 0`).
"""
function _tfim_chi_F_infinite(J::Float64, h::Float64; rtol::Float64=1e-10)
    if abs(abs(h) - abs(J)) < 1e-14
        throw(
            DomainError(
                h,
                "Fidelity susceptibility diverges at the critical point |h| = J; " *
                "query at h = J ± ε instead.",
            ),
        )
    end
    integrand = k -> begin
        ε2 = J^2 + h^2 - 2 * J * h * cos(k)
        (J * sin(k))^2 / ε2^2
    end
    val, _ = quadgk(integrand, 0.0, π; rtol=rtol)
    return val / (8 * π)
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch dispatch
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::TFIM, ::FidelitySusceptibility, bc::OBC;
          per_site::Bool=false, kwargs...) -> Float64

Ground-state fidelity susceptibility χ_F of the OBC TFIM with `N = bc.N`
sites with respect to the transverse field `h`:

    χ_F(h) = Σ_{n ≠ 0} |⟨n | ∂_h H | 0⟩|² / (E_n - E_0)²,

evaluated in closed form via the Bogoliubov diagonalisation (no
numerical differentiation).  Cost O(N³).

`per_site=true` returns `χ_F / N`.

References: Gu, [Gu2010](@cite); Damski,
PRB 87, 165101 (2013).
"""
function fetch(
    model::TFIM, ::FidelitySusceptibility, bc::OBC; per_site::Bool=false, kwargs...
)
    N = _bc_size(bc, kwargs)
    χ = _tfim_chi_F_obc(N, model.J, model.h)
    return per_site ? χ / N : χ
end

"""
    fetch(model::TFIM, ::FidelitySusceptibility, ::Infinite;
          rtol::Float64=1e-10, kwargs...) -> Float64

Per-site fidelity susceptibility χ_F / L of the infinite TFIM with
respect to the transverse field `h`, computed by Gauss–Kronrod
quadrature of the closed-form Bogoliubov-vacuum overlap integral.

Closed-form values (h ≠ J):

    χ_F / L = 1 / (16 (J² − h²))            (ordered, h < J)
    χ_F / L = J² / (16 h² (h² − J²))         (disordered, h > J)

Both branches diverge as `1 / |J − h|` at the critical point
`|h| = J` — a `DomainError` is thrown if `||h| − |J||` is below
`1e-14`.

References: Gu, [Gu2010](@cite); Damski,
PRB 87, 165101 (2013).
"""
function fetch(
    model::TFIM, ::FidelitySusceptibility, ::Infinite; rtol::Float64=1e-10, kwargs...
)
    return _tfim_chi_F_infinite(model.J, model.h; rtol=rtol)
end
