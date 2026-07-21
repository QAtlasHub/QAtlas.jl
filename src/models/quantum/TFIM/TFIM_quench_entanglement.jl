# ─────────────────────────────────────────────────────────────────────────────
# Transverse Field Ising Model — post-quench entanglement entropy via the
# Calabrese–Cardy quasi-particle picture (J. Stat. Mech. P04010 (2005)).
#
# Setup.  Prepare the system in the ground state |Ψ_0⟩ of an *initial*
# TFIM Hamiltonian H_0 = TFIM(J_0, h_0).  At t = 0 quench the parameters
# to a *final* Hamiltonian H_f = TFIM(J_f, h_f) and let the state evolve
# unitarily,
#
#   |Ψ(t)⟩ = exp(-i H_f t) |Ψ_0⟩.
#
# For a contiguous block A = {site 1, …, site ℓ} of the OBC chain, the
# Gaussian time-evolved state has reduced density matrix
#
#   ρ_A(t) = Tr_B |Ψ(t)⟩⟨Ψ(t)|,
#
# whose von Neumann entropy is computed by Peschel's correlation-matrix
# trick once the time-evolved Majorana covariance Σ(t) is known.
#
# Implementation.  In the JW + Majorana representation used throughout
# `TFIM_dynamics.jl` (see that file's header for sign conventions),
#
#   1. Σ_0 = `_majorana_covariance_gs(h_mat_0)` — initial GS covariance
#      built from H_0's Majorana matrix `h_mat_0 = _majorana_ham(N, J_0, h_0)`.
#
#   2. R(t) = `_majorana_evolution(h_mat_f, t)` ∈ SO(2N) is the Heisenberg
#      Majorana evolution under H_f, where
#      `h_mat_f = _majorana_ham(N, J_f, h_f)`.
#
#   3. Σ(t) = R(t) Σ_0 R(t)ᵀ — congruence transformation of the real
#      antisymmetric covariance.
#
#   4. Σ_A(t) = Σ(t)[1:2ℓ, 1:2ℓ] — restriction to the 2ℓ Majoranas on
#      sites 1..ℓ.
#
#   5. The Hermitian matrix `i Σ_A(t)` has eigenvalues in ± pairs
#      {±ν_k(t)}_{k=1..ℓ} with 0 ≤ ν_k(t) ≤ 1; the per-mode entropy
#      `_peschel_mode_entropy(ν_k)` (defined in `TFIM_entanglement.jl`)
#      gives
#
#        S(ℓ, t) = Σ_{k=1}^{ℓ} s(ν_k(t)).
#
# The same JW-factorisation argument that justifies the equilibrium
# `S^{(f)}_A = S^{(s)}_A` (Fagotti–Calabrese 2010) carries over to the
# time-evolved state because contiguous bipartitions preserve the JW
# string commutation that pulls the parity factor out of ρ_A(t).
#
# Cost.  O(N³) for `_majorana_evolution` (full matrix exponential) plus
# O(ℓ³) for the Hermitian eigendecomposition of `i Σ_A(t)`.  Reusing R(t)
# across ℓ would amortise to O(ℓ³); we don't expose that path yet
# because the existing per-call API is already cheap (~10 ms for N = 64).
#
# Universal scaling (Calabrese–Cardy 2005).  For a quench of a free
# fermion chain into the gapless TFIM, the entropy grows linearly,
#
#   S(ℓ, t) ≈ (c/3) v_E t           for t < ℓ / (2 v_E),
#
# and saturates at a volume-law value `(c/3) v_E (ℓ / 2 v_E) = (c/6) ℓ`
# (modulo non-universal constants set by H_0).  The crossover scale
# t* = ℓ / (2 v_E) is the time at which a quasi-particle pair created
# at the centre of A first reaches both endpoints.
#
# References.
#   * P. Calabrese, J. Cardy, J. Stat. Mech. P04010 (2005) — quasi-
#     particle picture for entanglement spreading after a global quench.
#   * I. Peschel, J. Phys. A 36, L205 (2003) — correlation-matrix method.
#   * M. Fagotti, P. Calabrese, J. Stat. Mech. P04016 (2008) — explicit
#     free-fermion construction (Σ(t) congruence form).
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra: eigvals, Hermitian

"""
    fetch(model_f::TFIM, ::QuenchEntanglementEntropy, bc::OBC;
          initial::TFIM, ℓ::Int, t::Real, kwargs...) -> Float64

Post-quench von Neumann entanglement entropy of the first `ℓ` sites of
the N-site OBC TFIM.

Prepare the chain in the ground state of `initial::TFIM` and quench
instantly to `model_f::TFIM`; this method returns

    S(ℓ, t) = -Tr ρ_A(t) log ρ_A(t)

evaluated by Peschel's correlation-matrix method on the time-evolved
Majorana covariance Σ(t) = R(t) Σ_0 R(t)ᵀ.  See the file header
(`TFIM_quench_entanglement.jl`) for the full derivation and the
Calabrese–Cardy quasi-particle picture for the expected linear-growth
behaviour.

Required keyword arguments
- `initial::TFIM`         — initial-Hamiltonian model whose ground state
                            is the t = 0 state.
- `ℓ::Int`                — subsystem length, `1 ≤ ℓ ≤ N - 1`.
- `t::Real`               — post-quench time.

`N` is read from `OBC(N)` (or legacy `kwargs[:N]`).  At `t = 0` the
result coincides with the equilibrium `VonNeumannEntropy`
of the *initial* model — this is the back-compat sanity check exercised
in `test/standalone/test_tfim_quench_entanglement.jl`.

Cost: `O(N³)` from the matrix exponential plus `O(ℓ³)` from the Peschel
eigendecomposition.

References: Calabrese–Cardy J. Stat. Mech. P04010 (2005); Peschel J.
Phys. A 36, L205 (2003).
"""
function fetch(
    model_f::TFIM,
    ::QuenchEntanglementEntropy,
    bc::OBC;
    initial::TFIM,
    ℓ::Int,
    t::Real,
    kwargs...,
)
    N = _bc_size(bc, kwargs)
    1 ≤ ℓ ≤ N - 1 || throw(
        ArgumentError(
            "QuenchEntanglementEntropy: ℓ must satisfy 1 ≤ ℓ ≤ N - 1; got ℓ = $ℓ, N = $N.",
        ),
    )

    hmat_0 = _majorana_ham(N, initial.J, initial.h)
    hmat_f = _majorana_ham(N, model_f.J, model_f.h)

    Σ_0 = _majorana_covariance_gs(hmat_0)
    R = _majorana_evolution(hmat_f, t)

    Σ_t = R * Σ_0 * transpose(R)
    # Re-symmetrise to strict antisymmetry to cancel round-off drift —
    # the same hygiene `_majorana_covariance_gs` applies to its output.
    Σ_t = (Σ_t - transpose(Σ_t)) / 2

    Σ_A = Σ_t[1:(2ℓ), 1:(2ℓ)]
    λ = eigvals(Hermitian(im .* Σ_A))
    S = 0.0
    @inbounds for k in (ℓ + 1):(2ℓ)
        S += _peschel_mode_entropy(λ[k])
    end
    return S
end
