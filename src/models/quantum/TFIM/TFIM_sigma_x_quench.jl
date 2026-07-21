# ─────────────────────────────────────────────────────────────────────────────
# Transverse Field Ising Model — sudden quench of the local transverse
# magnetisation ⟨σˣ_i⟩(t).
#
# Setup (sudden quench):
#   * Initial state  |ψ_0⟩ = ground state of H(h_0) = -J Σ σᶻ σᶻ - h_0 Σ σˣ.
#   * For t > 0 the system evolves with H_f = H(h_f).
#   * Observable  ⟨σˣ_i⟩(t) = ⟨ψ_0| e^{i H_f t} σˣ_i e^{-i H_f t} |ψ_0⟩.
#
# OBC (Majorana covariance evolution, exact at finite N):
#   1. Build the initial Majorana covariance Σ_0 of GS(H(h_0)) via
#      `_majorana_covariance_gs(_majorana_ham(N, J, h_0))`.
#   2. Heisenberg-evolve with the post-quench generator,
#      R(t) = exp(h_f · t) = `_majorana_evolution(_majorana_ham(N, J, h_f), t)`.
#   3. Time-evolved covariance  Σ(t) = R(t) Σ_0 R(t)^T  (the operator
#      identity γ_a(t) = R(t)_{ab} γ_b combined with R orthogonal).
#   4. Local magnetisation  ⟨σˣ_i⟩(t) = -i ⟨γ_{2i-1}(t) γ_{2i}(t)⟩
#                                     = Σ(t)[2i-1, 2i]  (real).
#
# Infinite (translationally invariant, closed-form k-integral):
#   In the thermodynamic limit each k-mode pair (k, -k) decouples.  With the
#   Bogoliubov angle θ_k(h) defined by
#
#       ε_k(h) = h − J cos k,   Δ_k = J sin k,
#       Λ_k(h) = 2 √(ε_k² + Δ_k²),
#       2 θ_k(h) = atan2( 2 Δ_k, 2 ε_k ),
#
#   the standard derivation (Barouch–McCoy–Dresden, PRA 2 (1970);
#   Calabrese–Essler–Fagotti, J. Stat. Mech. P07016 (2012)) gives
#
#       ⟨σˣ⟩(t) = (1/π) ∫₀^π dk [
#                      cos(2 θ_k^f) · cos(2 Δθ_k)
#                    + sin(2 θ_k^f) · sin(2 Δθ_k) · cos(2 Λ_k^f t) ]
#
#   with  Δθ_k ≡ θ_k(h_f) − θ_k(h_0).  Sanity:
#     • t = 0          → cos(2 θ_k^f − 2 Δθ_k) = cos(2 θ_k(h_0))
#                        ⇒ ⟨σˣ⟩(0) = ground-state value at h_0.
#     • h_0 = h_f      → Δθ_k = 0  ⇒  time-independent = GS at h_f.
#     • t → ∞ time average → diagonal ensemble (GGE) value
#                          ⟨σˣ⟩_GGE = (1/π) ∫ cos(2 θ_k^f) cos(2 Δθ_k) dk.
# ─────────────────────────────────────────────────────────────────────────────

using QuadGK: quadgk

# ═══════════════════════════════════════════════════════════════════════════════
# Internal: closed-form k-integral for the infinite quench
# ═══════════════════════════════════════════════════════════════════════════════

@inline function _tfim_two_theta(h::Real, J::Real, k::Real)
    # 2 θ_k(h) ∈ (-π, π].  Built via atan2 so the full Brillouin zone is
    # covered without quadrant ambiguity at h = J cos k (the gap-closing
    # locus on the BZ for h_c = J).  The factor of 2 in numerator/denominator
    # cancels but is kept for symmetry with Λ = 2 √((...)² + (...)²).
    return atan(2 * J * sin(k), 2 * (h - J * cos(k)))
end

@inline function _tfim_lambda(h::Real, J::Real, k::Real)
    return 2 * sqrt(J^2 + h^2 - 2 * J * h * cos(k))
end

"""
    _tfim_sigma_x_quench_infinite(J, h_0, h_f, t) -> Float64

Closed-form integral

    ⟨σˣ⟩(t) = (1/π) ∫₀^π dk [
                 cos(2 θ_k^f) · cos(2 Δθ_k)
               + sin(2 θ_k^f) · sin(2 Δθ_k) · cos(2 Λ_k^f t) ]

with Δθ_k = θ_k(h_f) − θ_k(h_0) for the infinite TFIM after a sudden
quench h_0 → h_f at fixed J.  Adaptive Gauss–Kronrod quadrature with
`rtol = 1e-12` (the integrand is bounded, smooth on (0, π) for any
non-critical h, and has a finite limit at the endpoints).
"""
function _tfim_sigma_x_quench_infinite(J::Real, h_0::Real, h_f::Real, t::Real)
    integrand = function (k)
        two_theta_f = _tfim_two_theta(h_f, J, k)
        two_theta_0 = _tfim_two_theta(h_0, J, k)
        d_two_theta = two_theta_f - two_theta_0
        Λ_f = _tfim_lambda(h_f, J, k)
        return cos(two_theta_f) * cos(d_two_theta) +
               sin(two_theta_f) * sin(d_two_theta) * cos(2 * Λ_f * t)
    end
    val, _ = quadgk(integrand, 0.0, π; rtol=1e-12)
    return val / π
end

# ═══════════════════════════════════════════════════════════════════════════════
# Internal: OBC Majorana covariance evolution — single (i, t) point
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _tfim_sigma_x_quench_obc(N, J, h_0, h_f, i, t) -> Float64

Exact `⟨σˣ_i⟩(t)` for the OBC TFIM after a sudden quench `h_0 → h_f`
at fixed `J`, computed from the Majorana covariance of the initial
ground state propagated by the post-quench Hamiltonian:

    Σ_0   = _majorana_covariance_gs(_majorana_ham(N, J, h_0))
    R(t)  = _majorana_evolution(_majorana_ham(N, J, h_f), t)
    Σ(t)  = R(t) Σ_0 R(t)^T
    ⟨σˣ_i⟩(t) = Σ(t)[2i-1, 2i]   (real).

A single `(i, t)` point costs one 2N × 2N eigendecomposition for Σ_0
plus one matrix exponential for R(t).  Sweeps in `i` should reuse the
returned matrix; sweeps in `t` should be hoisted to a custom loop that
caches `Σ_0` and recomputes only `R(t)`.
"""
function _tfim_sigma_x_quench_obc(N::Int, J::Real, h_0::Real, h_f::Real, i::Int, t::Real)
    (1 ≤ i ≤ N) || throw(ArgumentError("site index out of range: i = $i, N = $N"))
    h0_mat = _majorana_ham(N, Float64(J), Float64(h_0))
    hf_mat = _majorana_ham(N, Float64(J), Float64(h_f))
    Σ_0 = _majorana_covariance_gs(h0_mat)
    R = _majorana_evolution(hf_mat, t)
    Σt = R * Σ_0 * transpose(R)
    return Σt[2i - 1, 2i]
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch dispatch
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model_f::TFIM, ::QuenchLocalMagnetization{:x}, bc::OBC;
          initial::TFIM, i::Int, t::Real, kwargs...) -> Float64

Time-evolved local transverse magnetisation `⟨σˣ_i⟩(t)` of the OBC
TFIM after a sudden quench.

* `model_f` is the post-quench model (sets `h_f`, `J`).
* `initial` is the pre-quench TFIM whose ground state `|ψ_0⟩` is the
  initial state.  Both models must share the same `J`; mismatch raises
  an `ArgumentError` (the quench is *not* defined for a `J → J'` jump
  in the current implementation).
* `i ∈ 1:N`, `t ∈ ℝ`, `N` from `bc.N` (or kwargs).

Implementation: Majorana covariance evolution Σ(t) = R(t) Σ_0 R(t)^T
with Σ_0 = GS covariance under H(h_0) and R(t) = exp(h_f · t).
Cost per call: one 2N × 2N eigendecomposition + one matrix exponential.

Sanity checks (covered by `test/standalone/test_tfim_sigma_x_quench.jl`):
  * t = 0  → equilibrium `⟨σˣ_i⟩` of GS(h_0).
  * h_0 = h_f → time-independent (= equilibrium GS at h_0).
  * Large-N central-site → matches the `Infinite()` closed form.

References: Barouch–McCoy–Dresden, Phys. Rev. A **2** (1970) 1075;
Calabrese–Essler–Fagotti, J. Stat. Mech. **P07016** (2012).
"""
function fetch(
    model_f::TFIM,
    ::QuenchLocalMagnetization{:x},
    bc::OBC;
    initial::TFIM,
    i::Int,
    t::Real,
    kwargs...,
)
    N = _bc_size(bc, kwargs)
    initial.J == model_f.J || throw(
        ArgumentError(
            "TFIM σˣ quench: J mismatch (initial.J = $(initial.J), final.J = $(model_f.J)). " *
            "Quenches in J are not implemented; only h_0 → h_f is supported.",
        ),
    )
    return _tfim_sigma_x_quench_obc(N, model_f.J, initial.h, model_f.h, i, Float64(t))
end

"""
    fetch(model_f::TFIM, ::QuenchLocalMagnetization{:x}, ::Infinite;
          initial::TFIM, t::Real, kwargs...) -> Float64

Translationally-invariant `⟨σˣ⟩(t)` for the infinite TFIM after a
sudden quench from `H(initial.h)` to `H(model_f.h)` (`initial.J ==
model_f.J` required).  Closed-form k-integral:

    ⟨σˣ⟩(t) = (1/π) ∫₀^π dk [ cos(2 θ_k^f) cos(2 Δθ_k)
                            + sin(2 θ_k^f) sin(2 Δθ_k) cos(2 Λ_k^f t) ]

with Δθ_k ≡ θ_k(h_f) − θ_k(h_0), evaluated by adaptive Gauss–Kronrod
quadrature.

References: Barouch–McCoy–Dresden, PRA **2** (1970); Calabrese–Essler–
Fagotti, J. Stat. Mech. P07016 (2012).
"""
function fetch(
    model_f::TFIM,
    ::QuenchLocalMagnetization{:x},
    ::Infinite;
    initial::TFIM,
    t::Real,
    kwargs...,
)
    initial.J == model_f.J || throw(
        ArgumentError(
            "TFIM σˣ quench: J mismatch (initial.J = $(initial.J), final.J = $(model_f.J)). " *
            "Quenches in J are not implemented; only h_0 → h_f is supported.",
        ),
    )
    return _tfim_sigma_x_quench_infinite(model_f.J, initial.h, model_f.h, Float64(t))
end
