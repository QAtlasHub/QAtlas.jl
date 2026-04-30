# ─────────────────────────────────────────────────────────────────────────────
# Transverse Field Ising Model — static (equal-time) σˣσˣ correlators (Tier 2)
#
# Companion to `TFIM_dynamics.jl`.  All the heavy lifting (Majorana
# Hamiltonian assembly, BdG-thermal Majorana 2-point function, Wick / Pfaffian
# evaluation) is already provided there; this file only adds equal-time
# (`t = 0`) `fetch` dispatches for `XXCorrelation{:static}` /
# `XXCorrelation{:connected}` and an `Infinite`-proxy variant that re-uses
# the OBC routine at a large N.
#
# Conventions:
#   - σˣ_k = -i γ_{2k-1} γ_{2k}, so the Pfaffian dimension is fixed at 4×4
#     (independent of |i - j|), c.f. the σᶻσᶻ case in TFIM_dynamics.jl.
#   - At equal time the Wick / Pfaffian evaluation produces a real number
#     up to floating-point round-off; we project to `Float64` via `real(...)`.
#   - At i = j we have ⟨σˣ_i² ⟩ = ⟨I⟩ = 1 exactly; the connected piece
#     reduces to the on-site σˣ variance, 1 − ⟨σˣ_i⟩².
#   - YY (`σʸσʸ`) is *not* implemented here — the JW string for σʸ leaves
#     residual γ-string operators between the two endpoints in the OBC
#     setting, requiring a different Majorana index list than σˣ; deferred
#     to a follow-up issue.
#
# `Infinite` proxy:
#   - As with `SusceptibilityZZ`/`ZZStructureFactor` at `Infinite()` (see
#     TFIM_zaxis.jl), the thermodynamic-limit value is delivered as the
#     OBC value at a user-tunable proxy size `N_proxy` (default 80).  The
#     caller is responsible for supplying bulk-friendly indices i, j
#     (e.g. centred around `N_proxy / 2`).
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# OBC: ⟨σˣ_i σˣ_j⟩_β  (static)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::TFIM, ::XXCorrelation{:static}, bc::OBC;
          beta::Real = Inf, i::Int, j::Int, kwargs...) -> Float64

Static (equal-time) thermal transverse correlator
`⟨σˣ_i σˣ_j⟩_β` for the OBC TFIM at inverse temperature `beta`.
`beta = Inf` (the default) gives the ground-state value.

Implementation: re-uses `_sx_sx_corr(N, J, h, i, j, 0.0; β=beta)` from
`TFIM_dynamics.jl` — a 4×4 Pfaffian over the four Majoranas
`(γ_{2i-1}, γ_{2i}, γ_{2j-1}, γ_{2j})`.  The result is real at equal
time, so the imaginary residue (round-off) is dropped via `real(...)`.

At `i = j`, `⟨(σˣ)²⟩ = ⟨I⟩ = 1` regardless of (β, J, h).
"""
function fetch(
    model::TFIM,
    ::XXCorrelation{:static},
    bc::OBC;
    beta::Real=Inf,
    i::Int,
    j::Int,
    kwargs...,
)
    N = _bc_size(bc, kwargs)
    (1 ≤ i ≤ N) || throw(ArgumentError("i = $i out of range 1:$N"))
    (1 ≤ j ≤ N) || throw(ArgumentError("j = $j out of range 1:$N"))
    return real(_sx_sx_corr(N, model.J, model.h, i, j, 0.0; β=beta))
end

"""
    fetch(model::TFIM, ::XXCorrelation{:connected}, bc::OBC;
          beta::Real = Inf, i::Int, j::Int, kwargs...) -> Float64

Connected static thermal transverse correlator
`⟨σˣ_i σˣ_j⟩_β − ⟨σˣ_i⟩_β ⟨σˣ_j⟩_β` for the OBC TFIM.

Unlike `ZZCorrelation{:connected}` (where ⟨σᶻ⟩ = 0 by Z₂ on OBC and the
connected and bare correlators coincide), `⟨σˣ⟩_β ≠ 0` in general — σˣ
is the field-coupled order parameter and acquires a non-zero
expectation `Σ[2i-1, 2i]` from the BdG-thermal covariance.

For `i = j` the on-site σˣ variance simplifies to
`1 − ⟨σˣ_i⟩²` since `(σˣ)² = I`.

Internals: a single 2N × 2N BdG diagonalisation gives both the σˣ
expectation values (read off the covariance Σ) and, with the identity
evolution `R = I` at `t = 0`, the 4×4 Pfaffian for `⟨σˣ_i σˣ_j⟩_β`.
"""
function fetch(
    model::TFIM,
    ::XXCorrelation{:connected},
    bc::OBC;
    beta::Real=Inf,
    i::Int,
    j::Int,
    kwargs...,
)
    N = _bc_size(bc, kwargs)
    (1 ≤ i ≤ N) || throw(ArgumentError("i = $i out of range 1:$N"))
    (1 ≤ j ≤ N) || throw(ArgumentError("j = $j out of range 1:$N"))

    hmat = _majorana_ham(N, model.J, model.h)
    Σ = _majorana_thermal_covariance(hmat, beta)
    m_i = _sx_expect(Σ, i)

    if i == j
        # ⟨(σˣ)²⟩ = 1 exactly; connected = 1 − ⟨σˣ⟩².
        return 1.0 - m_i^2
    end

    R = _majorana_evolution(hmat, 0.0)  # = I, exact at t = 0
    c2 = real(_sx_sx_corr_from_cached(Σ, R, i, j))
    m_j = _sx_expect(Σ, j)
    return c2 - m_i * m_j
end

# ═══════════════════════════════════════════════════════════════════════════════
# Infinite-proxy variants (large-N OBC)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::TFIM, ::XXCorrelation{:static}, ::Infinite;
          beta::Real = Inf, i::Int, j::Int,
          N_proxy::Int = 80, kwargs...) -> Float64

Static `⟨σˣ_i σˣ_j⟩_β` in the thermodynamic limit, delivered as the
OBC value at proxy size `N_proxy` (default `80`) — same compromise as
`SusceptibilityZZ`/`ZZStructureFactor` at `Infinite()`, see
`TFIM_zaxis.jl` for the rationale.

The caller is responsible for picking bulk-friendly indices, e.g.
`N_proxy / 4 ≤ i, j ≤ 3 N_proxy / 4`.  At those interior sites the
boundary contamination decays exponentially with the distance to the
nearest edge in the gapped phase, and as `1/distance` at criticality.

Raise `N_proxy` if more accuracy is needed; the cost is the single
`2 N_proxy × 2 N_proxy` BdG diagonalisation.
"""
function fetch(
    model::TFIM,
    ::XXCorrelation{:static},
    ::Infinite;
    beta::Real=Inf,
    i::Int,
    j::Int,
    N_proxy::Int=80,
    kwargs...,
)
    return fetch(model, XXCorrelation{:static}(), OBC(N_proxy); beta=beta, i=i, j=j)
end

"""
    fetch(model::TFIM, ::XXCorrelation{:connected}, ::Infinite;
          beta::Real = Inf, i::Int, j::Int,
          N_proxy::Int = 80, kwargs...) -> Float64

Connected static `⟨σˣ_i σˣ_j⟩_β,c` in the thermodynamic limit, via the
same OBC large-N proxy as `XXCorrelation{:static}, Infinite()`.
"""
function fetch(
    model::TFIM,
    ::XXCorrelation{:connected},
    ::Infinite;
    beta::Real=Inf,
    i::Int,
    j::Int,
    N_proxy::Int=80,
    kwargs...,
)
    return fetch(model, XXCorrelation{:connected}(), OBC(N_proxy); beta=beta, i=i, j=j)
end
