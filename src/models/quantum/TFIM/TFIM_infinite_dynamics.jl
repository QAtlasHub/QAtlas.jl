# ─────────────────────────────────────────────────────────────────────────────
# Transverse Field Ising Model — Infinite-volume dynamic / kinematic
# observables.
#
# This module exposes three quantities that complete TFIM's Tier-2 dynamic
# coverage on `Infinite()`:
#
#   1. `tfim_quasiparticle_dispersion(model, k)`  — closed-form Bogoliubov
#      single-quasiparticle band Λ(k) = 2 √(J² + h² − 2 J h cos k).  Public
#      helper, exported from `QAtlas`.
#
#   2. `tfim_two_spinon_dos(model, ω; q_total = 0)`  — two-spinon density of
#      states ρ_2(ω; q_total) at fixed total momentum, computed by the
#      standard 1/|f'(k*)| sum at the saddle points where
#      Λ(k) + Λ(q_total − k) = ω.  Returns 0 outside the 2-spinon
#      continuum.  Public helper, exported from `QAtlas`.
#
#   3. `fetch(model, ZZStructureFactor(), Infinite(); β, q, ω, …)` — dynamic
#      longitudinal structure factor S_zz(q, ω; β), realised as the
#      large-N OBC proxy of the time- and space-Fourier transform of
#      ⟨σᶻ_i(t) σᶻ_j(0)⟩_β.  σᶻ in TFIM is non-local after Jordan–Wigner
#      so a closed-form free-fermion form-factor expression is intricate;
#      QAtlas defers it to the OBC Pfaffian routine in `TFIM_dynamics.jl`,
#      reusing the cached Majorana covariance / evolution matrices.
#
# The static branch (`ω === nothing`) of `ZZStructureFactor, Infinite` is
# already provided in `TFIM_zaxis.jl`.  This file overrides that method
# with a router that dispatches on `ω` so both branches coexist.
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Single-quasiparticle dispersion
# ═══════════════════════════════════════════════════════════════════════════════

"""
    tfim_quasiparticle_dispersion(model::TFIM, k::Real) -> Float64

Single-quasiparticle Bogoliubov dispersion of the infinite TFIM,

    Λ(k) = 2 √(J² + h² − 2 J h cos k),

at momentum `k ∈ [0, π]`.  Useful for plotting band structure and as
the kinematic input to two-spinon density-of-states / structure factor
calculations.

Special values:
* `Λ(0) = 2 |J − h|`     — equal to the gap `Δ` in either phase.
* `Λ(π) = 2 (J + h)`      — top of the band.
* min over `k ∈ [0, π]` is the mass gap [`MassGap`](@ref) at `Infinite`.
"""
function tfim_quasiparticle_dispersion(model::TFIM, k::Real)
    J = float(model.J)
    h = float(model.h)
    return 2 * sqrt(J^2 + h^2 - 2 * J * h * cos(k))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Two-spinon density of states
# ═══════════════════════════════════════════════════════════════════════════════

"""
    tfim_two_spinon_dos(model::TFIM, ω::Real; q_total::Real = 0.0) -> Float64

Two-spinon density of states at total momentum `q_total` and frequency
`ω` for the infinite TFIM:

    ρ_2(ω; q_total) = (1/π) ∫₀^π dk  δ(ω − Λ(k) − Λ(q_total − k)).

For `q_total = 0` the two-spinon continuum is supported on
`[2 Δ, Λ(0) + Λ(π)] = [2 |J − h|, 2 (J + h)]`; outside this window
the routine returns `0.0`.

Computation: numerical root-finding for `k* ∈ (0, π)` where
`f(k) := Λ(k) + Λ(q_total − k) = ω`, then

    ρ_2(ω; q_total) = (1/π) Σ_{k*} 1 / |f'(k*)|.

Roots are isolated by a brute-force scan with `Nscan = 4096` samples
followed by 50 bisection refinements per bracket — sufficient for
~12-digit precision in `k*` and ~`1/Nscan` precision in counting roots
near van Hove singularities (where `|f'|` vanishes and ρ_2 diverges
integrably).  Right at a van Hove point the returned value is finite
because the bisected `k*` is offset from the exact saddle.
"""
function tfim_two_spinon_dos(model::TFIM, ω::Real; q_total::Real=0.0)
    J = float(model.J)
    h = float(model.h)
    Λ(k) = 2 * sqrt(J^2 + h^2 - 2 * J * h * cos(k))
    Λp(k) = (2 * J * h * sin(k)) / sqrt(J^2 + h^2 - 2 * J * h * cos(k))
    f(k) = Λ(k) + Λ(q_total - k)
    fp(k) = Λp(k) - Λp(q_total - k)

    Nscan = 4096
    ks = range(1e-8, π - 1e-8; length=Nscan)
    fs = [f(k) for k in ks]

    # Quick reject: ω outside the scanned support.
    if all(fs .> ω + 1e-14) || all(fs .< ω - 1e-14)
        return 0.0
    end

    ρ = 0.0
    for i in 1:(Nscan - 1)
        if (fs[i] - ω) * (fs[i + 1] - ω) < 0
            a = ks[i]
            b = ks[i + 1]
            fa = fs[i]
            for _ in 1:50
                m = (a + b) / 2
                fm = f(m)
                if abs(fm - ω) < 1e-12
                    a = b = m
                    break
                end
                if (fa - ω) * (fm - ω) < 0
                    b = m
                else
                    a = m
                    fa = fm
                end
            end
            kstar = (a + b) / 2
            denom = abs(fp(kstar))
            # Skip near-van-Hove zeros of f' to avoid spurious blow-up; the
            # contribution at exact saddle is integrable but pointwise infinite.
            if denom > 1e-10
                ρ += 1 / denom
            end
        end
    end
    return ρ / π
end

# ═══════════════════════════════════════════════════════════════════════════════
# Dynamic ZZ structure factor at Infinite — OBC large-N proxy
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _tfim_zz_structure_factor_dynamic_proxy(J, h, β, q, ω, N_proxy, t_max, dt) -> Float64

OBC large-N proxy implementation of the dynamic longitudinal structure
factor S_zz(q, ω; β).  Used by the `Infinite()` `fetch` router; the
heavy lifting is performed here so the public method stays a thin
wrapper.

The Fourier integral is replaced by a finite Riemann sum on
`t ∈ [-t_max, t_max]` with spacing `dt`, and `i, j` run over the
central bulk window `[N/4, 3N/4]` of an `N_proxy`-site OBC chain to
suppress boundary contamination.  The Majorana Hamiltonian and
thermal covariance are computed once; the evolution matrix `R(t)` is
recomputed at every time step (single 2N × 2N `expm`).
"""
function _tfim_zz_structure_factor_dynamic_proxy(
    J::Real, h::Real, β::Real, q::Real, ω::Real, N_proxy::Int, t_max::Real, dt::Real
)
    Jf = Float64(J)
    hf = Float64(h)
    N = N_proxy
    bulk_lo = max(1, N ÷ 4)
    bulk_hi = min(N, 3 * N ÷ 4)
    Nb = bulk_hi - bulk_lo + 1

    hmat = _majorana_ham(N, Jf, hf)
    Σ = _majorana_thermal_covariance(hmat, β)

    ts = collect((-t_max):dt:t_max)
    S = 0.0 + 0.0im

    for t in ts
        R = _majorana_evolution(hmat, t)
        cij_t = 0.0 + 0.0im
        for i in bulk_lo:bulk_hi, j in bulk_lo:bulk_hi
            c = _sz_sz_corr_from_cached(Σ, R, i, j)
            cij_t += c * exp(-im * q * (i - j))
        end
        cij_t /= Nb
        S += exp(im * ω * t) * cij_t * dt
    end
    return real(S)
end

"""
    fetch(model::TFIM, ::ZZStructureFactor, ::Infinite;
          beta::Real, q::Real, ω::Union{Real,Nothing} = nothing,
          N_proxy::Int = 64, t_max::Real = 20.0, dt::Real = 0.1, kwargs...)
        -> Float64

Longitudinal structure factor `S_zz(q, [ω,] β)` of the infinite TFIM.
This router dispatches on `ω`:

* `ω === nothing` (default) → static `S_zz(q, β)`, computed by the
  large-N OBC proxy `_zz_static_structure_factor` (default
  `N_proxy = 80`; same path as defined in `TFIM_zaxis.jl`).

* `ω::Real` → dynamic `S_zz(q, ω; β)`, computed as the OBC large-N
  proxy of the time- and space-Fourier transform of
  `⟨σᶻ_i(t) σᶻ_j(0)⟩_β`,

      S_zz(q, ω; β) = ∫dt e^{iωt} · (1/N_b) Σ_{i,j ∈ bulk}
                            e^{-iq(i-j)} ⟨σᶻ_i(t) σᶻ_j(0)⟩_β,

  with `t ∈ [-t_max, t_max]` discretised at spacing `dt` and `(i, j)`
  restricted to the central bulk window `[N/4, 3N/4]` of an
  `N_proxy`-site OBC chain.

  Default `N_proxy = 64`, `t_max = 20.0`, `dt = 0.1` is a balance of
  precision and cost; the dominant errors are
  (a) ω-resolution `~ π/t_max ≈ 0.157`,
  (b) UV cutoff `~ π/dt ≈ 31.4`,
  (c) finite-size finite-bulk corrections `~ exp(−(N_proxy − 4 ξ)/ξ)`.
  Raise the appropriate parameter to tighten any of these.

Performance: `O(|ts| · N_b² · M³)` Pfaffians per `(q, ω)` point, where
`M ≈ 2 (i + j) − 2`.  At default settings this is ~1 sec/point on a
single core after the Majorana eigendecomposition is amortised.
Multiple `(q, ω)` points should be batched by writing a custom loop
that reuses `Σ` and the per-`t` evolution matrices `R(t)`.

Static path remains the recommended one for any equilibrium sum-rule
work; the dynamic path is intended primarily as a benchmark for
TPQMPS / DMRG dynamic structure factor reference values.
"""
function fetch(
    model::TFIM,
    ::ZZStructureFactor,
    ::Infinite;
    beta::Real,
    q::Real,
    ω::Union{Real,Nothing}=nothing,
    N_proxy::Int=ω === nothing ? 80 : 64,
    t_max::Real=20.0,
    dt::Real=0.1,
    kwargs...,
)
    if ω === nothing
        # Static branch — preserve the behaviour of the method previously
        # defined in TFIM_zaxis.jl by routing to the same proxy.
        return _zz_static_structure_factor(N_proxy, model.J, model.h, beta, q)
    end
    return _tfim_zz_structure_factor_dynamic_proxy(
        model.J, model.h, beta, q, ω, N_proxy, t_max, dt
    )
end
