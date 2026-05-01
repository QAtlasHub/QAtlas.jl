# ─────────────────────────────────────────────────────────────────────────────
# Classical 2D Ising — finite-temperature thermodynamic potentials
#
# Two paths share the same fetch surface:
#
#   * Finite Lx × Ly torus (PBC):  log Z = log tr(T^Lx) via the existing
#     `_ising_sq_transfer_matrix`.  Energy / SpecificHeat / Entropy fall
#     out by `ForwardDiff` of `log Z` with respect to `β` (the transfer
#     matrix is generic enough to carry `Dual` numbers).
#
#   * Thermodynamic limit (`Infinite()`):  Onsager 1944 closed form for
#     `-βf` per site,
#
#         -β f(K) = log(√2 · cosh(2K))
#                   + (1/(2π)) ∫₀^π log[(1 + √(1 − g² sin²q)) / 2] dq,
#         g(K)    = 2 sinh(2K) / cosh²(2K),     K = βJ.
#
#     `Energy(:per_site)` and `SpecificHeat` follow by automatic
#     differentiation of `_onsager_log_z_per_site` with respect to β.
#     The `c_v` second derivative diverges logarithmically at the
#     Onsager critical point `K_c = log(1+√2)/2` (i.e. `g = 1`), so test
#     calls stay away from that slice.
#
# Granularity convention follows the rest of QAtlas:
#   `Energy{:per_site}`  is native at every BC (per-site is the only
#   finite quantity in the thermodynamic limit, and per-site is the
#   natural transfer-matrix output for the finite torus too).
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra: tr
using QuadGK: quadgk

# ─── Central finite-difference helpers ────────────────────────────────
# QAtlas keeps `ForwardDiff` in `[extras]` (test-only); the closed-form
# energy / specific-heat below are derived from `log Z` by symmetric
# central differences whose truncation error sits at `O(δ²)` and balances
# `O(eps/δ)` round-off when `δ ≈ eps^{1/3} ≈ 6e-6` for the first
# derivative and `δ ≈ eps^{1/4} ≈ 1.2e-4` for the second.

@inline function _cd1(f, β::Real; δ::Real=1e-5)
    return (f(β + δ) - f(β - δ)) / (2δ)
end

@inline function _cd2(f, β::Real; δ::Real=1e-3)
    return (f(β + δ) - 2 * f(β) + f(β - δ)) / δ^2
end

# ═══════════════════════════════════════════════════════════════════════════════
# Granularity declaration
# ═══════════════════════════════════════════════════════════════════════════════

native_energy_granularity(::IsingSquare, ::PBC) = :per_site
native_energy_granularity(::IsingSquare, ::Infinite) = :per_site

# ═══════════════════════════════════════════════════════════════════════════════
# Onsager closed form (Infinite)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _onsager_log_z_per_site(K) -> Float64

Onsager (1944) thermodynamic-limit `log Z / N = -βf` per site of the
classical 2D Ising model on the square lattice, parameterised by the
reduced coupling `K = βJ`.  The bond-counting convention matches
`_ising_sq_transfer_matrix(Lx, Ly, β, J)` (each bond enumerated once
for `Lx, Ly ≥ 3`):

    -β f(K) = log(2) + (1/(2π)) ∫₀^π dφ · log[(A(φ) + √(A(φ)² − B²))/2]
    A(φ)    = cosh²(2K) − sinh(2K) · cos(φ),
    B       = sinh(2K).

The integrand is finite for all `K ≥ 0` (`A ≥ B` strictly except at the
critical point `K_c = ln(1+√2)/2`, where `A − B → 0` at `φ = 0`).
QuadGK's adaptive quadrature handles the integrable square-root edge
without intervention.

Reference: L. Onsager, Phys. Rev. 65, 117 (1944), Eq. (105).
"""
function _onsager_log_z_per_site(K::Real)
    B = sinh(2K)
    integrand = function (φ)
        A = cosh(2K)^2 - B * cos(φ)
        disc = A^2 - B^2
        # `disc` is non-negative analytically; clamp to 0 to absorb round-off
        # at criticality (φ → 0, K = K_c).
        return log((A + sqrt(max(disc, zero(K)))) / 2)
    end
    val, _ = quadgk(integrand, 0.0, π; rtol=1e-10)
    return log(2) + val / (2π)
end

"""
    fetch(m::IsingSquare, ::FreeEnergy, ::Infinite; beta, J=m.J) -> Float64

Per-site Helmholtz free energy `f(β) = -β⁻¹ · log Z / N` of the
classical 2D Ising model on the infinite square lattice (Onsager 1944).
"""
function fetch(m::IsingSquare, ::FreeEnergy, ::Infinite; beta::Real, J::Real=m.J)
    return -_onsager_log_z_per_site(beta * J) / beta
end

"""
    fetch(m::IsingSquare, ::Energy{:per_site}, ::Infinite; beta, J=m.J) -> Float64

Per-site thermal energy `ε(β) = -∂(log Z / N)/∂β` from the Onsager
closed form via `ForwardDiff`.  At low `T` (β → ∞) `ε → -2J` — the
ferromagnetic ground state has every spin aligned, contributing
`-J` per bond and 2 bonds per site.
"""
function fetch(m::IsingSquare, ::Energy{:per_site}, ::Infinite; beta::Real, J::Real=m.J)
    return -_cd1(b -> _onsager_log_z_per_site(b * J), beta)
end

"""
    fetch(m::IsingSquare, ::SpecificHeat, ::Infinite; beta, J=m.J) -> Float64

Per-site specific heat `c_v(β) = β² · ∂²(log Z / N)/∂β²` via
`ForwardDiff` (twice).  Diverges logarithmically at the Onsager
critical point `K_c = ln(1+√2)/2` (i.e. `T_c = 2J/ln(1+√2)`).  Caller
is responsible for staying off that slice; finite values everywhere
else.
"""
function fetch(m::IsingSquare, ::SpecificHeat, ::Infinite; beta::Real, J::Real=m.J)
    return beta^2 * _cd2(b -> _onsager_log_z_per_site(b * J), beta)
end

"""
    fetch(m::IsingSquare, ::ThermalEntropy, ::Infinite; beta, J=m.J) -> Float64

Per-site Gibbs entropy `s(β) = β · (ε − f)` from the Onsager
free-energy and energy paths.
"""
function fetch(m::IsingSquare, ::ThermalEntropy, ::Infinite; beta::Real, J::Real=m.J)
    ε = fetch(m, Energy(:per_site), Infinite(); beta=beta, J=J)
    f = fetch(m, FreeEnergy(), Infinite(); beta=beta, J=J)
    return beta * (ε - f)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Finite Lx × Ly torus (PBC) — transfer matrix + ForwardDiff
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _ising_sq_log_z(m, Lx, Ly, β, J) -> Real

`log tr(T^Lx)` for the `Lx × Ly` torus — the finite-N partition
function in log form, generic in `β` and `J` so `ForwardDiff` Duals
propagate through the transfer matrix.
"""
function _ising_sq_log_z(::IsingSquare, Lx::Integer, Ly::Integer, β::Real, J::Real)
    T = _ising_sq_transfer_matrix(Ly, β, J)
    return log(tr(T^Lx))
end

"""
    fetch(m::IsingSquare, ::FreeEnergy, ::PBC; beta, Lx=m.Lx, Ly=m.Ly, J=m.J) -> Float64

Per-site Helmholtz free energy `f(β) = -log Z / (β · Lx · Ly)` for the
`Lx × Ly` torus.  Builds the `2^Ly × 2^Ly` transfer matrix and
trace-cubes; cost is `O(2^{3 Ly})`, so practical Ly ≤ 10–12.
"""
function fetch(
    m::IsingSquare,
    ::FreeEnergy,
    ::PBC;
    beta::Real,
    Lx::Integer=m.Lx,
    Ly::Integer=m.Ly,
    J::Real=m.J,
)
    Lx > 0 && Ly > 0 || error("IsingSquare PBC FreeEnergy: positive Lx, Ly required.")
    return -_ising_sq_log_z(m, Lx, Ly, beta, J) / (beta * Lx * Ly)
end

"""
    fetch(m::IsingSquare, ::Energy{:per_site}, ::PBC; beta, Lx, Ly, J) -> Float64

Per-site thermal energy `ε(β) = -∂(log Z)/∂β / (Lx · Ly)` for the
finite torus, via `ForwardDiff` over the transfer-matrix `log Z`.
"""
function fetch(
    m::IsingSquare,
    ::Energy{:per_site},
    ::PBC;
    beta::Real,
    Lx::Integer=m.Lx,
    Ly::Integer=m.Ly,
    J::Real=m.J,
)
    Lx > 0 && Ly > 0 || error("IsingSquare PBC Energy: positive Lx, Ly required.")
    dlogZ = _cd1(b -> _ising_sq_log_z(m, Lx, Ly, b, J), beta)
    return -dlogZ / (Lx * Ly)
end

"""
    fetch(m::IsingSquare, ::SpecificHeat, ::PBC; beta, Lx, Ly, J) -> Float64

Per-site specific heat `c_v(β) = β² · Var(H) / (Lx · Ly)` for the
finite torus, via `ForwardDiff` (twice) on `log Z`.
"""
function fetch(
    m::IsingSquare,
    ::SpecificHeat,
    ::PBC;
    beta::Real,
    Lx::Integer=m.Lx,
    Ly::Integer=m.Ly,
    J::Real=m.J,
)
    Lx > 0 && Ly > 0 || error("IsingSquare PBC SpecificHeat: positive Lx, Ly required.")
    d2logZ = _cd2(b -> _ising_sq_log_z(m, Lx, Ly, b, J), beta)
    return beta^2 * d2logZ / (Lx * Ly)
end

"""
    fetch(m::IsingSquare, ::ThermalEntropy, ::PBC; beta, Lx, Ly, J) -> Float64

Per-site Gibbs entropy `s(β) = β · (ε − f)` for the finite torus.
"""
function fetch(
    m::IsingSquare,
    ::ThermalEntropy,
    ::PBC;
    beta::Real,
    Lx::Integer=m.Lx,
    Ly::Integer=m.Ly,
    J::Real=m.J,
)
    ε = fetch(m, Energy(:per_site), PBC(0); beta=beta, Lx=Lx, Ly=Ly, J=J)
    f = fetch(m, FreeEnergy(), PBC(0); beta=beta, Lx=Lx, Ly=Ly, J=J)
    return beta * (ε - f)
end
