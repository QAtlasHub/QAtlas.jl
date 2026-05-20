# ─────────────────────────────────────────────────────────────────────────────
# Kitaev (2001) 1D p-wave superconducting wire — exact BdG solution.
#
# Hamiltonian (spinless fermions, real-space, lattice spacing a = 1):
#
#   H = -μ Σᵢ c†ᵢ cᵢ
#       - t Σᵢ (c†ᵢ c_{i+1} + h.c.)
#       + Δ Σᵢ (cᵢ c_{i+1} + h.c.)
#
# • μ — chemical potential (uniform onsite)
# • t — nearest-neighbour hopping (real)
# • Δ — p-wave pairing amplitude (real)
#
# After Jordan-Wigner / Bogoliubov-de Gennes diagonalisation the spectrum is
# free; on a PBC ring the closed-form dispersion is
#
#   E(k) = √( (2t cos k + μ)² + 4 Δ² sin² k ).
#
# Phase diagram (for Δ ≠ 0, t ≠ 0):
#
#   |μ| < 2|t|  → topological phase (winding ν = -1, two Majorana zero
#                 modes at OBC chain ends with hybridisation
#                 energy ~ e^{-N/ξ}).
#   |μ| > 2|t|  → trivial phase (ν = +1, gapped, no edge modes).
#   |μ| = 2|t|  → critical line (E(0) = 0 or E(π) = 0).
#
# The TFIM is the special case (μ = -2h, t = J, Δ = J): the BdG spectrum
# of `Kitaev1D(μ=-2h, t=J, Δ=J)` matches that of `TFIM(J=J, h=h)` at OBC,
# providing a built-in cross-check (`_kitaev1d_bdg_spectrum` is a strict
# generalisation of `_tfim_bdg_spectrum`).
#
# References:
#   - A. Yu. Kitaev, "Unpaired Majorana fermions in quantum wires",
#     Phys.-Usp. 44, 131 (2001).
#   - J. Alicea, "New directions for the pursuit of Majorana fermions in
#     solid state systems", Rep. Prog. Phys. 75, 076501 (2012).
#   - J. K. Asbóth, L. Oroszlány, A. Pályi, "A Short Course on Topological
#     Insulators", Lect. Notes Phys. 919 (2016) — Pfaffian invariant.
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Pauli σ (this file)
#   Observable:  Spin S = σ/2  (QAtlas-wide spin convention; see docs/src/conventions.md)

using LinearAlgebra: eigvals, Symmetric
using QuadGK: quadgk

"""
    Kitaev1D(; μ = 0.0, t = 1.0, Δ = 1.0) <: AbstractQAtlasModel

Kitaev (2001) one-dimensional p-wave superconducting wire,

```math
H = -\\mu \\sum_i c_i^{\\dagger} c_i
    - t   \\sum_i (c_i^{\\dagger} c_{i+1} + \\text{h.c.})
    + \\Delta \\sum_i (c_i c_{i+1} + \\text{h.c.}).
```

`μ` is the chemical potential, `t` the hopping, `Δ` the p-wave pairing.
`|μ| < 2|t|` is the topological phase (Majorana edge modes);
`|μ| > 2|t|` is the trivial phase; `|μ| = 2|t|` is the gapless critical
line.

This is the 1D Majorana wire and is **distinct** from the 2D
[`KitaevHoneycomb`](@ref) spin model.

The TFIM is a special case (`μ = -2h`, `t = J`, `Δ = J`); the BdG
spectrum of `Kitaev1D(μ=-2h, t=J, Δ=J)` agrees exactly with that of
`TFIM(J=J, h=h)` at OBC.
"""
struct Kitaev1D <: AbstractQAtlasModel
    μ::Float64
    t::Float64
    Δ::Float64
end
function Kitaev1D(; μ::Real=0.0, t::Real=1.0, Δ::Real=1.0)
    return Kitaev1D(Float64(μ), Float64(t), Float64(Δ))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Internal: BdG quasiparticle spectrum (OBC, finite N)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _kitaev1d_bdg_spectrum(N, μ, t, Δ) -> Vector{Float64}

Return the `N` non-negative BdG quasiparticle energies of the OBC
Kitaev1D chain with `N` sites and parameters `(μ, t, Δ)`.

The BdG matrix in Nambu basis is `H_BdG = [[A, B]; [-B, -A]]` (the same
real, antisymmetric 2N × 2N convention used by `_tfim_bdg_spectrum`),
with

    A_{ii}   = -μ
    A_{i,i+1} = A_{i+1,i} = -t
    B_{i,i+1} = +Δ,  B_{i+1,i} = -Δ.

Eigenvalues come in ± pairs; the function returns the non-negative half,
sorted ascending.  The BdG zero modes of the topological phase (lifted
by exponentially small N⁻¹ corrections) appear as the smallest entry.

This is a strict generalisation of `_tfim_bdg_spectrum`: at
`(μ, t, Δ) = (-2h, J, J)` the matrix coincides with the TFIM BdG matrix
and the returned spectrum matches `_tfim_bdg_spectrum(N, J, h)`.
"""
function _kitaev1d_bdg_spectrum(N::Int, μ::Float64, t::Float64, Δ::Float64)::Vector{Float64}
    A = zeros(N, N)
    @inbounds for i in 1:N
        A[i, i] = -μ
    end
    @inbounds for i in 1:(N - 1)
        A[i, i + 1] = -t
        A[i + 1, i] = -t
    end

    B = zeros(N, N)
    @inbounds for i in 1:(N - 1)
        B[i, i + 1] = Δ
        B[i + 1, i] = -Δ
    end

    H_bdg = [A B; -B -A]
    vals = eigvals(Symmetric(H_bdg))
    # Drop the negative half (eigenvalues come in ± pairs); keep the
    # non-negative branch even if one is numerically ~0 (Majorana zero
    # mode in the topological phase) — TFIM filters at 1e-10, but here
    # the topological zero mode is *the* observable we need to expose.
    pos = filter(v -> v >= -1e-12, vals)
    sort!(pos)
    return pos[1:N]
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy granularity convention (see src/core/quantities.jl)
# ═══════════════════════════════════════════════════════════════════════════════

native_energy_granularity(::Kitaev1D, ::Infinite) = :per_site

# ═══════════════════════════════════════════════════════════════════════════════
# ExactSpectrum (OBC) — full BdG quasiparticle spectrum
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::Kitaev1D, ::ExactSpectrum, bc::OBC; N::Int) -> Vector{Float64}

Return the `N` non-negative BdG quasiparticle energies of the OBC
Kitaev1D chain, sorted ascending.

`N` is read from `bc.N` (`OBC(N)` / `OBC(; N)`) or, as a legacy
fallback, from the `N` keyword argument.

In the topological phase (`|μ| < 2|t|`, `Δ ≠ 0`) the lowest entry is the
exponentially-small Majorana edge-mode hybridisation energy
`~ e^{-N/ξ}`.  In the trivial phase the lowest entry approaches the bulk
gap `min(|2t + μ|, |2t - μ|)`.
"""
function fetch(model::Kitaev1D, ::ExactSpectrum, bc::OBC; kwargs...)
    N = _bc_size(bc, kwargs)
    return _kitaev1d_bdg_spectrum(N, model.μ, model.t, model.Δ)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy: thermodynamic limit (Infinite, T = 0)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::Kitaev1D, ::Energy{:per_site}, ::Infinite) -> Float64

Ground-state energy per site of the infinite Kitaev1D chain at `T = 0`,

```math
\\varepsilon_0 = -\\frac{1}{2\\pi} \\int_{-\\pi}^{\\pi} \\frac{E(k)}{2}\\, dk,
\\qquad
E(k) = \\sqrt{(2t\\cos k + \\mu)^2 + 4\\Delta^2 \\sin^2 k}.
```

(The factor `1/2` accounts for the BdG particle-hole doubling: only the
negative-energy band is filled.)

Computed by adaptive Gauss-Kronrod quadrature.
"""
function fetch(model::Kitaev1D, ::Energy{:per_site}, ::Infinite; kwargs...)
    μ = model.μ
    t = model.t
    Δ = model.Δ
    result, _ = quadgk(
        k -> begin
            Ek = sqrt((2t * cos(k) + μ)^2 + 4Δ^2 * sin(k)^2)
            -Ek / 2
        end, -π, π; rtol=1e-10
    )
    return result / (2π)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Mass gap (lowest single-quasiparticle energy in the bulk)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::Kitaev1D, ::MassGap, ::Infinite) -> Float64

Bulk single-quasiparticle gap of the infinite Kitaev1D chain,

```math
\\Delta_{\\mathrm{gap}}
  = \\min_k \\sqrt{(2t\\cos k + \\mu)^2 + 4\\Delta^2 \\sin^2 k}.
```

Closed form (for `t ≠ 0`, `Δ ≠ 0`):

- `|μ| ≥ 2|t|`: minimum at `k = 0` or `k = π`, giving
  `Δ_gap = ||μ| - 2|t||`.
- `|μ| < 2|t|`: stationary point at `cos k* = -μ t / (2(t² - Δ²))`
  when `|t| ≠ |Δ|` and `|cos k*| ≤ 1`; otherwise the minimum is at
  `k = 0` or `k = π`.

`Δ = 0` gives the gapless metal whenever `|μ| < 2|t|`, and the routine
returns `0.0` in that case.

`MassGap` at `OBC` is provided as the smallest non-negative BdG
eigenvalue (numerically equal to [`EdgeModeEnergy`](@ref) at `OBC`).
"""
function fetch(model::Kitaev1D, ::MassGap, ::Infinite; kwargs...)
    μ = model.μ
    t = model.t
    Δ = model.Δ
    # Δ = 0: gapless metal whenever |μ| < 2|t|
    if Δ == 0.0
        return abs(μ) >= 2 * abs(t) ? abs(μ) - 2 * abs(t) : 0.0
    end
    # General case: minimise (2t cos k + μ)² + 4Δ² sin² k over k ∈ [0, π].
    # Stationarity in c = cos k:
    #   d/dc[(2t c + μ)² + 4Δ²(1 - c²)] = 4t(2tc + μ) - 8Δ² c = 0
    #   ⇒  c* = -μ t / (2(t² - Δ²))     (t ≠ ±Δ)
    a = t^2 - Δ^2
    if a == 0.0
        # |t| = |Δ|: dispersion is monotone in cos k; minima at k = 0 or π.
        return min(abs(2t + μ), abs(2t - μ))
    end
    c_star = -μ * t / (2 * a)
    # Always include k = 0 and k = π:
    candidates = Float64[abs(2t + μ), abs(2t - μ)]
    if -1.0 <= c_star <= 1.0
        s_sq = 1 - c_star^2
        push!(candidates, sqrt((2t * c_star + μ)^2 + 4Δ^2 * s_sq))
    end
    return minimum(candidates)
end

"""
    fetch(model::Kitaev1D, ::MassGap, bc::OBC; N::Int) -> Float64

Single-quasiparticle gap of the `N`-site OBC Kitaev1D chain — the
smallest non-negative BdG eigenvalue of the 2N × 2N BdG matrix.

In the topological phase this is the Majorana edge-mode energy
`~ e^{-N/ξ}` (use [`EdgeModeEnergy`](@ref) for the same value under a
boundary-mode-explicit name).  In the trivial phase it converges to the
bulk gap as `N → ∞`.
"""
function fetch(model::Kitaev1D, ::MassGap, bc::OBC; kwargs...)
    if model.Δ == 0.0 && abs(model.μ) < 2 * abs(model.t)
        return error(
            "Kitaev1D MassGap@OBC: Δ = 0 with |μ| < 2|t| is the gapless metal " *
            "regime — the dispersion has zeros at k_F = ±arccos(-μ/2t) and the " *
            "BdG spectrum lowest level is a finite-size remnant of those Fermi " *
            "points, not a physical gap.  Refusing to silently mask this with a " *
            "misleading number; re-evaluate with Δ ≠ 0.",
        )
    end
    N = _bc_size(bc, kwargs)
    Λ = _kitaev1d_bdg_spectrum(N, model.μ, model.t, model.Δ)
    return Λ[1]
end

# ═══════════════════════════════════════════════════════════════════════════════
# Edge-mode energy (OBC) — same value as MassGap@OBC, named for boundary modes
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::Kitaev1D, ::EdgeModeEnergy, bc::OBC; N::Int) -> Float64

Energy of the lowest-lying boundary mode on an `N`-site OBC Kitaev1D
chain — the smallest non-negative BdG eigenvalue.

In the topological phase (`|μ| < 2|t|`, `Δ ≠ 0`) the two end-localised
Majorana modes hybridise into a single complex fermion with
exponentially-small splitting `~ e^{-N/ξ}` where `ξ ~ 1/log(2|t|/|μ|)`
for `|μ| ≪ 2|t|`.

Numerically equal to `fetch(model, MassGap(), OBC(N))`; the two methods
exist as separate names so call sites can be explicit about which
physical interpretation they have in mind.
"""
function fetch(model::Kitaev1D, ::EdgeModeEnergy, bc::OBC; kwargs...)
    if model.Δ == 0.0 && abs(model.μ) < 2 * abs(model.t)
        return error(
            "Kitaev1D EdgeModeEnergy@OBC: Δ = 0 with |μ| < 2|t| is the gapless " *
            "metal regime — there is no SC pairing, no topological invariant, " *
            "and no Majorana edge modes; the BdG-spectrum lowest level is a " *
            "finite-size Fermi-point remnant.  Refusing to mask this with a " *
            "misleading number; re-evaluate with Δ ≠ 0.",
        )
    end
    N = _bc_size(bc, kwargs)
    Λ = _kitaev1d_bdg_spectrum(N, model.μ, model.t, model.Δ)
    return Λ[1]
end

# ═══════════════════════════════════════════════════════════════════════════════
# Correlation length (Infinite, T = 0)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::Kitaev1D, ::CorrelationLength, ::Infinite) -> Float64

`T = 0` correlation length of the infinite Kitaev1D chain, set by the
inverse bulk gap,

```math
\\xi = \\frac{1}{\\Delta_{\\mathrm{gap}}}.
```

Returns `Inf` on the gapless line `|μ| = 2|t|` (and on the gapless metal
`Δ = 0`, `|μ| < 2|t|`).  In QAtlas convention `ξ` is dimensionless (in
units of the lattice spacing).
"""
function fetch(model::Kitaev1D, ::CorrelationLength, ::Infinite; kwargs...)
    gap = fetch(model, MassGap(), Infinite())
    return gap <= 0.0 ? Inf : 1 / gap
end

# ═══════════════════════════════════════════════════════════════════════════════
# Topological invariant (Pfaffian sign at k = 0 and k = π)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _kitaev1d_majorana_bloch(μ, t, k::Symbol) -> Matrix{Float64}

Return the 2 × 2 real antisymmetric Majorana Bloch matrix `A(k)` of the
infinite Kitaev1D chain at the time-reversal-invariant momenta
`k ∈ {:zero, :pi}`.  At those momenta the pairing amplitude `2Δ sin k`
vanishes, so the matrix depends only on `μ` and `t`:

    A(k=0)  = [ 0  μ + 2t; -(μ + 2t)  0 ]   ⇒   Pf A(k=0)  = μ + 2t
    A(k=π)  = [ 0  μ - 2t; -(μ - 2t)  0 ]   ⇒   Pf A(k=π)  = μ - 2t

(Convention: the diagonal Majorana Bloch matrix at TR-invariant momenta
has off-diagonal entry equal to the on-site mass term `μ + 2t cos k`.)
This is the matrix on which the Kitaev (2001) Z₂ Pfaffian invariant is
evaluated.
"""
function _kitaev1d_majorana_bloch(μ::Float64, t::Float64, k::Symbol)::Matrix{Float64}
    cosk = k === :zero ? 1.0 : (k === :pi ? -1.0 : error("k must be :zero or :pi"))
    m = μ + 2t * cosk
    return [0.0 m; -m 0.0]
end

"""
    fetch(model::Kitaev1D, ::TopologicalInvariant, ::Infinite) -> Int

Pfaffian Z₂ invariant of the infinite Kitaev1D chain (Kitaev 2001),

```math
\\nu = \\operatorname{sgn}\\bigl[\\operatorname{Pf} A(k=0)
                                  \\cdot \\operatorname{Pf} A(k=\\pi)\\bigr]
     = \\operatorname{sgn}\\bigl[(\\mu + 2t)(\\mu - 2t)\\bigr]
     = \\operatorname{sgn}(\\mu^2 - 4t^2).
```

Returns `-1` in the topological phase (`|μ| < 2|t|`) and `+1` in the
trivial phase (`|μ| > 2|t|`).  Throws on the gapless line `|μ| = 2|t|`
(Pfaffian vanishes; invariant ill-defined) and on the gapless metal
`Δ = 0` with `|μ| < 2|t|`.

The two 2 × 2 Pfaffians are computed using the generic
[`pfaffian`](@ref) routine in `src/core/pfaffian.jl`, exercising the
same numerical machinery used elsewhere in QAtlas for free-fermion Wick
contractions.
"""
function fetch(model::Kitaev1D, ::TopologicalInvariant, ::Infinite; kwargs...)
    μ = model.μ
    t = model.t
    Δ = model.Δ
    # The Z₂ invariant is well-defined only when the bulk is gapped.
    if Δ == 0.0 && abs(μ) < 2 * abs(t)
        error(
            "TopologicalInvariant: bulk is gapless (Δ = 0 with |μ| < 2|t|); " *
            "the Pfaffian invariant is ill-defined.",
        )
    end
    pf0 = pfaffian(_kitaev1d_majorana_bloch(μ, t, :zero))
    pfπ = pfaffian(_kitaev1d_majorana_bloch(μ, t, :pi))
    prod = pf0 * pfπ
    if prod == 0.0
        error(
            "TopologicalInvariant: Pfaffian vanishes at |μ| = 2|t| = $(2 * abs(t)); " *
            "the bulk is gapless and the invariant is ill-defined " *
            "(μ = $μ, t = $t).",
        )
    end
    return prod > 0 ? 1 : -1
end
