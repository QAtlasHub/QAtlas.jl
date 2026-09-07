# ─────────────────────────────────────────────────────────────────────────────
# Rice-Mele (1982) 1D dimerised chain with a staggered on-site potential —
# exact solution.
#
# Hamiltonian (spinless fermions, two sites A/B per unit cell, N unit cells):
#
#   H = Σᵢ [ v c†_{i,A} c_{i,B} + w c†_{i,B} c_{i+1,A} + h.c. ]
#     + Δ Σᵢ ( n_{i,A} − n_{i,B} )
#
# • v — intracell hopping     • w — intercell hopping     • Δ — ±Δ on A/B
#
# Bloch h(k) = [[Δ, q(k)]; [q*(k), −Δ)]] with q(k) = v + w e^{ik}, so
#
#   E_±(k) = ± √(v² + w² + 2 v w cos k + Δ²).
#
# Δ = 0 is [`SSH`](@ref). Δ ≠ 0 breaks the BULK chiral symmetry, so there is no
# `TopologicalInvariant` row: SSH's W ∈ {0,1} is protected by exactly that
# symmetry. What is quantised here is the Thouless pump — a Zak-phase change
# around a loop in (v−w, Δ), a property of a PATH and not of a point.
#
# It does NOT break the finite-chain ± symmetry. A chain of complete cells has a
# palindromic bond sequence, so inversion R negates the staggered diagonal while
# the sublattice sign Γ negates the hopping, and S = ΓR anticommutes with H.
# `test_ricemele.jl` measures this, and measures that an incomplete cell loses it.
#
# The `(hx, hy, hz)` form of the nonlinear-response literature —
# `d(k) = (−hx cos(k/2), −hy sin(k/2), −hz)` on `k ∈ [−2π, 2π)` — is this model
# with v = (hx+hy)/2, w = (hx−hy)/2, Δ = hz.
#
# References:
#   - M. J. Rice and E. J. Mele, "Elementary Excitations of a Linearly
#     Conjugated Diatomic Polymer", Phys. Rev. Lett. 49, 1455 (1982).
#   - D. J. Thouless, "Quantization of particle transport",
#     Phys. Rev. B 27, 6083 (1983).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: tight-binding hopping amplitudes and a staggered on-site energy (this file)
#   Spinless fermions — there is no spin observable (this is a charge model).
#   Δ is the on-site energy itself, not half of it: the A/B splitting is 2Δ.

using LinearAlgebra: eigvals, Symmetric
using QuadGK: quadgk

"""
    RiceMele(; v::Real = 1.0, w::Real = 1.0, Δ::Real = 0.0) <: AbstractQAtlasModel

Rice-Mele dimerised chain with a staggered on-site potential: intracell hopping
`v`, intercell hopping `w`, and `±Δ` alternating on the two sublattices.

`Δ = 0` recovers [`SSH`](@ref) exactly.  `Δ ≠ 0` breaks the chiral symmetry that
protects SSH's winding number, so there is **no** `TopologicalInvariant` row, and
the gap `√((|v|−|w|)² + Δ²)` never closes.

`E_±(k) = ±√(v² + w² + 2vw cos k + Δ²)`.  As in [`SSH`](@ref), `MassGap` is the
Fermi-level-to-band-edge gap `min_k E_+(k)`; the band gap is twice it.

Currently registered fetches:

| Quantity                    | BC                 | Coverage                                                       |
| --------------------------- | ------------------ | -------------------------------------------------------------- |
| [`ExactSpectrum`](@ref)     | `OBC`              | All `2N` single-particle energies by dense diagonalization     |
| [`Energy`](@ref)`{:per_site}` | `Infinite`       | Half-filled energy density by Gauss-Kronrod dispersion integral |
| [`Energy`](@ref)`{:per_site}` | `OBC`            | Half-filled energy density, `N` lowest levels over `2N` sites   |
| [`MassGap`](@ref)           | `Infinite` / `OBC` | Fermi-level-to-band-edge gap                                    |
| [`CorrelationLength`](@ref) | `Infinite`         | Inverse of the single-particle gap                              |
"""
struct RiceMele <: AbstractQAtlasModel
    v::Float64
    w::Float64
    Δ::Float64
end

function RiceMele(; v::Real=1.0, w::Real=1.0, Δ::Real=0.0)
    (isfinite(v) && isfinite(w) && isfinite(Δ)) || throw(
        ArgumentError("RiceMele: v, w and Δ must be finite; got v = $v, w = $w, Δ = $Δ")
    )
    return RiceMele(Float64(v), Float64(w), Float64(Δ))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy granularity convention (see src/core/quantities.jl)
# ═══════════════════════════════════════════════════════════════════════════════

native_energy_granularity(::RiceMele, ::Infinite) = :per_site
native_energy_granularity(::RiceMele, ::OBC) = :per_site

# ═══════════════════════════════════════════════════════════════════════════════
# Internal: dispersion and the finite-chain single-particle spectrum
# ═══════════════════════════════════════════════════════════════════════════════

"`E_+(k) = √(v² + w² + 2vw cos k + Δ²)`, the upper band."
function _rice_mele_dispersion(k, v::Float64, w::Float64, Δ::Float64)
    return sqrt(v^2 + w^2 + 2 * v * w * cos(k) + Δ^2)
end

"""
    _rice_mele_obc_spectrum(N, v, w, Δ) -> Vector{Float64}

All `2N` single-particle energies of the OBC Rice-Mele chain with `N` unit cells,
sorted ascending.

The `2N × 2N` Hamiltonian is real symmetric tridiagonal: off-diagonals
`v, w, v, w, …, v` (odd bonds `v`, even bonds `w`) and diagonal `+Δ, −Δ, +Δ, …`.

All `2N`, sorted ascending.
"""
function _rice_mele_obc_spectrum(N::Int, v::Float64, w::Float64, Δ::Float64)
    N >= 1 || throw(ArgumentError("RiceMele: need N ≥ 1 unit cells; got N = $N"))
    return _rice_mele_obc_spectrum_sites(2N, v, w, Δ)
end

"""
    _rice_mele_obc_spectrum_sites(n, v, w, Δ) -> Vector{Float64}

The same chain by SITE count, so an incomplete cell can be built. Not reachable
through `fetch`, which speaks only in cells; it exists so the `±`-symmetry claim
has a test that can see it fail.
"""
function _rice_mele_obc_spectrum_sites(n::Int, v::Float64, w::Float64, Δ::Float64)
    n >= 1 || throw(ArgumentError("RiceMele: need n ≥ 1 sites; got n = $n"))
    H = zeros(n, n)
    @inbounds for j in 1:n
        H[j, j] = isodd(j) ? Δ : -Δ          # site j: odd ⇒ sublattice A (+Δ)
    end
    @inbounds for j in 1:(n - 1)
        t = isodd(j) ? v : w                 # bond (j, j+1): odd ⇒ intracell v
        H[j, j + 1] = t
        H[j + 1, j] = t
    end
    return eigvals(Symmetric(H))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Exact spectrum (OBC)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::RiceMele, ::ExactSpectrum, bc::OBC; N::Int) -> Vector{Float64}

All `2N` single-particle energies of the OBC Rice-Mele chain (`N` unit cells,
`2N` sites), sorted ascending.

The whole spectrum, unlike [`SSH`](@ref)'s non-negative half — returning half
would hide its loss from a caller who asks for an incomplete cell.
"""
function fetch(model::RiceMele, ::ExactSpectrum, bc::OBC; kwargs...)
    N = _bc_size(bc, kwargs)
    return _rice_mele_obc_spectrum(N, model.v, model.w, model.Δ)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy (half filling)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::RiceMele, ::Energy{:per_site}, ::Infinite) -> Float64

Ground-state energy per site of the infinite chain at `T = 0` and half filling
(lower band full),

```math
\\varepsilon_0 = -\\frac{1}{4\\pi} \\int_{-\\pi}^{\\pi} E_+(k)\\, dk,
\\qquad E_+(k) = \\sqrt{v^2 + w^2 + 2vw\\cos k + \\Delta^2}.
```

The `1/4π` rather than `1/2π` divides the per-unit-cell band energy by the two
sites in a cell.  Computed by adaptive Gauss-Kronrod quadrature.
"""
function fetch(model::RiceMele, ::Energy{:per_site}, ::Infinite; kwargs...)
    v, w, Δ = model.v, model.w, model.Δ
    result, _ = quadgk(k -> _rice_mele_dispersion(k, v, w, Δ), -π, π; rtol=1e-10)
    return -result / (4π)
end

"""
    fetch(model::RiceMele, ::Energy{:per_site}, bc::OBC; N::Int) -> Float64

Ground-state energy per site of the `N`-cell OBC chain at half filling: the `N`
lowest single-particle levels over `2N` sites.

Carries the open-boundary correction, so it approaches `Infinite` as `O(1/N)`
rather than equalling it — which is what a finite-chain variational calculation
should be compared against.
"""
function fetch(model::RiceMele, ::Energy{:per_site}, bc::OBC; kwargs...)
    N = _bc_size(bc, kwargs)
    ε = _rice_mele_obc_spectrum(N, model.v, model.w, model.Δ)
    return sum(view(ε, 1:N)) / (2N)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Mass gap
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::RiceMele, ::MassGap, ::Infinite) -> Float64

Fermi-level-to-band-edge gap of the infinite chain,

```math
\\Delta_{\\mathrm{gap}} = \\min_k E_+(k) = \\sqrt{(|v| - |w|)^2 + \\Delta^2},
```

the minimum sitting at `k = π` for `vw > 0` and `k = 0` for `vw < 0`.  The
particle-hole *band* gap is twice this, matching [`SSH`](@ref)'s convention, to
which it reduces at `Δ = 0`.

Never zero for `Δ ≠ 0`: the staggered potential gaps out SSH's `|v| = |w|` Dirac
point.
"""
function fetch(model::RiceMele, ::MassGap, ::Infinite; kwargs...)
    return sqrt((abs(model.v) - abs(model.w))^2 + model.Δ^2)
end

"""
    fetch(model::RiceMele, ::MassGap, bc::OBC; N::Int) -> Float64

Fermi-level-to-band-edge gap of the `N`-cell OBC chain: half the HOMO-LUMO
difference `(ε_{N+1} − ε_N)/2` at half filling — the `Infinite` convention.

Written by locating the Fermi level rather than assuming it at zero, which is
[`SSH`](@ref)'s phrasing. The two agree at every `Δ`.

Does NOT converge to `MassGap` at `Infinite` when `|w| > |v|`: there it is the
edge-mode scale, flat in `N` from about `N = 10`.  Same caveat as [`SSH`](@ref).
"""
function fetch(model::RiceMele, ::MassGap, bc::OBC; kwargs...)
    N = _bc_size(bc, kwargs)
    ε = _rice_mele_obc_spectrum(N, model.v, model.w, model.Δ)
    return (ε[N + 1] - ε[N]) / 2
end

# ═══════════════════════════════════════════════════════════════════════════════
# Correlation length
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::RiceMele, ::CorrelationLength, ::Infinite) -> Float64

`T = 0` correlation length of the infinite chain, the inverse single-particle
gap `ξ = 1/√((|v|−|w|)² + Δ²)`.

Finite for every `Δ ≠ 0`, unlike [`SSH`](@ref) on the `|v| = |w|` line.
"""
function fetch(model::RiceMele, ::CorrelationLength, ::Infinite; kwargs...)
    gap = fetch(model, MassGap(), Infinite())
    return gap <= 0.0 ? Inf : 1 / gap
end
