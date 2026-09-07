# ─────────────────────────────────────────────────────────────────────────────
# Su-Schrieffer-Heeger (1979) 1D dimerised tight-binding chain — exact solution.
#
# Hamiltonian (spinless fermions, two sites A/B per unit cell, N unit cells):
#
#   H = Σᵢ [ v c†_{i,A} c_{i,B} + w c†_{i,B} c_{i+1,A} + h.c. ]
#
# • v — intracell hopping (A ↔ B within a cell)
# • w — intercell hopping (B ↔ A of the next cell)
#
# Particle-conserving (no pairing): the 2N × 2N single-particle Hamiltonian is a
# real symmetric tridiagonal hopping matrix with alternating off-diagonals
# v, w, v, w, …, v (N copies of v, N−1 of w) and zero on-site energy.
#
# Bloch Hamiltonian h(k) = [[0, q(k)]; [q*(k), 0]] with q(k) = v + w e^{ik},
# so the two bands are
#
#   E_±(k) = ± |q(k)| = ± √(v² + w² + 2 v w cos k),
#
# gapped everywhere except k = π when v = w.  Phase diagram (chiral symmetry
# class BDI):
#
#   |w| > |v|  → topological (winding W = 1; two near-zero edge modes at OBC
#                 ends, splitting ~ e^{−N/ξ}; exactly zero at the v = 0 sweet
#                 spot for any N).
#   |w| < |v|  → trivial (winding W = 0; no edge modes).
#   |w| = |v|  → gapless Dirac point (k = π if vw>0, k = 0 if vw<0; winding ill-defined).
#
# References:
#   - W. P. Su, J. R. Schrieffer, A. J. Heeger,
#     [SSH1979](@cite).
#   - J. K. Asbóth, L. Oroszlány, A. Pályi, "A Short Course on Topological
#     Insulators", Lect. Notes Phys. 919 (2016) Ch. 1.
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: tight-binding hopping amplitudes (this file)
#   Spinless fermions — there is no spin observable (this is a charge model).

using LinearAlgebra: eigvals, Symmetric
using QuadGK: quadgk

"""
    SSH(; v = 1.0, w = 1.0) <: AbstractQAtlasModel

Su-Schrieffer-Heeger (1979) one-dimensional dimerised tight-binding chain,

```math
H = \\sum_i \\left( v\\, c_{i,A}^{\\dagger} c_{i,B}
                  + w\\, c_{i,B}^{\\dagger} c_{i+1,A} + \\text{h.c.} \\right).
```

`v` is the intracell hopping and `w` the intercell hopping.  `|w| > |v|` is the
topological phase (winding `W = 1`, edge modes); `|w| < |v|` is the trivial
phase (`W = 0`); `|w| = |v|` is the gapless Dirac point.

The two-band dispersion is `E_±(k) = ±√(v² + w² + 2 v w cos k)`; the
single-particle gap (QAtlas `MassGap`) is `min_k|q(k)| = ||v| − |w||` (`|v − w|`
for same-sign hoppings) and the band gap `E_+ − E_−` is twice that.  This is the
particle-conserving cousin of the [`Kitaev1D`](@ref) Majorana wire (both have a
chiral sublattice symmetry — SSH is class BDI for all `v,w`, Kitaev1D at its
symmetric point — and protected edge modes), without superconducting pairing.

Currently registered fetches:

| Quantity                   | BC                 | Coverage                                                              |
| -------------------------- | ------------------ | --------------------------------------------------------------------- |
| [`ExactSpectrum`](@ref)    | `OBC`              | Non-negative single-particle energies via dense diagonalization       |
| [`Energy`](@ref)           | `Infinite`         | Ground-state energy density via Gauss-Kronrod dispersion integral     |
| [`MassGap`](@ref)          | `Infinite` / `OBC` | Single-particle gap                                                   |
| [`CorrelationLength`](@ref)| `Infinite`         | Inverse of single-particle gap                                        |
| [`TopologicalInvariant`](@ref) | `Infinite`     | Winding number of the chiral off-diagonal loop                        |
| [`FreeEnergy`](@ref)        | `Infinite`         | Thermal Helmholtz free energy density                                 |
| [`ThermalEntropy`](@ref)    | `Infinite`         | Thermal entropy density                                               |
| [`SpecificHeat`](@ref)      | `Infinite`         | Thermal specific heat density                                         |
| [`UniversalityClass`](@ref) | `Infinite`         | `:XY` universality class on the critical/gapless line `|v| = |w|`     |
"""
struct SSH <: AbstractQAtlasModel
    v::Float64
    w::Float64
end
function SSH(; v::Real=1.0, w::Real=1.0)
    (isfinite(v) && isfinite(w)) ||
        throw(ArgumentError("SSH: v and w must be finite; got v = $v, w = $w"))
    return SSH(Float64(v), Float64(w))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Internal: single-particle spectrum (OBC, finite N unit cells = 2N sites)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _ssh_obc_spectrum(N, v, w) -> Vector{Float64}

Return the `N` non-negative single-particle energies of the OBC SSH chain with
`N` unit cells (`2N` sites) and hoppings `(v, w)`, sorted ascending.

The `2N × 2N` Hamiltonian is the real symmetric tridiagonal hopping matrix with
zero diagonal and off-diagonal entries `v, w, v, w, …, v` (odd bonds `= v`,
even bonds `= w`).  Chiral (sublattice) symmetry makes the spectrum symmetric
about zero, so the non-negative half (the `N` particle-branch energies) fully
determines it.  In the topological phase the smallest entry is the
exponentially-small edge-mode splitting `~ e^{−N/ξ}` (exactly `0` at `v = 0`).
"""
function _ssh_obc_spectrum(N::Int, v::Float64, w::Float64)::Vector{Float64}
    N >= 1 || throw(ArgumentError("SSH: need N ≥ 1 unit cells; got N = $N"))
    n = 2N
    H = zeros(n, n)
    @inbounds for j in 1:(n - 1)
        t = isodd(j) ? v : w          # bond (j, j+1): odd ⇒ intracell v, even ⇒ intercell w
        H[j, j + 1] = t
        H[j + 1, j] = t
    end
    vals = eigvals(Symmetric(H))
    sort!(vals)
    # The UPPER half. Not `filter(non-negative)[1:N]`, which at v = 0 keeps both members
    # of the exact zero pair and drops a particle level.
    half = vals[(N + 1):(2N)]
    # Taking a half only means anything if the spectrum is ± symmetric, which chiral
    # symmetry gives structurally. Assert it, so a future matrix that breaks it says so
    # instead of returning a plausible wrong set.
    scale = max(1.0, maximum(abs, vals))
    maximum(abs, vals[1:N] .+ reverse(half)) < 1.0e-8 * scale ||
        error("_ssh_obc_spectrum: spectrum is not ± symmetric; v = $v, w = $w, N = $N")
    return half
end

# |q(k)| = √(v² + w² + 2 v w cos k), the upper band energy at momentum k.
_ssh_dispersion(k::Real, v::Float64, w::Float64) = sqrt(v^2 + w^2 + 2 * v * w * cos(k))

# ═══════════════════════════════════════════════════════════════════════════════
# Energy granularity convention (see src/core/quantities.jl)
# ═══════════════════════════════════════════════════════════════════════════════

native_energy_granularity(::SSH, ::Infinite) = :per_site

# ═══════════════════════════════════════════════════════════════════════════════
# ExactSpectrum (OBC) — non-negative single-particle energies
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::SSH, ::ExactSpectrum, bc::OBC; N::Int) -> Vector{Float64}

The `N` non-negative single-particle energies of the OBC SSH chain (`N` unit
cells, `2N` sites), sorted ascending.  By chiral symmetry the full spectrum is
`±` this set.  In the topological phase (`|w| > |v|`) the lowest entry is the
edge-mode splitting `~ e^{−N/ξ}` (exactly `0` at `v = 0`).
"""
function fetch(model::SSH, ::ExactSpectrum, bc::OBC; kwargs...)
    N = _bc_size(bc, kwargs)
    return _ssh_obc_spectrum(N, model.v, model.w)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy: thermodynamic limit (Infinite, T = 0, half filling)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::SSH, ::Energy{:per_site}, ::Infinite) -> Float64

Ground-state energy per site of the infinite SSH chain at `T = 0` and half
filling (the lower band fully occupied),

```math
\\varepsilon_0 = -\\frac{1}{4\\pi} \\int_{-\\pi}^{\\pi} |q(k)|\\, dk,
\\qquad |q(k)| = \\sqrt{v^2 + w^2 + 2 v w \\cos k}.
```

The `1/4π` (rather than `1/2π`) divides the per-unit-cell band energy by the
two sites per cell.  At a fully dimerised sweet spot (`v = 0` or `w = 0`) the
band is flat and `ε₀ = −max(|v|, |w|)/2`.  Computed by adaptive Gauss-Kronrod
quadrature.
"""
function fetch(model::SSH, ::Energy{:per_site}, ::Infinite; kwargs...)
    v = model.v
    w = model.w
    result, _ = quadgk(k -> _ssh_dispersion(k, v, w), -π, π; rtol=1e-10)
    return -result / (4π)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Mass gap (bulk single-particle gap)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::SSH, ::MassGap, ::Infinite) -> Float64

Single-particle gap of the infinite SSH chain — the lowest positive
single-particle energy `min_k |q(k)|`:

```math
\\Delta_{\\mathrm{gap}} = \\min_k |q(k)| = \\bigl| |v| - |w| \\bigr|
```

(the minimum sits at `k = π` when `vw > 0` and at `k = 0` when `vw < 0`; for
same-sign hoppings this reduces to `|v − w|`).  This Fermi-level-to-band-edge
gap equals the smallest non-negative OBC eigenvalue ([`MassGap`](@ref) at
`OBC`); the full particle-hole *band* gap
`E_+ − E_−` is twice this.  Vanishes on the gapless line `|v| = |w|`.
"""
fetch(model::SSH, ::MassGap, ::Infinite; kwargs...) = abs(abs(model.v) - abs(model.w))

"""
    fetch(model::SSH, ::MassGap, bc::OBC; N::Int) -> Float64

Single-particle gap of the `N`-cell OBC SSH chain — the smallest non-negative
single-particle eigenvalue.

In the topological phase (`|w| > |v|`) this IS the edge-mode splitting: the two
end-localised zero modes hybridise as `~ e^{−N/ξ}`, and exactly `0` at the
`v = 0` sweet spot, where the end sites decouple.  In the trivial phase the same
number converges to the single-particle gap `|v − w|` as `N → ∞`.

The boundary-mode reading is an *interpretation of this value in a phase*, not a
second quantity — it is recorded on the registry row (#816).
"""
function fetch(model::SSH, ::MassGap, bc::OBC; kwargs...)
    N = _bc_size(bc, kwargs)
    return _ssh_obc_spectrum(N, model.v, model.w)[1]
end

# ═══════════════════════════════════════════════════════════════════════════════
# Correlation length (Infinite, T = 0)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::SSH, ::CorrelationLength, ::Infinite) -> Float64

`T = 0` correlation length of the infinite SSH chain, set by the inverse
single-particle gap, `ξ = 1 / ||v| − |w||`.  Returns `Inf` on the gapless line
`|v| = |w|`.
"""
function fetch(model::SSH, ::CorrelationLength, ::Infinite; kwargs...)
    gap = fetch(model, MassGap(), Infinite())
    return gap <= 0.0 ? Inf : 1 / gap
end

# ═══════════════════════════════════════════════════════════════════════════════
# Topological invariant (winding number of q(k) = v + w e^{ik})
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::SSH, ::TopologicalInvariant, ::Infinite) -> Int

Winding number of the chiral off-diagonal `q(k) = v + w e^{ik}` around the
origin as `k : 0 → 2π`,

```math
W = \\frac{1}{2\\pi} \\oint \\mathrm{Im}\\frac{q'(k)}{q(k)}\\, dk
  = \\begin{cases} 1 & |w| > |v| \\;\\text{(topological)} \\\\
                   0 & |w| < |v| \\;\\text{(trivial)} \\end{cases}.
```

Computed by Gauss-Kronrod quadrature of the argument-derivative `Im(q'/q)`
(`q' = i w e^{ik}`); the magnitude `|W| ∈ {0, 1}` is returned (the integral's
sign flips with `sign(w)` and is not physical here).  Throws on the gapless line
`|v| = |w|` — where `q` passes through the origin (`q(π) = 0` if `vw > 0`,
`q(0) = 0` if `vw < 0`) and the winding is ill-defined — and if the integral
fails to land near an integer (near-gapless parameters).
"""
function fetch(model::SSH, ::TopologicalInvariant, ::Infinite; kwargs...)
    v = model.v
    w = model.w
    gap = abs(abs(v) - abs(w))
    scale = max(abs(v), abs(w), 1.0)
    gap <= 1e-8 * scale && error(
        "SSH TopologicalInvariant: |v| ≈ |w| (gap $(gap) ≪ scale $(scale)) — q passes " *
        "through the origin and the winding number is ill-defined (v = $v, w = $w).",
    )
    # W = (1/2π) ∮ Im(q'/q) dk with q = v + w e^{ik}, q' = i w e^{ik}.
    integral, _ = quadgk(k -> imag((im * w * cis(k)) / (v + w * cis(k))), 0, 2π; rtol=1e-10)
    raw = integral / (2π)
    abs(raw - round(raw)) > 0.25 && error(
        "SSH TopologicalInvariant: winding integral ($(raw)) is not near an integer — " *
        "likely near-gapless parameters (v = $v, w = $w).",
    )
    return abs(round(Int, raw))
end

# =========================================================================
# Finite-T Thermodynamics
# =========================================================================

# Numerically stable log(1 + exp(y))
@inline function _ssh_log1pexp(y::Real)
    return y > 0 ? y + log1p(exp(-y)) : log1p(exp(y))
end

# Numerically stable Fermi-Dirac occupation n_F(x) = 1/(1 + e^x)
@inline function _ssh_nF(x::Real)
    return x > 0 ? exp(-x) / (1 + exp(-x)) : 1 / (1 + exp(x))
end

# ── Type-dispatched integrands ───────────────────────────────────────────────
#
# The quantity reaches these kernels as a type in the AbstractQAtlas vocabulary.
# Selecting the integrand from a `Symbol` discarded that, leaving `integrand` a
# union of anonymous closure types at the `quadgk` call site, so the adaptive
# quadrature machinery was instantiated once per branch of the union.  A named
# integrand struct dispatched on the quantity keeps the static information the
# vocabulary already carries.
#
# The integrand expressions and the `quadgk` calls are unchanged, so the values
# are unchanged.

# Each field keeps its OWN type, exactly as the closures this replaced captured
# them.  Promoting them to a common type would drag the model parameters up to
# whatever `beta` is — and when `beta` is a `ForwardDiff.Dual` (the derived-input
# suppliers differentiate through these), `v`/`w` become Duals too and no longer
# match `_ssh_dispersion(::Real, ::Float64, ::Float64)`.  Even where the
# dispersion is generic enough not to error, promoting runs the whole quadrature
# in Dual arithmetic for parameters that carry no derivative information.
struct _SSHIntegrand{Q,V<:Real,W<:Real,B<:Real}
    v::V
    w::W
    beta::B
end

function _SSHIntegrand{Q}(v::Real, w::Real, beta::Real) where {Q}
    return _SSHIntegrand{Q,typeof(v),typeof(w),typeof(beta)}(v, w, beta)
end

@inline function (g::_SSHIntegrand{FreeEnergy})(k)
    y = g.beta * _ssh_dispersion(k, g.v, g.w)
    return y + 2 * _ssh_log1pexp(-y)
end

@inline function (g::_SSHIntegrand{ThermalEntropy})(k)
    y = g.beta * _ssh_dispersion(k, g.v, g.w)
    return _ssh_log1pexp(-y) + y * _ssh_nF(y)
end

@inline function (g::_SSHIntegrand{SpecificHeat})(k)
    y = g.beta * _ssh_dispersion(k, g.v, g.w)
    return (y / 2)^2 * sech(y / 2)^2
end

_ssh_quad(g::_SSHIntegrand) = first(quadgk(g, 0.0, pi; rtol=1e-10))

"""
    _ssh_thermo_infinite(quantity, v, w, beta) -> Real

Per-site thermodynamic potential of the infinite SSH chain, dispatched on the
concrete quantity type.
"""
function _ssh_thermo_infinite end

function _ssh_thermo_infinite(::FreeEnergy, v::Real, w::Real, beta::Real)
    return -_ssh_quad(_SSHIntegrand{FreeEnergy}(v, w, beta)) / (2 * pi * beta)
end

function _ssh_thermo_infinite(::ThermalEntropy, v::Real, w::Real, beta::Real)
    return _ssh_quad(_SSHIntegrand{ThermalEntropy}(v, w, beta)) / pi
end

function _ssh_thermo_infinite(::SpecificHeat, v::Real, w::Real, beta::Real)
    return _ssh_quad(_SSHIntegrand{SpecificHeat}(v, w, beta)) / pi
end

"""
    fetch(m::SSH, ::FreeEnergy, ::Infinite; beta::Real, v=m.v, w=m.w, kwargs...) -> Float64

Per-site grand-potential density of the infinite SSH chain at inverse temperature `beta`.
"""
function fetch(
    m::SSH, q::FreeEnergy, ::Infinite; beta::Real, v::Real=m.v, w::Real=m.w, kwargs...
)
    beta > 0 ||
        throw(DomainError(beta, "SSH FreeEnergy requires beta > 0; got beta = $beta."))
    return _ssh_thermo_infinite(q, v, w, beta)
end

"""
    fetch(m::SSH, ::ThermalEntropy, ::Infinite; beta::Real, v=m.v, w=m.w, kwargs...) -> Float64

Per-site thermodynamic entropy of the infinite SSH chain at inverse temperature `beta`.
"""
function fetch(
    m::SSH, q::ThermalEntropy, ::Infinite; beta::Real, v::Real=m.v, w::Real=m.w, kwargs...
)
    beta > 0 ||
        throw(DomainError(beta, "SSH ThermalEntropy requires beta > 0; got beta = $beta."))
    return _ssh_thermo_infinite(q, v, w, beta)
end

"""
    fetch(m::SSH, ::SpecificHeat, ::Infinite; beta::Real, v=m.v, w=m.w, kwargs...) -> Float64

Per-site specific heat of the infinite SSH chain at inverse temperature `beta`.
"""
function fetch(
    m::SSH, q::SpecificHeat, ::Infinite; beta::Real, v::Real=m.v, w::Real=m.w, kwargs...
)
    beta > 0 ||
        throw(DomainError(beta, "SSH SpecificHeat requires beta > 0; got beta = $beta."))
    return _ssh_thermo_infinite(q, v, w, beta)
end
