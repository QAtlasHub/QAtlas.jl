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
#     Phys. Rev. Lett. 42, 1698 (1979).
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
    pos = filter(e -> e >= -1e-12, vals)
    sort!(pos)
    length(pos) >= N || error(
        "_ssh_obc_spectrum: chiral symmetry gave $(length(pos)) non-negative eigenvalues, " *
        "expected ≥ $N (numerical degeneracy near a zero mode? v = $v, w = $w, N = $N).",
    )
    return pos[1:N]
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
gap equals the smallest non-negative OBC eigenvalue ([`MassGap`](@ref) /
[`EdgeModeEnergy`](@ref) at `OBC`); the full particle-hole *band* gap
`E_+ − E_−` is twice this.  Vanishes on the gapless line `|v| = |w|`.
"""
fetch(model::SSH, ::MassGap, ::Infinite; kwargs...) = abs(abs(model.v) - abs(model.w))

"""
    fetch(model::SSH, ::MassGap, bc::OBC; N::Int) -> Float64

Single-particle gap of the `N`-cell OBC SSH chain — the smallest non-negative
single-particle eigenvalue.  In the topological phase this is the edge-mode
splitting `~ e^{−N/ξ}` (numerically equal to [`EdgeModeEnergy`](@ref) at `OBC`);
in the trivial phase it converges to the single-particle gap `|v − w|` as
`N → ∞`.
"""
function fetch(model::SSH, ::MassGap, bc::OBC; kwargs...)
    N = _bc_size(bc, kwargs)
    return _ssh_obc_spectrum(N, model.v, model.w)[1]
end

# ═══════════════════════════════════════════════════════════════════════════════
# Edge-mode energy (OBC) — same value as MassGap@OBC, named for boundary modes
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::SSH, ::EdgeModeEnergy, bc::OBC; N::Int) -> Float64

Energy of the lowest-lying boundary mode on an `N`-cell OBC SSH chain — the
smallest non-negative single-particle eigenvalue.

In the topological phase (`|w| > |v|`) the two end-localised zero modes
hybridise with exponentially-small splitting `~ e^{−N/ξ}`, exactly `0` at the
`v = 0` sweet spot (the two end sites decouple).  Numerically equal to
`fetch(model, MassGap(), OBC(N))`; the two names exist so call sites can be
explicit about the boundary-mode interpretation.
"""
function fetch(model::SSH, ::EdgeModeEnergy, bc::OBC; kwargs...)
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

function _ssh_thermo_infinite(quantity::Symbol, v::Real, w::Real, beta::Real)
    if quantity === :free_energy
        integrand = k -> begin
            lambda_val = _ssh_dispersion(k, v, w)
            y = beta * lambda_val
            y + 2 * _ssh_log1pexp(-y)
        end
        val, _ = quadgk(integrand, 0.0, pi; rtol=1e-10)
        return -val / (2 * pi * beta)
    elseif quantity === :entropy
        integrand = k -> begin
            lambda_val = _ssh_dispersion(k, v, w)
            y = beta * lambda_val
            _ssh_log1pexp(-y) + y * _ssh_nF(y)
        end
        val, _ = quadgk(integrand, 0.0, pi; rtol=1e-10)
        return val / pi
    elseif quantity === :specific_heat
        integrand = k -> begin
            lambda_val = _ssh_dispersion(k, v, w)
            y = beta * lambda_val
            (y / 2)^2 * sech(y / 2)^2
        end
        val, _ = quadgk(integrand, 0.0, pi; rtol=1e-10)
        return val / pi
    else
        error("Unknown SSH thermal quantity: $quantity")
    end
end

"""
    fetch(m::SSH, ::FreeEnergy, ::Infinite; beta::Real, v=m.v, w=m.w, kwargs...) -> Float64

Per-site grand-potential density of the infinite SSH chain at inverse temperature `beta`.
"""
function fetch(
    m::SSH, ::FreeEnergy, ::Infinite; beta::Real, v::Real=m.v, w::Real=m.w, kwargs...
)
    beta > 0 ||
        throw(DomainError(beta, "SSH FreeEnergy requires beta > 0; got beta = $beta."))
    return _ssh_thermo_infinite(:free_energy, v, w, beta)
end

"""
    fetch(m::SSH, ::ThermalEntropy, ::Infinite; beta::Real, v=m.v, w=m.w, kwargs...) -> Float64

Per-site thermodynamic entropy of the infinite SSH chain at inverse temperature `beta`.
"""
function fetch(
    m::SSH, ::ThermalEntropy, ::Infinite; beta::Real, v::Real=m.v, w::Real=m.w, kwargs...
)
    beta > 0 ||
        throw(DomainError(beta, "SSH ThermalEntropy requires beta > 0; got beta = $beta."))
    return _ssh_thermo_infinite(:entropy, v, w, beta)
end

"""
    fetch(m::SSH, ::SpecificHeat, ::Infinite; beta::Real, v=m.v, w=m.w, kwargs...) -> Float64

Per-site specific heat of the infinite SSH chain at inverse temperature `beta`.
"""
function fetch(
    m::SSH, ::SpecificHeat, ::Infinite; beta::Real, v::Real=m.v, w::Real=m.w, kwargs...
)
    beta > 0 ||
        throw(DomainError(beta, "SSH SpecificHeat requires beta > 0; got beta = $beta."))
    return _ssh_thermo_infinite(:specific_heat, v, w, beta)
end
