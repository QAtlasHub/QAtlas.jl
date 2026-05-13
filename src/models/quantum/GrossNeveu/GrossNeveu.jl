# ─────────────────────────────────────────────────────────────────────────────
# GrossNeveu — 1+1-D asymptotically-free 4-fermion model (Gross-Neveu 1974).
#
# Lagrangian (O(2N) flavour symmetric):
#
#     ℒ = ψ̄ i γ^μ ∂_μ ψ + (g² / 2) (ψ̄ ψ)²,
#
# with `ψ` an O(2N)-multiplet of Majorana fermions (or, equivalently,
# an N-flavour Dirac doublet).  At `g = 0` the theory is free with
# central charge `c = N` (N free Dirac fermions, each contributing
# c = 1 in 1+1-D).  At any `g ≠ 0` dynamical chiral symmetry breaking
# generates an `O(2N)`-invariant non-perturbative mass
#
#     m_dyn ∝ Λ_UV · exp(-π / (N g²))      (large-N gap equation;
#                                           Gross-Neveu 1974),
#
# manifestly non-perturbative (the prefactor sets a renormalisation
# scheme) — the celebrated 2-D analog of QCD asymptotic freedom.
# Andrei-Lowenstein (1979) solved the model exactly by the Bethe
# ansatz.
#
# Phase 1 registered the **UV free-fermion central charge** `c = N`
# at `g = 0`.  Phase 2 (#247) adds the large-N dynamical mass
# `m_F = Λ exp(-π/(N g²))` (Gross-Neveu 1974) via `MassGap`.
# The Andrei-Lowenstein exact S-matrix and full RG-flow handling for
# `CentralCharge` at `g ≠ 0` remain tracked as Phase 3.
#
# References:
#   - D. J. Gross, A. Neveu, Phys. Rev. D 10, 3235 (1974).
#   - N. Andrei, J. H. Lowenstein, Phys. Rev. Lett. 43, 1698 (1979).
# ─────────────────────────────────────────────────────────────────────────────

"""
    GrossNeveu(; N::Integer = 1, g::Real = 0.0) <: AbstractQAtlasModel

1+1-D Gross-Neveu model (Gross-Neveu 1974) with `N` Dirac flavours
(equivalently `O(2N)` Majorana symmetry) and four-fermion coupling
`g ∈ ℝ`.  Asymptotically free for `N ≥ 1`.

Phase 1 registered the **UV free-fermion central charge** at
`g = 0`.  Phase 2 (#247) adds the large-N dynamical mass
`m_F = Λ exp(-π/(N g²))` via `MassGap`.  Andrei-Lowenstein exact
S-matrix and the chiral condensate remain tracked as Phase 3.

Quantities registered:

| Quantity                       | BC         | Method                                |
| ------------------------------ | ---------- | ------------------------------------- |
| [`CentralCharge`](@ref)        | `Infinite` | analytic (`c = N` at the free point)  |
| [`MassGap`](@ref)              | `Infinite` | analytic (`Λ exp(-π/(N g²))`, large-N)|

# References

- D. J. Gross, A. Neveu, *Phys. Rev. D* **10**, 3235 (1974).
- N. Andrei, J. H. Lowenstein, *Phys. Rev. Lett.* **43**, 1698 (1979).
"""
struct GrossNeveu <: AbstractQAtlasModel
    N::Int
    g::Float64
    function GrossNeveu(N::Integer, g::Real)
        N ≥ 1 ||
            throw(DomainError(N, "GrossNeveu requires N ≥ 1 Dirac flavours; got N = $N."))
        return new(Int(N), Float64(g))
    end
end
GrossNeveu(; N::Integer=1, g::Real=0.0) = GrossNeveu(N, g)

# ═══════════════════════════════════════════════════════════════════════════════
# UV free-fermion central charge
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::GrossNeveu, ::CentralCharge, ::Infinite; N=m.N, g=m.g) -> Int

Free-fermion UV central charge of the Gross-Neveu model:

    c(g = 0) = N

(N Dirac fermions in 1+1-D, each c = 1).  At any `g ≠ 0` the
theory is **massive in the IR** (dynamical chiral symmetry breaking;
Gross-Neveu 1974) and the IR central charge is zero; this Phase 1
exposes only the UV `g = 0` value and raises `DomainError` for
non-zero coupling, deferring the RG-flow case to Phase 2.

# References

- D. J. Gross, A. Neveu, *Phys. Rev. D* **10**, 3235 (1974).
"""
function fetch(
    m::GrossNeveu, ::CentralCharge, ::Infinite; N::Integer=m.N, g::Real=m.g, kwargs...
)
    iszero(g) || throw(
        DomainError(
            g,
            "GrossNeveu CentralCharge currently exposes only the UV free-fermion point g = 0 (c = N); the asymptotically-free RG flow to a massive IR is tracked as Phase 2.  Got g = $g.",
        ),
    )
    return N
end

# ═══════════════════════════════════════════════════════════════════════════════
# Dynamic fermion mass (large-N, Phase 2)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::GrossNeveu, ::MassGap, ::Infinite; Λ::Real, N=m.N, g=m.g) -> Float64

Dynamically generated fermion mass `m_F` in the large-N Gross-Neveu
model (Gross-Neveu 1974):

    m_F = Λ · exp(−π / (N · g²)).

`Λ` is the renormalisation scheme's UV cutoff / dimensional-transmutation
scale and is a **required keyword argument** — it is intentionally not
stored on the `GrossNeveu` struct because the choice of scheme is up to
the caller.  The exponential suppression at weak coupling is the
textbook signature of asymptotic-freedom-driven mass-gap formation.

`Λ > 0`, `g > 0`, `N ≥ 1` required (`DomainError` otherwise).

# Examples

- `fetch(GrossNeveu(; N=1, g=1.0), MassGap(), Infinite(); Λ=1.0)` → `exp(-π) ≈ 0.0432139`
- `fetch(GrossNeveu(; N=2, g=1.0), MassGap(), Infinite(); Λ=1.0)` → `exp(-π/2) ≈ 0.2078796`

# References

- D. J. Gross, A. Neveu, *Phys. Rev. D* **10**, 3235 (1974).
- N. Andrei, J. H. Lowenstein, *Phys. Rev. Lett.* **43**, 1698 (1979).
"""
function fetch(
    m::GrossNeveu, ::MassGap, ::Infinite; Λ::Real, N::Integer=m.N, g::Real=m.g, kwargs...
)
    Λ > 0 || throw(DomainError(Λ, "GrossNeveu MassGap requires Λ > 0; got Λ = $Λ."))
    g > 0 || throw(DomainError(g, "GrossNeveu MassGap requires g > 0; got g = $g."))
    N ≥ 1 || throw(DomainError(N, "GrossNeveu MassGap requires N ≥ 1; got N = $N."))
    return Λ * exp(-π / (N * g^2))
end
