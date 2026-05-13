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
# This Phase-1 entry registers only the **UV free-fermion central
# charge** `c = N` at `g = 0`.  The Andrei-Lowenstein exact S-matrix
# and the renormalisation-scheme-dependent dynamical mass require
# dedicated `SMatrix` / `BetaFunction` quantity types and are
# tracked as Phase 2.
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

Phase 1 registers only the **UV free-fermion central charge** at
`g = 0`.  Dynamical mass generation, the Andrei-Lowenstein exact
S-matrix and the chiral condensate are tracked as Phase 2.

Quantities registered:

| Quantity                       | BC         | Method                                |
| ------------------------------ | ---------- | ------------------------------------- |
| [`CentralCharge`](@ref)        | `Infinite` | analytic (`c = N` at the free point)  |

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
