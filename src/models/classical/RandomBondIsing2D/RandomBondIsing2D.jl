# ─────────────────────────────────────────────────────────────────────────────
# RandomBondIsing2D — 2D ±J Random-Bond Ising / Edwards-Anderson 1975.
#
# Hamiltonian:
#
#     H = -Σ_{⟨i,j⟩} J_{ij} σ^z_i σ^z_j,
#
# on a 2D square lattice, with quenched binary couplings
#
#     J_{ij} = +J  with probability p,
#            = -J  with probability 1-p,
#
# (Edwards-Anderson 1975).  The (T, p) phase diagram (Nishimori 1981;
# Honecker-Picco-Pujol 2001) has three distinguished features:
#
#   • Pure ferromagnet (p = 1): 2D Ising critical at
#         K_c = (1/2) log(1 + √2) ≈ 0.4407       (Onsager 1944),
#     2D Ising universality class, central charge  c = 1/2.
#
#   • Nishimori line  p_N(K) = (1 + tanh K) / 2: exact spin-glass /
#     FM-paramagnet crossover, gauge-symmetric  (Nishimori 1981).
#
#   • Multicritical Nishimori point  (T_N, p_N) ≈ (0.953(1), 0.1093(2))
#     (Honecker-Picco-Pujol 2001 Monte-Carlo); universality unsettled
#     (logarithmic CFT?  c ≈ 0.464 numerically; Ohzeki-Kawamura 2008).
#
# Phase 1 (this entry) exposes only the clean, universality-tagged
# piece — the 2D Ising universality of the FM critical line at p = 1
# (c = 1/2), delegated to `MinimalModel(4, 3)`.  Nishimori-line and
# multicritical-point physics (logarithmic CFT subtleties) are
# deferred to Phase 2.
#
# References:
#   - S. F. Edwards, P. W. Anderson, J. Phys. F 5, 965 (1975).
#   - H. Nishimori, Prog. Theor. Phys. 66, 1169 (1981).
#   - A. Honecker, M. Picco, P. Pujol, Phys. Rev. Lett. 87, 047201 (2001).
#   - L. Onsager, Phys. Rev. 65, 117 (1944).
# ─────────────────────────────────────────────────────────────────────────────

"""
    RandomBondIsing2D(; J::Real = 1.0, p::Real = 1.0) <: AbstractQAtlasModel

2D square-lattice ±J random-bond Ising model (Edwards-Anderson 1975):

    H = -Σ_{⟨i,j⟩} J_{ij} σ^z_i σ^z_j,
    J_{ij} = +J  (prob. p)  or  -J  (prob. 1-p).

The default `p = 1` is the pure ferromagnet — Onsager's 2D Ising
critical point with central charge `c = 1/2`.

Quantities registered (Phase 1):

| Quantity              | BC         | Method                                                    |
| --------------------- | ---------- | --------------------------------------------------------- |
| [`CentralCharge`](@ref) | `Infinite` | delegation → `MinimalModel(4, 3)` (Ising) at `p = 1`     |

Nishimori-line and multicritical Nishimori-point universality
(logarithmic CFT, Honecker-Picco-Pujol 2001) are Phase 2.

# References

- S. F. Edwards, P. W. Anderson, *J. Phys. F* **5**, 965 (1975).
- H. Nishimori, *Prog. Theor. Phys.* **66**, 1169 (1981).
- A. Honecker, M. Picco, P. Pujol, *Phys. Rev. Lett.* **87**, 047201 (2001).
"""
struct RandomBondIsing2D <: AbstractQAtlasModel
    J::Float64
    p::Float64
    function RandomBondIsing2D(J::Real, p::Real)
        J > 0 || throw(DomainError(J, "RandomBondIsing2D requires J > 0; got J = $J."))
        (0 ≤ p ≤ 1) ||
            throw(DomainError(p, "RandomBondIsing2D requires p ∈ [0, 1]; got p = $p."))
        return new(Float64(J), Float64(p))
    end
end
RandomBondIsing2D(; J::Real=1.0, p::Real=1.0) = RandomBondIsing2D(J, p)

# ═══════════════════════════════════════════════════════════════════════════════
# Central charge
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::RandomBondIsing2D, ::CentralCharge, ::Infinite; p=m.p) -> Rational

Central charge of the 2D ±J random-bond Ising model along the
phase-diagram slice selected by `p`.

Phase 1 exposes only the pure ferromagnetic critical line `p = 1`,
where the model is in the 2D Ising universality class (Onsager 1944)
with `c = 1/2`, delegated to `MinimalModel(4, 3)`.

For `0 ≤ p < 1` the Nishimori-line ferromagnet/paramagnet crossover
and the multicritical Nishimori point at `p ≈ 0.1093`
(Honecker-Picco-Pujol 2001) have unsettled / logarithmic-CFT
universality and are deferred to Phase 2 — these calls raise
`DomainError`.

# References

- L. Onsager, *Phys. Rev.* **65**, 117 (1944).
- A. Honecker, M. Picco, P. Pujol, *Phys. Rev. Lett.* **87**, 047201 (2001).
"""
function fetch(m::RandomBondIsing2D, ::CentralCharge, ::Infinite; p::Real=m.p, kwargs...)
    if p == 1
        # Pure ferromagnetic critical line: 2D Ising universality.
        return QAtlas.fetch(QAtlas.MinimalModel(4, 3), CentralCharge())
    elseif p ≥ 0 && p ≤ 1
        throw(
            DomainError(
                p,
                "RandomBondIsing2D CentralCharge: Phase 1 exposes only the pure ferromagnetic " *
                "limit p = 1 (2D Ising universality, c = 1/2). The Nishimori-line " *
                "ferromagnet-paramagnet crossover and the multicritical Nishimori point " *
                "(Honecker-Picco-Pujol 2001) have unsettled / logarithmic-CFT universality " *
                "and are deferred to Phase 2. Got p = $p.",
            ),
        )
    else
        throw(DomainError(p, "RandomBondIsing2D requires 0 ≤ p ≤ 1; got p = $p."))
    end
end
