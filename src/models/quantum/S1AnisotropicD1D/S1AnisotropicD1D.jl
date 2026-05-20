# ─────────────────────────────────────────────────────────────────────────────
# S1AnisotropicD1D — Spin-1 Heisenberg chain with single-ion (uniaxial)
# anisotropy.
#
# Hamiltonian:
#
#     H = J Σᵢ  Sᵢ · Sᵢ₊₁  +  D Σᵢ (Sᵢ^z)²,    S = 1,  J > 0,  D ∈ ℝ.
#
# Three distinguished regimes of the (J > 0) phase diagram:
#
#   * D = 0  : pure spin-1 antiferromagnetic Heisenberg chain (Haldane
#              chain) — gapped symmetry-protected topological (SPT)
#              Haldane phase with bulk gap Δ ≈ 0.41048 J
#              (White-Huse 1993, DMRG).
#   * D → +∞ : easy-plane "large-D" phase — gapped, trivial product
#              ground state ⊗ᵢ |Sᵢ^z = 0⟩.  Separated from the Haldane
#              phase by a Gaussian (c = 1) transition near D ≈ 0.97 J
#              (Chen-Roncaglia 2008; Tzeng-Yang-Hsu 2017).
#   * D → −∞ : Ising-like easy-axis Néel order — gapped, doubly
#              degenerate symmetry-broken ground states |…+1,−1,+1…⟩
#              and |…−1,+1,−1…⟩.  Separated from Haldane by an Ising
#              (c = 1/2) transition near D ≈ −0.31 J.
#
# Phase 1 of this model exposes ONLY the closed numerical reference at
# D = 0, by delegation to the existing [`S1Heisenberg1D`](@ref) entry.
# Any non-zero D raises DomainError pointing to the DMRG follow-up
# (Chen-Roncaglia 2008; Tzeng-Yang-Hsu 2017).
#
# References:
#   - F. D. M. Haldane, Phys. Lett. A 93, 464 (1983); PRL 50, 1153 (1983).
#   - S. R. White, D. A. Huse, Phys. Rev. B 48, 3844 (1993) — Δ ≈ 0.41048 J.
#   - W. Chen, K. Hida, B. C. Sanctuary, Phys. Rev. B 67, 104401 (2003).
#   - Y.-C. Chen, R. Roncaglia, J. Stat. Mech. P10024 (2008) — D-driven SPT
#     ↔ trivial transition in the S=1 chain with single-ion D.
#   - Y.-D. Tzeng, H.-H. Hung, Y.-C. Chen, M.-F. Yang, Phys. Rev. B 96,
#     205104 (2017) — phase boundary refinement & string order parameter.
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Spin S (this file)
#   Observable:  Spin S         (QAtlas-wide spin convention; see docs/src/conventions.md)

"""
    S1AnisotropicD1D(; J::Real = 1.0, D::Real = 0.0) <: AbstractQAtlasModel

Spin-1 Heisenberg chain with uniaxial single-ion anisotropy,

    H = J Σᵢ Sᵢ · Sᵢ₊₁ + D Σᵢ (Sᵢ^z)²,

with `J > 0` (antiferromagnetic) and `D ∈ ℝ` arbitrary single-ion
anisotropy strength.  Same 3-dimensional local Hilbert space as
[`S1Heisenberg1D`](@ref) but with the on-site easy-axis (`D < 0`) or
easy-plane (`D > 0`) term added.

Quantities registered (Phase 1):

| Quantity              | BC         | Method                                    |
| --------------------- | ---------- | ----------------------------------------- |
| [`MassGap`](@ref)     | `Infinite` | delegated to S1Heisenberg1D at D = 0      |

Phase 1 exposes only the `D = 0` Haldane-chain reference point,
delegated to the existing [`S1Heisenberg1D`](@ref).  Any `D ≠ 0`
raises `DomainError` — finite-D crosses Gaussian (large-D) and Ising
(Néel) transitions whose closed forms are unknown and are deferred to
the DMRG-based Phase 2 (Chen-Roncaglia 2008; Tzeng-Yang-Hsu 2017).

# References

- F. D. M. Haldane, *Phys. Lett. A* **93**, 464 (1983).
- S. R. White, D. A. Huse, *Phys. Rev. B* **48**, 3844 (1993).
- Y.-C. Chen, R. Roncaglia, *J. Stat. Mech.* P10024 (2008).
- Y.-D. Tzeng, H.-H. Hung, Y.-C. Chen, M.-F. Yang, *Phys. Rev. B* **96**, 205104 (2017).
"""
struct S1AnisotropicD1D <: AbstractQAtlasModel
    J::Float64
    D::Float64
    function S1AnisotropicD1D(J::Real, D::Real)
        J > 0 || throw(DomainError(J, "S1AnisotropicD1D requires J > 0; got J = $J."))
        return new(Float64(J), Float64(D))
    end
end
S1AnisotropicD1D(; J::Real=1.0, D::Real=0.0) = S1AnisotropicD1D(J, D)

# ═══════════════════════════════════════════════════════════════════════════════
# Mass gap — D = 0 delegation to S1Heisenberg1D
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::S1AnisotropicD1D, ::MassGap, ::Infinite;
          J = m.J, D = m.D, kwargs...) -> Float64

Bulk gap of the spin-1 anisotropic chain at the **D = 0 reference
point**, delegated to [`S1Heisenberg1D`](@ref): Haldane gap
`Δ ≈ 0.41048 J` (White-Huse 1993 DMRG, numerical-exact).

`D ≠ 0` raises `DomainError` — Phase 2.  At finite `D` the system
crosses

  * a Gaussian (c = 1) transition near `D ≈ 0.97 J` separating Haldane
    from the large-D phase (Chen-Roncaglia 2008), and
  * an Ising (c = 1/2) transition near `D ≈ −0.31 J` separating Haldane
    from the easy-axis Néel phase (Tzeng-Yang-Hsu 2017);

neither phase boundary has a closed-form expression in `D/J`.

# References

- S. R. White, D. A. Huse, *Phys. Rev. B* **48**, 3844 (1993).
- Y.-C. Chen, R. Roncaglia, *J. Stat. Mech.* P10024 (2008).
- Y.-D. Tzeng, H.-H. Hung, Y.-C. Chen, M.-F. Yang, *Phys. Rev. B* **96**, 205104 (2017).
"""
function fetch(
    m::S1AnisotropicD1D, ::MassGap, ::Infinite; J::Real=m.J, D::Real=m.D, kwargs...
)
    J > 0 || throw(DomainError(J, "S1AnisotropicD1D MassGap requires J > 0; got J = $J."))
    iszero(D) || throw(
        DomainError(
            D,
            "S1AnisotropicD1D MassGap: closed-form Haldane-gap reference available " *
            "only at D = 0 (pure spin-1 Heisenberg, White-Huse 1993 DMRG " *
            "Δ ≈ 0.41048 J). Generic D introduces easy-axis (D < 0) Néel or " *
            "large-D (D > 0) trivial phases separated from Haldane by Ising/Gaussian " *
            "transitions; required DMRG references (Chen-Roncaglia 2008, " *
            "Tzeng-Yang-Hsu 2017) are deferred to Phase 2. Got D = $D.",
        ),
    )
    return QAtlas.fetch(QAtlas.S1Heisenberg1D(; J=J), MassGap(), Infinite())
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy per site — D = 0 delegation to S1Heisenberg1D
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::S1AnisotropicD1D, ::Energy{:per_site}, ::Infinite;
          J = m.J, D = m.D, kwargs...) -> Float64

Ground-state energy per site of the spin-1 anisotropic chain at the
**D = 0 reference point**, delegated to [`S1Heisenberg1D`](@ref):
e₀ ≈ -1.40148403897 J (White-Huse 1993 DMRG, numerical-exact).

`D ≠ 0` raises `DomainError` — same Phase 2 gate as `MassGap`.

# References

- S. R. White, D. A. Huse, *Phys. Rev. B* **48**, 3844 (1993).
- Y.-C. Chen, R. Roncaglia, *J. Stat. Mech.* P10024 (2008).
- Y.-D. Tzeng, H.-H. Hung, Y.-C. Chen, M.-F. Yang, *Phys. Rev. B* **96**, 205104 (2017).
"""
function fetch(
    m::S1AnisotropicD1D,
    ::Energy{:per_site},
    ::Infinite;
    J::Real=m.J,
    D::Real=m.D,
    kwargs...,
)
    J > 0 || throw(
        DomainError(J, "S1AnisotropicD1D Energy{:per_site} requires J > 0; got J = $J.")
    )
    iszero(D) || throw(
        DomainError(
            D,
            "S1AnisotropicD1D Energy{:per_site}: closed-form reference available " *
            "only at D = 0 (pure spin-1 Heisenberg, White-Huse 1993 DMRG " *
            "e₀ ≈ -1.40148 J). Generic D introduces easy-axis (D < 0) Néel or " *
            "large-D (D > 0) trivial phases separated from Haldane by Ising/Gaussian " *
            "transitions; required DMRG references (Chen-Roncaglia 2008, " *
            "Tzeng-Yang-Hsu 2017) are deferred to Phase 2. Got D = $D.",
        ),
    )
    return QAtlas.fetch(QAtlas.S1Heisenberg1D(; J=J), Energy{:per_site}(), Infinite())
end
