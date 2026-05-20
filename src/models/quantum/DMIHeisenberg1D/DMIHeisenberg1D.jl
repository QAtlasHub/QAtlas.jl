# ─────────────────────────────────────────────────────────────────────────────
# DMIHeisenberg1D — spin-½ Heisenberg chain with Dzyaloshinskii-Moriya (issue #298)
#
# Hamiltonian:
#
#   H = J Σ_i Sᵢ · Sᵢ₊₁ + D ẑ · Σ_i (Sᵢ × Sᵢ₊₁),    J > 0,  D ∈ ℝ.
#
# Expanding the DM cross-product term along ẑ:
#
#   D ẑ · (Sᵢ × Sᵢ₊₁) = D (Sᵢˣ Sᵢ₊₁ʸ − Sᵢʸ Sᵢ₊₁ˣ)
#
# breaks SU(2) → U(1) (rotations about ẑ) and introduces a chiral
# preference for spiral order.  Two regimes:
#
#   * D = 0 — pure Heisenberg chain.  Bethe (1931) / Hulthén (1938):
#       E₀ / N = J · (1/4 − ln 2) ≈ −0.4431 J
#     Gapless, SU(2)-symmetric, des Cloizeaux-Pearson spinon continuum.
#
#   * D ≠ 0 — a site-dependent gauge rotation about ẑ,
#       Sᵢ⁺ → Sᵢ⁺ e^{i α i},   tan α = D/J
#     maps the model onto a twisted XXZ chain with
#       J_xy = √(J² + D²),  J_z = J,  twist α per bond.
#     The ground state is a spiral (Affleck-Oshikawa 1999); closed-form
#     energy density is available via Bethe ansatz on the gauged XXZ
#     but is Phase-1-too-complex.  Deferred to Phase 2.
#
# Phase 1 implementation strategy: delegate the D = 0 closed-form point to
# the existing `Heisenberg1D` model, using its legacy
# `GroundStateEnergyDensity` fetch (which is what `Heisenberg1D` exposes
# at `Infinite`).  The public surface of `DMIHeisenberg1D` itself uses
# the modern `Energy{:per_site}` axis-explicit convention.
#
# References:
#   - H. Bethe, Z. Physik 71, 205 (1931).
#   - L. Hulthén, Ark. Mat. Astron. Fys. 26A, No. 11 (1938).
#   - I. E. Dzyaloshinskii, J. Phys. Chem. Solids 4, 241 (1958).
#   - T. Moriya, Phys. Rev. 120, 91 (1960).
#   - I. Affleck, M. Oshikawa, Phys. Rev. B 60, 1038 (1999) — spiral
#     ground state and gauge transformation to twisted XXZ for the
#     spin-½ Heisenberg-DM chain.
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Spin S (this file)
#   Observable:  Spin S         (QAtlas-wide spin convention; see docs/src/conventions.md)

"""
    DMIHeisenberg1D(; J::Real = 1.0, D::Real = 0.0) <: AbstractQAtlasModel

Spin-½ Heisenberg chain with Dzyaloshinskii-Moriya interaction along ẑ:

    H = J Σ_i Sᵢ · Sᵢ₊₁ + D Σ_i (Sᵢˣ Sᵢ₊₁ʸ − Sᵢʸ Sᵢ₊₁ˣ),    J > 0,  D ∈ ℝ.

Phase 1 implements only the closed-form `D = 0` point (pure Heisenberg),
delegating to `Heisenberg1D`.  Generic `D ≠ 0` introduces a spiral
ground state via a gauge mapping to a twisted XXZ chain (Affleck-
Oshikawa 1999); closed-form energy is technically Bethe-ansatz-tractable
on the gauged model but is deferred to Phase 2, raising `DomainError`
from `fetch(..., Energy{:per_site}(), Infinite())`.

The default constructor `DMIHeisenberg1D()` lands on `J = 1`, `D = 0` —
the Bethe-Hulthén point, the only Phase-1 closed-form case.

# Fields

- `J::Float64` — Heisenberg exchange (J > 0 antiferromagnetic).
- `D::Float64` — DM coupling along ẑ (any real).

# References

- H. Bethe, Z. Physik **71**, 205 (1931).
- L. Hulthén, Ark. Mat. Astron. Fys. **26A**, No. 11 (1938).
- I. E. Dzyaloshinskii, J. Phys. Chem. Solids **4**, 241 (1958).
- T. Moriya, Phys. Rev. **120**, 91 (1960).
- I. Affleck, M. Oshikawa, Phys. Rev. B **60**, 1038 (1999).
"""
struct DMIHeisenberg1D <: AbstractQAtlasModel
    J::Float64
    D::Float64
    function DMIHeisenberg1D(J::Real, D::Real)
        J > 0 || throw(DomainError(J, "DMIHeisenberg1D requires J > 0; got J = $J."))
        return new(Float64(J), Float64(D))
    end
end
DMIHeisenberg1D(; J::Real=1.0, D::Real=0.0) = DMIHeisenberg1D(J, D)

# ═══════════════════════════════════════════════════════════════════════════════
# Energy granularity trait
# ═══════════════════════════════════════════════════════════════════════════════

native_energy_granularity(::DMIHeisenberg1D, ::Infinite) = :per_site

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: ground-state energy density — Phase-1 D = 0 delegation
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::DMIHeisenberg1D, ::Energy{:per_site}, ::Infinite;
          J=m.J, D=m.D) -> Float64

Ground-state energy density of the spin-½ Heisenberg-DM chain in the
thermodynamic limit.  Phase 1 supports only the closed-form `D = 0`
point:

- `D = 0` → Bethe-Hulthén: `E/N = J · (1/4 − ln 2)`.
  Delegates to `fetch(Heisenberg1D(), GroundStateEnergyDensity(),
  Infinite(); J=J)`.

- `D ≠ 0` → `DomainError`: spiral order via gauge rotation to twisted
  XXZ (Affleck-Oshikawa 1999); closed-form Bethe-ansatz energy
  technically available but deferred to Phase 2.

Floating-point tolerance for the `D = 0` match is `atol = 1e-12`.

Note on delegation: `Heisenberg1D` currently exposes its thermodynamic-
limit energy density via the legacy `GroundStateEnergyDensity` quantity,
not the modern `Energy{:per_site}` axis.  This wrapper bridges to the
modern axis on the public surface while keeping the closed-form constant
single-source-of-truth in the `Heisenberg1D` delegate.

# References

- L. Hulthén, *Ark. Mat. Astron. Fys.* **26A**, No. 11 (1938) —
  Bethe-Hulthén ground-state energy density.
- I. E. Dzyaloshinskii, *J. Phys. Chem. Solids* **4**, 241 (1958).
- T. Moriya, *Phys. Rev.* **120**, 91 (1960).
- I. Affleck, M. Oshikawa, *Phys. Rev. B* **60**, 1038 (1999) —
  spiral / twisted-XXZ mapping (deferred to Phase 2).
"""
function fetch(
    m::DMIHeisenberg1D, ::Energy{:per_site}, ::Infinite; J::Real=m.J, D::Real=m.D, kwargs...
)
    J > 0 || throw(DomainError(J, "DMIHeisenberg1D Energy requires J > 0; got J = $J."))
    if !iszero(D)
        throw(
            DomainError(
                D,
                "DMIHeisenberg1D Energy: closed-form ground-state energy known " *
                "only at D = 0 (pure Heisenberg, Bethe-Hulthén 1938); D ≠ 0 " *
                "introduces spiral order via gauge rotation to a twisted XXZ " *
                "chain (Affleck-Oshikawa 1999) and is deferred to Phase 2. " *
                "Got D = $D.",
            ),
        )
    end
    # D = 0 → pure Heisenberg.  Heisenberg1D exposes the Bethe-Hulthén
    # density via the legacy GroundStateEnergyDensity tag, with J as a
    # kwarg (the struct itself carries no parameters).
    return fetch(Heisenberg1D(), GroundStateEnergyDensity(), Infinite(); J=J)
end
