# ─────────────────────────────────────────────────────────────────────────────
# J1J2Heisenberg1D — spin-½ J₁-J₂ Heisenberg chain (issue #297)
#
# Hamiltonian:
#
#   H = J₁ Σ_i Sᵢ · Sᵢ₊₁ + J₂ Σ_i Sᵢ · Sᵢ₊₂,    J₁ > 0,  J₂ ≥ 0.
#
# The frustration ratio j = J₂/J₁ controls the physics.  Two closed-form
# points are accessible analytically in Phase 1; all generic j > 0,
# j ≠ 1/2 require DMRG (deferred to Phase 2):
#
#   * j = 0   — pure Heisenberg, Bethe (1931) / Hulthén (1938):
#       E₀ / N = J₁ · (1/4 − ln 2) ≈ −0.4431 J₁
#     The ground state is gapless with des Cloizeaux-Pearson spinon
#     continuum.
#
#   * j = 1/2 — Majumdar-Ghosh (1969): exact dimer-product ground state,
#       E₀ / N = −3 J₁ / 8 = −0.375 J₁
#     Two-fold degenerate (even/odd dimer coverings), size-independent
#     on any even-length chain.
#
#   * 0 < j < 1/2 (Heisenberg-like, gapless, with logarithmic dimer
#     corrections) and j > 1/2 (gapped dimer phase incommensurate spiral
#     correlations for j > j_L ≈ 0.5206) — no closed form; require DMRG.
#     Deferred to Phase 2; raises a `DomainError` here.
#
# Phase 1 implementation strategy: delegate the two closed-form points to
# the existing `Heisenberg1D` (j = 0) and `MajumdarGhosh` (j = 1/2)
# models, using their `GroundStateEnergyDensity` fetch — the legacy tag
# is what those models currently expose at `Infinite`.  The public
# surface of `J1J2Heisenberg1D` itself uses the modern `Energy{:per_site}`
# axis-explicit convention.
#
# References:
#   - H. Bethe, Z. Physik 71, 205 (1931).
#   - L. Hulthén, Ark. Mat. Astron. Fys. 26A, No. 11 (1938).
#   - C. K. Majumdar, D. K. Ghosh, J. Math. Phys. 10, 1388 (1969).
#   - S. R. White, I. Affleck, Phys. Rev. B 54, 9862 (1996) — DMRG study
#     of the generic-j J₁-J₂ chain; identifies the dimer phase and the
#     Lifshitz point j_L ≈ 0.5206 separating commensurate / incommensurate
#     short-range order.
# ─────────────────────────────────────────────────────────────────────────────

"""
    J1J2Heisenberg1D(; J1::Real = 1.0, J2::Real = 0.5) <: AbstractQAtlasModel

Spin-½ J₁-J₂ Heisenberg chain:

    H = J₁ Σ_i Sᵢ · Sᵢ₊₁ + J₂ Σ_i Sᵢ · Sᵢ₊₂,    J₁ > 0,  J₂ ≥ 0.

Define the frustration ratio `j = J₂ / J₁`.  Phase 1 implements only the
two closed-form points:

- `j = 0`   → Bethe-Hulthén `E/N = J₁ (1/4 − ln 2)`,
- `j = 1/2` → Majumdar-Ghosh `E/N = −3 J₁ / 8`.

Generic `j` (0 < j < ∞, j ∉ {0, 1/2}) requires DMRG; deferred to Phase 2,
raises `DomainError` from `fetch(..., Energy{:per_site}(), Infinite())`.

The default constructor `J1J2Heisenberg1D()` lands on `J₁ = 1`, `J₂ = 1/2`
— the Majumdar-Ghosh point, the most physically interesting closed-form
case.

# Fields

- `J1::Float64` — nearest-neighbour antiferromagnetic exchange (J₁ > 0).
- `J2::Float64` — next-nearest-neighbour exchange (J₂ ≥ 0).

# References

- H. Bethe, Z. Physik **71**, 205 (1931).
- L. Hulthén, Ark. Mat. Astron. Fys. **26A**, No. 11 (1938).
- C. K. Majumdar, D. K. Ghosh, J. Math. Phys. **10**, 1388 (1969).
- S. R. White, I. Affleck, Phys. Rev. B **54**, 9862 (1996).
"""
struct J1J2Heisenberg1D <: AbstractQAtlasModel
    J1::Float64
    J2::Float64
    function J1J2Heisenberg1D(J1::Real, J2::Real)
        J1 > 0 || throw(DomainError(J1, "J1J2Heisenberg1D requires J1 > 0; got J1 = $J1."))
        J2 ≥ 0 ||
            throw(DomainError(J2, "J1J2Heisenberg1D requires J2 ≥ 0; got J2 = $J2."))
        return new(Float64(J1), Float64(J2))
    end
end
J1J2Heisenberg1D(; J1::Real=1.0, J2::Real=0.5) = J1J2Heisenberg1D(J1, J2)

# ═══════════════════════════════════════════════════════════════════════════════
# Energy granularity trait
# ═══════════════════════════════════════════════════════════════════════════════

native_energy_granularity(::J1J2Heisenberg1D, ::Infinite) = :per_site

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: ground-state energy density — Phase-1 closed-form delegation
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::J1J2Heisenberg1D, ::Energy{:per_site}, ::Infinite;
          J1=m.J1, J2=m.J2) -> Float64

Ground-state energy density of the spin-½ J₁-J₂ Heisenberg chain in the
thermodynamic limit.  Phase 1 supports only the two closed-form points:

- `j = J₂/J₁ = 0`   → Bethe-Hulthén: `E/N = J₁ · (1/4 − ln 2)`.
  Delegates to `fetch(Heisenberg1D(), GroundStateEnergyDensity(),
  Infinite(); J=J1)`.

- `j = 1/2`         → Majumdar-Ghosh dimer GS: `E/N = −3 J₁ / 8`.
  Delegates to `fetch(MajumdarGhosh(; J=J1), GroundStateEnergyDensity(),
  Infinite())`.

- otherwise         → `DomainError`: no closed form; numerical DMRG
  required, deferred to Phase 2.

Floating-point tolerance for the j = 0 and j = 1/2 matches is `atol =
1e-12`.

Note on delegation: `Heisenberg1D` and `MajumdarGhosh` currently expose
their thermodynamic-limit energy density via the legacy
`GroundStateEnergyDensity` quantity, not the modern `Energy{:per_site}`
axis.  This wrapper bridges to the modern axis on the public surface
while keeping the closed-form constants in a single source-of-truth
location (the delegate model).

# References

- L. Hulthén, *Ark. Mat. Astron. Fys.* **26A**, No. 11 (1938) —
  Bethe-Hulthén ground-state energy density.
- C. K. Majumdar, D. K. Ghosh, *J. Math. Phys.* **10**, 1388 (1969) —
  exact dimer ground state at j = 1/2.
- S. R. White, I. Affleck, *Phys. Rev. B* **54**, 9862 (1996) — DMRG
  study of generic j (deferred to Phase 2).
"""
function fetch(
    m::J1J2Heisenberg1D,
    ::Energy{:per_site},
    ::Infinite;
    J1::Real=m.J1,
    J2::Real=m.J2,
    kwargs...,
)
    J1 > 0 || throw(
        DomainError(J1, "J1J2Heisenberg1D Energy requires J1 > 0; got J1 = $J1.")
    )
    J2 ≥ 0 || throw(
        DomainError(J2, "J1J2Heisenberg1D Energy requires J2 ≥ 0; got J2 = $J2.")
    )
    j = J2 / J1
    if isapprox(j, 0.0; atol=1e-12)
        # j = 0 → pure Heisenberg.  Heisenberg1D exposes the Bethe-Hulthén
        # density via the legacy `GroundStateEnergyDensity` tag, with `J`
        # as a kwarg (the struct itself carries no parameters).
        return fetch(Heisenberg1D(), GroundStateEnergyDensity(), Infinite(); J=J1)
    elseif isapprox(j, 0.5; atol=1e-12)
        # j = 1/2 → Majumdar-Ghosh.  `MajumdarGhosh` is parameterised by
        # `J` as a struct field; route through its legacy
        # `GroundStateEnergyDensity` fetch.
        return fetch(MajumdarGhosh(; J=J1), GroundStateEnergyDensity(), Infinite())
    else
        throw(
            DomainError(
                j,
                "J1J2Heisenberg1D Energy: closed-form ground-state energy " *
                "known only at j = J2/J1 = 0 (pure Heisenberg, " *
                "Bethe-Hulthén 1938) or j = 1/2 (Majumdar-Ghosh 1969); " *
                "generic j is deferred to Phase 2 (DMRG, White-Affleck 1996). " *
                "Got j = $j.",
            ),
        )
    end
end
