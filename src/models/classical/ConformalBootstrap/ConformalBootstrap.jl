# ─────────────────────────────────────────────────────────────────────────────
# ConformalBootstrap — the 3D Ising critical point pinned by the modular
# conformal bootstrap (Kos-Poland-Simmons-Duffin 2014 methodology;
# Kos-Poland-Simmons-Duffin-Vichi 2016 "Precision Islands" — precise values;
# Simmons-Duffin 2017 lightcone bootstrap cross-check).
#
# The 3D Ising universality class governs the liquid-vapour critical point,
# uniaxial ferromagnets, and innumerable physical systems.  The bootstrap
# fixes the scaling dimensions of the lowest scalar (σ) and energy (ε)
# primaries with record precision:
#
#     Δ_σ = 0.5181489(10)            (KPSD 2016, "Precision Islands")
#     Δ_ε = 1.412625(10)             (KPSD 2016, "Precision Islands")
#
# (KPSD 2014 introduced the mixed-correlator bootstrap method but
# reported wider error bars; the precise (10) uncertainties above are
# from KPSD 2016 / arXiv:1603.04436.  Cross-checked by Simmons-Duffin
# 2017 lightcone bootstrap.)
#
# These two numbers yield the canonical 3D Ising critical exponents
#
#     η = 2Δ_σ − (d − 2) = 2Δ_σ − 1  ≈ 0.03629,
#     ν = 1 / (d − Δ_ε) = 1 / (3 − Δ_ε) ≈ 0.62997.
#
# Phase 1 ships only the 3D Ising universality class; other
# bootstrap-pinned CFTs (O(N) vector models, stress-tensor and higher-spin
# multiplets, fermionic CFTs) are deferred to Phase 2.
#
# References:
#   - F. Kos, D. Poland, D. Simmons-Duffin,
#     "Bootstrapping the O(N) vector models", JHEP 06 (2014) 091
#     (methodology; wider error bars).
#   - F. Kos, D. Poland, D. Simmons-Duffin, A. Vichi,
#     "Precision Islands in the Ising and O(N) Models",
#     JHEP 08 (2016) 036; arXiv:1603.04436 (precise Δ_σ, Δ_ε).
#   - D. Simmons-Duffin,
#     "The Lightcone Bootstrap and the Spectrum of the 3D Ising CFT",
#     JHEP 03 (2017) 086 (lightcone-bootstrap cross-check).
#   - M. Reehorst, S. Rychkov, D. Simmons-Duffin, B. Sirois, N. Su,
#     B. van Rees, "Navigator Function for the Conformal Bootstrap",
#     SciPost Phys. 11 (2021) 072.
# ─────────────────────────────────────────────────────────────────────────────

"""
    ConformalBootstrap() <: AbstractQAtlasModel

The 3D Ising critical point with high-precision conformal-bootstrap
scaling dimensions.  Methodology: Kos-Poland-Simmons-Duffin 2014
(mixed-correlator bootstrap).  Precise values: KPSD-Vichi 2016
"Precision Islands" (arXiv:1603.04436).  Parameterless: the 3D
Ising universality class is unique.

Phase 1 ships only the 3D Ising universality class; other
bootstrap-pinned CFTs (O(N) vector models, stress-tensor and
higher-spin multiplets, fermionic CFTs) are deferred to Phase 2.

Quantities registered (Phase 1):

| Quantity                       | BC         | Method                                  |
| ------------------------------ | ---------- | --------------------------------------- |
| [`ConformalWeights`](@ref)     | `Infinite` | bootstrap reference (KPSD 2016 islands) |

Available kwargs for `ConformalWeights`:

  - `field = :σ` — lowest scalar primary, `Δ_σ = 0.51814894`.
  - `field = :ε` — lowest energy primary, `Δ_ε = 1.41262528`.

# References

- F. Kos, D. Poland, D. Simmons-Duffin, *JHEP* **06** (2014) 091
  — methodology (mixed-correlator bootstrap).
- F. Kos, D. Poland, D. Simmons-Duffin, A. Vichi, *JHEP* **08**
  (2016) 036; arXiv:1603.04436 — "Precision Islands", precise Δ_σ, Δ_ε.
- D. Simmons-Duffin, *JHEP* **03** (2017) 086 — lightcone-bootstrap
  cross-check.
- M. Reehorst, S. Rychkov, D. Simmons-Duffin, B. Sirois, N. Su,
  B. van Rees, *SciPost Phys.* **11** (2021) 072.
"""
struct ConformalBootstrap <: AbstractQAtlasModel end
# CONVENTION
#   Hamiltonian: see file-header description above
#   Observable:  per src/core/quantities.jl (matches the dispatch tag)
#   Reference:   docs/src/conventions.md (project-wide convention policy)
#   STATUS:      backfilled by PR (audit gate); per-field domain content
#                left to a follow-up - see issue tracker for the model-specific
#                Hamiltonian sign / observable normalisation.

# ═══════════════════════════════════════════════════════════════════════════════
# ConformalWeights — Δ_σ, Δ_ε bootstrap reference values
# (KPSD 2014 methodology; precise (10) values from KPSD-Vichi 2016
# "Precision Islands" / arXiv:1603.04436).
# ═══════════════════════════════════════════════════════════════════════════════

function fetch(
    ::ConformalBootstrap, ::ConformalWeights, ::Infinite; field::Symbol=:σ, kwargs...
)
    if field == :σ
        return 0.51814894
    elseif field == :ε
        return 1.41262528
    else
        throw(
            DomainError(
                field,
                "ConformalBootstrap ConformalWeights: Phase 1 exposes only field=:σ (spin) " *
                "and field=:ε (energy). Higher-spin and stress-tensor multiplets " *
                "(Simmons-Duffin 2017) deferred to Phase 2. Got field=:$field.",
            ),
        )
    end
end
