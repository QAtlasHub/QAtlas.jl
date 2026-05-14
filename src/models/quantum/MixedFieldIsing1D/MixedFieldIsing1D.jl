# ─────────────────────────────────────────────────────────────────────────────
# Mixed-Field Ising Model — 1D chain with transverse AND longitudinal field
#
# Hamiltonian:
#   H = -J Σᵢ σᶻᵢσᶻᵢ₊₁  -  h_x Σᵢ σˣᵢ  -  h_z Σᵢ σᶻᵢ,    J > 0.
#
# Two qualitatively different regimes:
#
#   * h_z = 0:  recovers the standard transverse-field Ising model (TFIM),
#               Jordan-Wigner solvable, with closed-form mass gap
#               Δ_∞(J, h_x) = 2|h_x - J| (Pfeuty 1970).
#
#   * h_z ≠ 0:  the model is NON-INTEGRABLE and is the canonical minimal
#               non-integrable Ising chain used throughout the ETH /
#               thermalisation / chaos literature (McCoy-Wu 1978;
#               Banuls-Cirac-Hastings 2011).  No closed form for the gap.
#
# Phase 1 scope (this file):
#   * h_z = 0   → delegate the (Δ, Infinite) fetch to TFIM (closed form).
#   * h_z ≠ 0   → throw DomainError, deferring to Phase 2 (numerical ED /
#                 DMRG, to be implemented in a follow-up).
#
# References:
#   - Pfeuty, Ann. Phys. 57, 79 (1970)
#   - McCoy & Wu, Phys. Rev. B 18, 4886 (1978)
#   - Banuls, Cirac & Hastings, Phys. Rev. Lett. 106, 050405 (2011)
# ─────────────────────────────────────────────────────────────────────────────

"""
    MixedFieldIsing1D(; J = 1.0, h_x = 1.0, h_z = 0.0) <: AbstractQAtlasModel

The 1D mixed-field Ising chain

    H = -J Σ_i σᶻ_i σᶻ_{i+1} - h_x Σ_i σˣ_i - h_z Σ_i σᶻ_i,    J > 0.

The default `h_z = 0.0` sits on the Phase-1 closed-form point where the
model coincides with the standard TFIM (`Δ_∞ = 2|h_x − J|`, gap closing
at the quantum critical point `h_x = J`).

At `h_z ≠ 0` the model is **non-integrable** — fetch methods will throw
`DomainError` in Phase 1 and route to numerical solvers in Phase 2.
"""
struct MixedFieldIsing1D <: AbstractQAtlasModel
    J::Float64
    h_x::Float64
    h_z::Float64
    function MixedFieldIsing1D(J::Real, h_x::Real, h_z::Real)
        J > 0 || throw(DomainError(J, "MixedFieldIsing1D requires J > 0; got J = $J."))
        return new(Float64(J), Float64(h_x), Float64(h_z))
    end
end
function MixedFieldIsing1D(; J::Real=1.0, h_x::Real=1.0, h_z::Real=0.0)
    MixedFieldIsing1D(J, h_x, h_z)
end

# ═══════════════════════════════════════════════════════════════════════════════
# MassGap, Infinite — Phase 1: delegate at h_z = 0, DomainError otherwise.
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::MixedFieldIsing1D, ::MassGap, ::Infinite; J, h_x, h_z) -> Float64

Mass gap of the infinite mixed-field Ising chain.

* At `h_z = 0` (default): delegate to `TFIM(; J=J, h=h_x)` and return the
  closed-form `Δ = 2|h_x − J|` (Pfeuty 1970).
* At `h_z ≠ 0`: throw `DomainError` — the model is non-integrable
  (canonical ETH / thermalisation benchmark, McCoy-Wu 1978;
  Banuls-Cirac-Hastings 2011) and no closed-form gap exists.  Numerical
  routes (ED / DMRG) will land in Phase 2.

`J, h_x, h_z` may be overridden via keyword to evaluate the gap at a
point different from the struct's stored parameters without rebuilding
the model.
"""
function fetch(
    m::MixedFieldIsing1D,
    ::MassGap,
    ::Infinite;
    J::Real=m.J,
    h_x::Real=m.h_x,
    h_z::Real=m.h_z,
    kwargs...,
)
    J > 0 || throw(DomainError(J, "MixedFieldIsing1D MassGap requires J > 0; got J = $J."))
    if !iszero(h_z)
        throw(
            DomainError(
                h_z,
                "MixedFieldIsing1D MassGap: the model is non-integrable at h_z ≠ 0 " *
                "(no closed-form ground-state gap). Phase 1 supports only the " *
                "h_z = 0 limit (transverse-field Ising, delegated to TFIM). " *
                "Got h_z = $h_z. Numerical ED/DMRG routes will land in Phase 2.",
            ),
        )
    end
    # Phase-1 delegate: TFIM's MassGap at Infinite is read from its struct
    # fields (no J/h fetch-kwargs), so we build a TFIM at (J, h_x) and call.
    return fetch(TFIM(; J=J, h=h_x), MassGap(), Infinite())
end
