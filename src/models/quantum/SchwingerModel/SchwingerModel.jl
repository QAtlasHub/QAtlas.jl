# ─────────────────────────────────────────────────────────────────────────────
# SchwingerModel — 1+1-D quantum electrodynamics (Schwinger 1962).
#
# Lagrangian (Lorentz form):
#
#     ℒ = -(1/4) F_{μν} F^{μν} + ψ̄ (i γ^μ ∂_μ - m - e γ^μ A_μ) ψ.
#
# In 1+1-D the gauge coupling `e` is dimensionful and the model is
# super-renormalisable.  At fermion mass `m = 0` Schwinger (1962)
# showed by abelian bosonisation that the spectrum is exactly a free
# massive scalar with the celebrated "Schwinger mass"
#
#     m_γ = e / √π.
#
# At nonzero mass the model is dual to a massive sine-Gordon model
# (Coleman-Jackiw-Susskind 1975), giving direct contact with sine-
# Gordon entries elsewhere in the atlas.
#
# Phase 1 registered the massless `m_γ = e/√π` mass gap.  Phase 2
# (#246) adds the anomaly-induced chiral condensate
# `⟨ψ̄ψ⟩ = −exp(γ_E)·e/(2π^{3/2})` (Coleman-Jackiw-Susskind 1975)
# via `ChiralCondensate`.  θ-vacuum structure and massive
# sine-Gordon duality remain tracked as Phase 3.
#
# References:
#   - J. Schwinger, [Schwinger1962](@cite).
#   - S. Coleman, R. Jackiw, L. Susskind, [ColemanJackiwSusskind1975](@cite).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Fermion bilinears c†c
#   Observable:  Fermion (number n = c†c, bilinear ⟨c†_i c_j⟩); derived spin observables follow spin S = σ/2
#   Reference:   docs/src/conventions.md §Fermion convention

"""
    SchwingerModel(; e::Real = 1.0, m::Real = 0.0) <: AbstractQAtlasModel

1+1-D quantum electrodynamics (Schwinger 1962) with gauge coupling
`e > 0` and fermion mass `m ≥ 0`.  At the massless point (`m = 0`)
the spectrum is exactly a free massive scalar with the Schwinger
mass `m_γ = e/√π`.

Quantities registered (Phases 1+2):

| Quantity                       | BC         | Method                              |
| ------------------------------ | ---------- | ----------------------------------- |
| [`MassGap`](@ref)              | `Infinite` | analytic (`e/√π`, massless only)    |
| [`ChiralCondensate`](@ref)    | `Infinite` | analytic (`-exp(γ_E)·e/(2π^{3/2})`, massless)  |

Massive `m > 0` Schwinger maps to massive sine-Gordon
(Coleman-Jackiw-Susskind 1975) and is tracked as Phase 2.

# References

- J. Schwinger, *Phys. Rev.* **128**, 2425 (1962).
- S. Coleman, R. Jackiw, L. Susskind, *Annals Phys.* **93**, 267 (1975).
"""
struct SchwingerModel <: AbstractQAtlasModel
    e::Float64
    m::Float64
end
SchwingerModel(; e::Real=1.0, m::Real=0.0) = SchwingerModel(Float64(e), Float64(m))

# ═══════════════════════════════════════════════════════════════════════════════
# Mass gap — massless Schwinger
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::SchwingerModel, ::MassGap, ::Infinite; e=m.e, m=m.m) -> Float64

Mass gap of the 1+1-D Schwinger model.  At the **massless point**
`m = 0` the celebrated Schwinger 1962 closed form applies:

    m_γ = e / √π.

`e ≤ 0` raises `DomainError`; `m ≠ 0` (massive Schwinger) raises
`DomainError` pointing to the massive-sine-Gordon Phase-2
implementation.

# References

- J. Schwinger, *Phys. Rev.* **128**, 2425 (1962).
"""
function fetch(
    model::SchwingerModel,
    ::MassGap,
    ::Infinite;
    e::Real=model.e,
    m::Real=model.m,
    kwargs...,
)
    e > 0 || throw(DomainError(e, "SchwingerModel MassGap requires e > 0; got e = $e."))
    iszero(m) || throw(
        DomainError(
            m,
            "SchwingerModel MassGap currently exposes only the massless point (m = 0); the massive case is dual to massive sine-Gordon (Coleman-Jackiw-Susskind 1975) and tracked as Phase 2.  Got m = $m.",
        ),
    )
    return e / sqrt(π)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Chiral condensate ⟨ψ̄ψ⟩ — anomaly-induced (massless Schwinger, Phase 2)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::SchwingerModel, ::ChiralCondensate, ::Infinite; e=model.e) -> Float64

Anomaly-induced chiral condensate in the **massless** Schwinger model
(Schwinger 1962; Coleman-Jackiw-Susskind 1975):

    ⟨ψ̄ψ⟩ = − exp(γ_E) · e / (2π^{3/2}),

where `γ_E ≈ 0.5772156649` is the Euler-Mascheroni constant.  Negative
by convention; magnitude scales linearly with the gauge coupling.
The non-zero condensate is the canonical 1+1-D example of anomaly-
induced spontaneous chiral-symmetry breaking: the classical
Lagrangian is chirally symmetric, but the axial U(1) anomaly forces
`⟨ψ̄ψ⟩ ≠ 0` in the quantum vacuum.

`e > 0` required.

# References

- J. Schwinger, *Phys. Rev.* **128**, 2425 (1962).
- S. Coleman, R. Jackiw, L. Susskind, *Annals Phys.* **93**, 267 (1975).
"""
function fetch(
    model::SchwingerModel, ::ChiralCondensate, ::Infinite; e::Real=model.e, kwargs...
)
    e > 0 ||
        throw(DomainError(e, "SchwingerModel ChiralCondensate requires e > 0; got e = $e."))
    return -exp(MathConstants.eulergamma) * e / (2 * π^(3/2))
end
