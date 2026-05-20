# ─────────────────────────────────────────────────────────────────────────────
# XYh1D — anisotropic XY chain in a transverse field (Lieb-Schultz-Mattis 1961).
#
# Hamiltonian (LSM 1961 convention, ferromagnetic sign for the XY exchanges):
#
#     H = -Σ_i ( J_x σ^x_i σ^x_{i+1} + J_y σ^y_i σ^y_{i+1} )
#         - h Σ_i σ^z_i,        J_x, J_y > 0,  h ∈ ℝ.
#
# Via the Jordan-Wigner transformation the model maps to a quadratic
# fermion problem with Bogoliubov dispersion
#
#     ε(k) = √[ (h − (J_x + J_y) cos k)² + (J_x − J_y)² sin² k ]   (≥ 0),
#
# and the single-particle gap is `2 ε(k_*)` minimised over k
# (Lieb-Schultz-Mattis 1961; Pfeuty 1970).
#
# Phase 1 of this file only exposes the **isotropic XX limit
# J_x = J_y = J**, where the (J_x − J_y) sin k branch vanishes and the
# dispersion collapses to ε(k) = |h − 2J cos k|, giving the closed-form
#
#     MassGap = 2 · max(0, |h| − 2J).
#
# - |h| <  2J  ⇒ gapless (Fermi surface inside the band, Luttinger liquid),
# - |h| =  2J  ⇒ critical (Lifshitz / BKT-like point),
# - |h| >  2J  ⇒ gapped paramagnetic phase.
#
# The general anisotropic case J_x ≠ J_y requires piecewise minimisation
# of the full dispersion (the band-bottom k_* depends on the
# (J_x, J_y, h) regime) and is deferred to Phase 2 — calling fetch on
# the anisotropic case raises DomainError.
#
# References:
#   - E. Lieb, T. Schultz, D. Mattis, "Two soluble models of an
#     antiferromagnetic chain", Annals Phys. 16, 407 (1961).
#   - P. Pfeuty, "The one-dimensional Ising model with a transverse
#     field", Annals Phys. 57, 79 (1970).
#   - S. Sachdev, *Quantum Phase Transitions* (2nd ed., CUP 2011), §2.
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Pauli σ (this file)
#   Observable:  Spin S = σ/2  (QAtlas-wide spin convention; see docs/src/conventions.md)

"""
    XYh1D(; Jx::Real = 1.0, Jy::Real = 1.0, h::Real = 0.0) <: AbstractQAtlasModel

Anisotropic XY chain in a transverse field (Lieb-Schultz-Mattis 1961):

    H = -Σ_i ( Jx σ^x_i σ^x_{i+1} + Jy σ^y_i σ^y_{i+1} ) - h Σ_i σ^z_i.

Requires `Jx > 0` and `Jy > 0`.

Quantities registered (Phase 1):

| Quantity               | BC         | Method                                       |
| ---------------------- | ---------- | -------------------------------------------- |
| [`MassGap`](@ref)      | `Infinite` | closed-form XX limit `2·max(0, |h| − 2J)`    |

Phase 1 only exposes the **isotropic XX limit** `Jx = Jy = J`; any
anisotropic `Jx ≠ Jy` raises `DomainError` (Phase 2: general-k
minimisation of the LSM/Pfeuty dispersion).

# References

- E. Lieb, T. Schultz, D. Mattis, *Annals Phys.* **16**, 407 (1961).
- P. Pfeuty, *Annals Phys.* **57**, 79 (1970).
- S. Sachdev, *Quantum Phase Transitions* (2nd ed., CUP 2011), §2.
"""
struct XYh1D <: AbstractQAtlasModel
    Jx::Float64
    Jy::Float64
    h::Float64
    function XYh1D(Jx::Real, Jy::Real, h::Real)
        Jx > 0 || throw(DomainError(Jx, "XYh1D requires Jx > 0; got Jx = $Jx."))
        Jy > 0 || throw(DomainError(Jy, "XYh1D requires Jy > 0; got Jy = $Jy."))
        return new(Float64(Jx), Float64(Jy), Float64(h))
    end
end
XYh1D(; Jx::Real=1.0, Jy::Real=1.0, h::Real=0.0) = XYh1D(Jx, Jy, h)

# ═══════════════════════════════════════════════════════════════════════════════
# Mass gap — closed-form XX limit (Jx = Jy)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::XYh1D, ::MassGap, ::Infinite;
          Jx=m.Jx, Jy=m.Jy, h=m.h, kwargs...) -> Float64

Single-particle Bogoliubov gap of the LSM/Pfeuty XY chain in a
transverse field, **restricted to the isotropic XX limit**
`Jx = Jy = J`.  In that limit the dispersion collapses to
ε(k) = |h − 2J cos k|, and the gap closes whenever `|h| ≤ 2J`:

    MassGap = 2 · max(0, |h| − 2J).

Anisotropic `Jx ≠ Jy` raises `DomainError` — Phase 2.

# References

- E. Lieb, T. Schultz, D. Mattis, *Annals Phys.* **16**, 407 (1961).
- P. Pfeuty, *Annals Phys.* **57**, 79 (1970).
"""
function fetch(
    m::XYh1D, ::MassGap, ::Infinite; Jx::Real=m.Jx, Jy::Real=m.Jy, h::Real=m.h, kwargs...
)
    Jx > 0 || throw(DomainError(Jx, "XYh1D MassGap requires Jx > 0; got Jx = $Jx."))
    Jy > 0 || throw(DomainError(Jy, "XYh1D MassGap requires Jy > 0; got Jy = $Jy."))
    if !isapprox(Jx, Jy; atol=1e-12)
        throw(
            DomainError(
                (Jx, Jy),
                "XYh1D MassGap: anisotropic case (Jx ≠ Jy) requires general-k " *
                "minimisation of √[(h-(Jx+Jy)cos k)² + (Jx-Jy)²sin²k] and is " *
                "deferred to Phase 2.  Phase 1 supports only the isotropic XX " *
                "limit Jx = Jy.  Got (Jx, Jy) = ($Jx, $Jy).",
            ),
        )
    end
    J = (Jx + Jy) / 2  # = Jx = Jy in the XX slice
    return 2 * max(0.0, abs(h) - 2J)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy per site — closed-form XX limit (Jx = Jy), any h
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::XYh1D, ::Energy{:per_site}, ::Infinite;
          Jx=m.Jx, Jy=m.Jy, h=m.h, kwargs...) -> Float64

Ground-state energy per site of the LSM/Pfeuty XY chain in a
transverse field, **restricted to the isotropic XX limit**
`Jx = Jy = J`.  Via Jordan-Wigner the model maps to free spinless
fermions with single-particle dispersion `ξ(k) = 2h − 4J cos k`; the
ground state fills all modes with `ξ(k) < 0`.

In the thermodynamic limit, with the LSM 1961 sign convention
`H = −J Σ(σˣσˣ + σʸσʸ) − h Σ σᶻ`,

* `h ≥ 2J`   (empty Fermi sea, all spins ↑): `E/N = −h`
* `h ≤ −2J`  (full Fermi sea, all spins ↓): `E/N =  h`
* `|h| < 2J` (partially filled, `k_F = arccos(h/2J)`):

      E/N = −h + (2h/π) · arccos(h/(2J))
                − (4J/π) · √(1 − (h/(2J))²).

At `h = 0`: `E/N = −4J/π ≈ −1.27323954J` (Lieb-Schultz-Mattis 1961).

Anisotropic `Jx ≠ Jy` raises `DomainError` — Phase 2.

# References

- E. Lieb, T. Schultz, D. Mattis, *Annals Phys.* **16**, 407 (1961).
- P. Pfeuty, *Annals Phys.* **57**, 79 (1970).
"""
function fetch(
    m::XYh1D,
    ::Energy{:per_site},
    ::Infinite;
    Jx::Real=m.Jx,
    Jy::Real=m.Jy,
    h::Real=m.h,
    kwargs...,
)
    Jx > 0 ||
        throw(DomainError(Jx, "XYh1D Energy{:per_site} requires Jx > 0; got Jx = $Jx."))
    Jy > 0 ||
        throw(DomainError(Jy, "XYh1D Energy{:per_site} requires Jy > 0; got Jy = $Jy."))
    if !isapprox(Jx, Jy; atol=1e-12)
        throw(
            DomainError(
                (Jx, Jy),
                "XYh1D Energy{:per_site}: anisotropic case (Jx ≠ Jy) requires " *
                "integration of the full Bogoliubov dispersion " *
                "√[(h-(Jx+Jy)cos k)² + (Jx-Jy)²sin²k] and is deferred to Phase 2. " *
                "Phase 1 supports only the isotropic XX limit Jx = Jy.  " *
                "Got (Jx, Jy) = ($Jx, $Jy).",
            ),
        )
    end
    J = (Jx + Jy) / 2  # = Jx = Jy in the XX slice
    if h >= 2J
        return -float(h)
    elseif h <= -2J
        return float(h)
    else
        x = h / (2J)
        return -h + (2h / pi) * acos(x) - (4J / pi) * sqrt(1 - x * x)
    end
end
