# ─────────────────────────────────────────────────────────────────────────────
# HeisenbergXYZ — spin-½ chain with three independent exchange couplings.
#
# Hamiltonian:
#
#     H = Σ_i [ Jx Sˣᵢ Sˣᵢ₊₁ + Jy Sʸᵢ Sʸᵢ₊₁ + Jz Sᶻᵢ Sᶻᵢ₊₁ ],   S = 1/2.
#
# This is the most general 1D nearest-neighbour spin-½ chain still
# admitting Bethe-ansatz integrability — Baxter (1972) showed it is
# the row-to-row transfer matrix of the 8-vertex model whose ground
# state is expressed in elliptic theta functions.
#
# The general closed form involves elliptic integrals and is deferred
# to a later phase.  What this file *does* expose are the canonical
# *axis-aligned reductions*, each delegated to an already-implemented
# Atlas entry:
#
#   * Jx = Jy   →  XXZ1D(J = Jx, Δ = Jz / Jx) at Infinite
#                  (covers Heisenberg AF/FM at Jz = Jx, XX at Jz = 0,
#                  Ising-like Jz/Jx > 1 etc. — already handled by
#                  XXZ.jl + Yang-Yang single integral and the
#                  closed-form points)
#
# All other (Jx ≠ Jy) triples raise DomainError pointing to the
# Baxter-elliptic phase 2 implementation.
#
# References:
#   - R. J. Baxter, [Baxter1972a](@cite).
#   - L. D. Faddeev, L. A. Takhtajan, [FaddeevTakhtajan1984](@cite).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Spin S (this file)
#   Observable:  Spin S         (QAtlas-wide spin convention; see docs/src/conventions.md)

"""
    HeisenbergXYZ(; Jx::Real = 1.0, Jy::Real = 1.0, Jz::Real = 1.0)
        <: AbstractQAtlasModel

Spin-½ XYZ chain with three independent exchange couplings

    H = Σ_i [ Jx Sˣᵢ Sˣᵢ₊₁ + Jy Sʸᵢ Sʸᵢ₊₁ + Jz Sᶻᵢ Sᶻᵢ₊₁ ].

Most-general 1-D nearest-neighbour spin-½ integrable model
(Baxter 1972).  This release exposes only the *axis-aligned reductions*
`Jx = Jy` by delegating to the existing [`XXZ1D`](@ref) machinery:

    HeisenbergXYZ(Jx = J, Jy = J, Jz)   ≡   XXZ1D(J = J, Δ = Jz / J).

General XYZ ground-state energy (Baxter elliptic Bethe ansatz)
is tracked as a follow-up phase and raises `DomainError` here.

Quantities registered:

| Quantity                       | BC         | Method                 |
| ------------------------------ | ---------- | ---------------------- |
| [`Energy`](@ref) (`:per_site`) | `Infinite` | delegated to XXZ1D     |

# References

- R. J. Baxter, *Annals Phys.* **70**, 193 (1972).
- L. D. Faddeev, L. A. Takhtajan, *J. Soviet Math.* **24**, 241 (1984).
"""
struct HeisenbergXYZ <: AbstractQAtlasModel
    Jx::Float64
    Jy::Float64
    Jz::Float64
end
function HeisenbergXYZ(; Jx::Real=1.0, Jy::Real=1.0, Jz::Real=1.0)
    return HeisenbergXYZ(Float64(Jx), Float64(Jy), Float64(Jz))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Ground-state energy per site — axis-aligned (Jx = Jy) reduction
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::HeisenbergXYZ, ::Energy{:per_site}, ::Infinite;
          Jx=m.Jx, Jy=m.Jy, Jz=m.Jz, kwargs...) -> Float64

Ground-state energy per site of the spin-½ XYZ chain in the
thermodynamic limit, restricted to the axis-aligned reduction
`Jx = Jy`:

    HeisenbergXYZ(Jx = J, Jy = J, Jz)   ⟶   fetch(XXZ1D(J = J, Δ = Jz/J), …).

This routes the call through the XXZ1D Yang-Yang single integral
(general -1 < Δ < 1) and the three closed-form points Δ ∈ {-1, 0, 1}
already implemented in [`XXZ1D`](@ref).  General (Jx ≠ Jy) triples
require the Baxter elliptic Bethe ansatz and currently raise
`DomainError` — Phase 2.

`Jx > 0` is required (so `Delta = Jz/Jx` is well-defined and the delegation routes into the XXZ1D-tested domain; FM-exchange `Jx < 0` is Phase 2).

# References

- C. N. Yang, C. P. Yang, *Phys. Rev.* **150**, 327 (1966).
- R. J. Baxter, *Annals Phys.* **70**, 193 (1972).
"""
function fetch(
    m::HeisenbergXYZ,
    ::Energy{:per_site},
    ::Infinite;
    Jx::Real=m.Jx,
    Jy::Real=m.Jy,
    Jz::Real=m.Jz,
    kwargs...,
)
    if Jx == Jy
        Jx > 0 || throw(
            DomainError(
                Jx,
                "HeisenbergXYZ Energy(:per_site): Jx > 0 required in the axial Jx = Jy " *
                "branch (Jx == 0 is Ising-like; Jx < 0 is the FM-exchange XXZ regime " *
                "whose XXZ1D delegation is not yet tested). Got Jx = $Jx.",
            ),
        )
        Δ = Jz / Jx
        return QAtlas.fetch(QAtlas.XXZ1D(; J=Jx, Δ=Δ), Energy{:per_site}(), Infinite())
    elseif Jz == 0
        # XY anisotropic line (Jx != Jy, Jz = 0). Lieb-Schultz-Mattis closed form.
        return _heisenberg_xyz_gs_energy_xy_line(Jx, Jy)
    else
        throw(
            DomainError(
                (Jx, Jy, Jz),
                "HeisenbergXYZ Energy(:per_site): generic XYZ (Jx != Jy, Jz != 0) " *
                "requires the Baxter (1972) elliptic Bethe-ansatz machinery, deferred " *
                "to Phase 3. Supported now: Jx = Jy (XXZ delegation) and Jz = 0 (XY line). " *
                "Got (Jx, Jy, Jz) = ($Jx, $Jy, $Jz).",
            ),
        )
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Luttinger parameter at the isotropic Heisenberg point (Phase 2)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::HeisenbergXYZ, ::LuttingerParameter, ::Infinite;
          Jx=m.Jx, Jy=m.Jy, Jz=m.Jz) -> Float64

Luttinger-liquid parameter at the **isotropic Heisenberg point**
`Jx = Jy = Jz`, delegated to `XXZ1D(Δ = 1)`:

    K = 1/2          (SU(2)-symmetric AFM, Luther-Peschel 1975)

The delegation chain is `HeisenbergXYZ → XXZ1D(Δ=1)` (matching the
`Energy(:per_site)` reduction).  Once a dedicated
`fetch(::Heisenberg1D, ::LuttingerParameter, ::Infinite)` lands on
main (tracked separately), the intermediate `Heisenberg1D` hop can
be reinstated; the final `K` is unchanged.

For `Jx = Jy ≠ Jz` (XXZ axial anisotropy), use `XXZ1D` directly.
For generic XYZ (Baxter 1972 elliptic / theta-function regime),
defer to a later phase — this Phase-2 path throws `DomainError`
for non-isotropic couplings.

# References

- A. Luther, I. Peschel, *Phys. Rev. B* **12**, 3908 (1975).
- R. J. Baxter, *Ann. Phys.* **70**, 193 (1972) — elliptic XYZ.
"""
function fetch(
    m::HeisenbergXYZ,
    ::LuttingerParameter,
    ::Infinite;
    Jx::Real=m.Jx,
    Jy::Real=m.Jy,
    Jz::Real=m.Jz,
    kwargs...,
)
    if !(Jx == Jy == Jz)
        throw(
            DomainError(
                (Jx, Jy, Jz),
                "HeisenbergXYZ LuttingerParameter: Phase 2 supports only the isotropic " *
                "Heisenberg point Jx = Jy = Jz. Generic XYZ requires elliptic " *
                "(Baxter 1972) machinery, deferred. Got (Jx, Jy, Jz) = ($Jx, $Jy, $Jz).",
            ),
        )
    end
    # Delegate to XXZ1D(Δ=1) — same target as the Heisenberg1D path
    # (Heisenberg1D itself delegates to XXZ1D(Δ=1) once PR #347 lands).
    return QAtlas.fetch(QAtlas.XXZ1D(; J=1.0, Δ=1.0), LuttingerParameter(), Infinite())
end

# ==============================================================================
# Ground-state energy density (Phase 2: XY line + XXZ delegation)
# ==============================================================================

"""
    fetch(m::HeisenbergXYZ, ::GroundStateEnergyDensity, ::Infinite;
          Jx=m.Jx, Jy=m.Jy, Jz=m.Jz) -> Float64

Ground-state energy density (per site) of the spin-1/2 XYZ chain at the
thermodynamic limit. Supported parameter regimes:

| Regime              | Method                                                  |
|---------------------|---------------------------------------------------------|
| `Jx = Jy`           | delegated to `XXZ1D(J=Jx, D=Jz/Jx)` `Energy(:per_site)` |
| `Jz = 0` (XY line)  | closed-form Lieb-Schultz-Mattis (1961) free-fermion     |
| generic XYZ         | `DomainError` (Baxter 1972 elliptic deferred to Phase 3)|

# References

- E. Lieb, T. Schultz, D. Mattis, *Ann. Phys.* **16**, 407 (1961).
- C. N. Yang, C. P. Yang, *Phys. Rev.* **150**, 327 (1966).
- R. J. Baxter, *Ann. Phys.* **70**, 193 (1972) -- generic XYZ via elliptic Bethe ansatz.
"""
function fetch(
    m::HeisenbergXYZ,
    ::GroundStateEnergyDensity,
    ::Infinite;
    Jx::Real=m.Jx,
    Jy::Real=m.Jy,
    Jz::Real=m.Jz,
    kwargs...,
)
    if Jx == Jy
        # Axial XXZ reduction. Reuse the Energy(:per_site) target since XXZ1D
        # exposes it (its GS-energy semantics agree at Infinite).
        Jx > 0 || throw(
            DomainError(
                Jx,
                "HeisenbergXYZ GroundStateEnergyDensity at Jx = Jy: Jx > 0 required " *
                "(Jx = 0 is Ising-like; Jx < 0 is the FM-exchange XXZ regime whose " *
                "XXZ1D delegation is not yet tested). Got Jx = $Jx.",
            ),
        )
        Δ = Jz / Jx
        return QAtlas.fetch(QAtlas.XXZ1D(; J=Jx, Δ=Δ), Energy{:per_site}(), Infinite())
    elseif Jz == 0
        # XY anisotropic line. Lieb-Schultz-Mattis closed form.
        return _heisenberg_xyz_gs_energy_xy_line(Jx, Jy)
    else
        throw(
            DomainError(
                (Jx, Jy, Jz),
                "HeisenbergXYZ GroundStateEnergyDensity: generic XYZ (Jx != Jy, Jz != 0) " *
                "requires the Baxter (1972) elliptic Bethe-ansatz machinery, deferred " *
                "to Phase 3. Supported now: Jx = Jy (XXZ delegation) and Jz = 0 (XY line). " *
                "Got (Jx, Jy, Jz) = ($Jx, $Jy, $Jz).",
            ),
        )
    end
end

# ==============================================================================
# Mass gap (Phase 2: critical axial + XY line free-fermion)
# ==============================================================================

"""
    fetch(m::HeisenbergXYZ, ::MassGap, ::Infinite;
          Jx=m.Jx, Jy=m.Jy, Jz=m.Jz) -> Float64

Single-particle mass gap (smallest excitation energy above the ground
state) of the spin-1/2 XYZ chain at the thermodynamic limit.

| Regime                              | Method                                         |
|-------------------------------------|------------------------------------------------|
| `Jx = Jy`, `abs(Jz/Jx) <= 1`        | gapless critical XXZ: returns `0.0`            |
| `Jz = 0` (XY line, `Jx != Jy`)      | LSM dispersion minimum: `(1/4)*abs(Jx - Jy)`   |
| massive axial AFM (`Jz/Jx > 1`)     | `DomainError` (Yang-Yang gap deferred)         |
| generic XYZ (`Jx != Jy`, `Jz != 0`) | `DomainError` (Baxter elliptic Phase 3)        |

# References

- E. Lieb, T. Schultz, D. Mattis, *Ann. Phys.* **16**, 407 (1961).
- C. N. Yang, C. P. Yang, *Phys. Rev.* **150**, 327 (1966).
- A. Luther, I. Peschel, *Phys. Rev. B* **12**, 3908 (1975).
- R. J. Baxter, *Ann. Phys.* **70**, 193 (1972).
"""
function fetch(
    m::HeisenbergXYZ,
    ::MassGap,
    ::Infinite;
    Jx::Real=m.Jx,
    Jy::Real=m.Jy,
    Jz::Real=m.Jz,
    kwargs...,
)
    if Jx == Jy
        Jx > 0 || throw(
            DomainError(
                Jx, "HeisenbergXYZ MassGap axial branch requires Jx > 0; got Jx = $Jx."
            ),
        )
        if abs(Jz / Jx) <= 1
            return 0.0  # critical XXZ -1 <= Delta <= 1: gapless spinon spectrum
        else
            throw(
                DomainError(
                    Jz / Jx,
                    "HeisenbergXYZ MassGap at Jx = Jy: massive AFM regime " *
                    "Delta = Jz/Jx > 1 not yet implemented (Yang-Yang gap " *
                    "formula deferred); got Delta = $(Jz/Jx).",
                ),
            )
        end
    elseif Jz == 0
        # XY anisotropic line: single-fermion gap min_k epsilon(k) with
        #   epsilon(k) = (1/4) * sqrt((Jx+Jy)^2 cos^2(k) + (Jx-Jy)^2 sin^2(k))
        # Minimum at k = +- pi/2 gives (1/4)*|Jx - Jy|.
        return abs(Jx - Jy) / 4
    else
        throw(
            DomainError(
                (Jx, Jy, Jz),
                "HeisenbergXYZ MassGap: generic XYZ (Jx != Jy, Jz != 0) requires " *
                "the Baxter (1972) / Johnson-Krinsky-McCoy (1973) elliptic mass " *
                "spectrum, deferred to Phase 3. Got (Jx, Jy, Jz) = ($Jx, $Jy, $Jz).",
            ),
        )
    end
end

# ==============================================================================
# Correlation length (Phase 2: critical axial + XY line free-fermion)
# ==============================================================================

"""
    fetch(m::HeisenbergXYZ, ::CorrelationLength, ::Infinite;
          Jx=m.Jx, Jy=m.Jy, Jz=m.Jz) -> Float64

Asymptotic correlation length of the spin-1/2 XYZ chain at the
thermodynamic limit.

| Regime                              | Method                                                              |
|-------------------------------------|---------------------------------------------------------------------|
| `Jx = Jy`, `abs(Jz/Jx) <= 1`        | gapless critical XXZ: returns `Inf`                                 |
| `Jz = 0` (XY line, `Jx != Jy`)      | `xi = 1/asinh(abs(Jx-Jy) / (2*sqrt(Jx*Jy)))` (dispersion zero)      |
| massive axial AFM (`Jz/Jx > 1`)     | `DomainError`                                                       |
| generic XYZ                         | `DomainError` (Baxter elliptic Phase 3)                             |

# References

- E. Lieb, T. Schultz, D. Mattis, *Ann. Phys.* **16**, 407 (1961).
- B. M. McCoy, T. T. Wu, *Phys. Rev.* **174**, 546 (1968).
"""
function fetch(
    m::HeisenbergXYZ,
    ::CorrelationLength,
    ::Infinite;
    Jx::Real=m.Jx,
    Jy::Real=m.Jy,
    Jz::Real=m.Jz,
    kwargs...,
)
    if Jx == Jy
        Jx > 0 || throw(
            DomainError(
                Jx,
                "HeisenbergXYZ CorrelationLength axial branch requires Jx > 0; got Jx = $Jx.",
            ),
        )
        if abs(Jz / Jx) <= 1
            return Inf
        else
            throw(
                DomainError(
                    Jz / Jx,
                    "HeisenbergXYZ CorrelationLength at Jx = Jy: massive AFM " *
                    "Delta = Jz/Jx > 1 not yet implemented; got Delta = $(Jz/Jx).",
                ),
            )
        end
    elseif Jz == 0
        prod = Jx * Jy
        prod > 0 || throw(
            DomainError(
                (Jx, Jy),
                "HeisenbergXYZ CorrelationLength on XY line requires Jx Jy > 0; " *
                "got Jx*Jy = $prod.",
            ),
        )
        return 1 / asinh(abs(Jx - Jy) / (2 * sqrt(prod)))
    else
        throw(
            DomainError(
                (Jx, Jy, Jz),
                "HeisenbergXYZ CorrelationLength: generic XYZ requires Baxter " *
                "(1972) elliptic correlation, deferred to Phase 3. Got (Jx, Jy, Jz) = ($Jx, $Jy, $Jz).",
            ),
        )
    end
end
# ==============================================================================
# Spontaneous magnetization (Phase 2: XY anisotropic line)
# ==============================================================================

"""
    fetch(m::HeisenbergXYZ, ::SpontaneousMagnetization, ::Infinite;
          Jx=m.Jx, Jy=m.Jy, Jz=m.Jz) -> Float64

Spontaneous magnetization (sublattice order parameter) of the spin-1/2
XYZ chain at the thermodynamic limit.

| Regime                              | Method                                                            |
|-------------------------------------|-------------------------------------------------------------------|
| `Jx = Jy`, `abs(Jz/Jx) <= 1`        | gapless critical XXZ, no LRO in 1D: returns `0.0`                 |
| `Jz = 0` (XY line, `Jx != Jy`)      | McCoy-Wu / Pfeuty 1/8: `M = (1 - (Jmin/Jmax)^2)^(1/8)`            |
| massive axial AFM (`Jz/Jx > 1`)     | `DomainError` (Neel staggered M Baxter 1971 deferred)             |
| generic XYZ                         | `DomainError` (Baxter 1973 PRL elliptic theta deferred)           |

# References

- B. M. McCoy, T. T. Wu, *Phys. Rev.* **174**, 546 (1968).
- P. Pfeuty, *Ann. Phys.* **57**, 79 (1970).
- R. J. Baxter, *Phys. Rev. Lett.* **31**, 1294 (1973) -- generic XYZ.
"""
function fetch(
    m::HeisenbergXYZ,
    ::SpontaneousMagnetization,
    ::Infinite;
    Jx::Real=m.Jx,
    Jy::Real=m.Jy,
    Jz::Real=m.Jz,
    kwargs...,
)
    if Jx == Jy
        Jx > 0 || throw(
            DomainError(
                Jx,
                "HeisenbergXYZ SpontaneousMagnetization axial branch requires Jx > 0; got Jx = $Jx.",
            ),
        )
        if abs(Jz / Jx) <= 1
            return 0.0
        else
            throw(
                DomainError(
                    Jz / Jx,
                    "HeisenbergXYZ SpontaneousMagnetization at Jx = Jy: massive AFM " *
                    "(Neel staggered M) deferred to Phase 3; got Delta = $(Jz/Jx).",
                ),
            )
        end
    elseif Jz == 0
        ax, ay = abs(Jx), abs(Jy)
        (ax > 0 && ay > 0) || throw(
            DomainError(
                (Jx, Jy),
                "HeisenbergXYZ SpontaneousMagnetization on XY line requires both " *
                "|Jx| > 0 and |Jy| > 0; got (Jx, Jy) = ($Jx, $Jy).",
            ),
        )
        jmin, jmax = minmax(ax, ay)
        return (1 - (jmin / jmax)^2)^(1 / 8)
    else
        throw(
            DomainError(
                (Jx, Jy, Jz),
                "HeisenbergXYZ SpontaneousMagnetization: generic XYZ requires Baxter " *
                "(1973) elliptic theta-function staggered polarization, deferred " *
                "to Phase 3. Got (Jx, Jy, Jz) = ($Jx, $Jy, $Jz).",
            ),
        )
    end
end
