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
#   - R. J. Baxter, Annals Phys. 70, 193 (1972).
#   - L. D. Faddeev, L. A. Takhtajan, J. Soviet Math. 24, 241 (1984).
# ─────────────────────────────────────────────────────────────────────────────

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
    HeisenbergXYZ(Float64(Jx), Float64(Jy), Float64(Jz))
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
    Jx == Jy || throw(
        DomainError(
            (Jx, Jy),
            "HeisenbergXYZ Energy(:per_site) currently supports only the axis-aligned reduction Jx = Jy (delegated to XXZ1D); general XYZ via Baxter elliptic Bethe ansatz is tracked as Phase 2.  Got (Jx, Jy) = ($Jx, $Jy).",
        ),
    )
    # Restrict to Jx > 0 so that Delta = Jz/Jx has unambiguous sign and the
    # delegation lands in the XXZ1D-tested domain.  Jx < 0 (FM exchange) would
    # flip Delta's sign into XXZ1D regions whose Yang-Yang analytic
    # continuation is not exercised by the current XXZ1D test suite.
    Jx > 0 || throw(
        DomainError(
            Jx,
            "HeisenbergXYZ Energy(:per_site): Jx > 0 required (Jx == 0 is Ising-like; Jx < 0 is the FM-exchange XXZ regime whose XXZ1D delegation is not yet tested).  Got Jx = $Jx.",
        ),
    )
    Δ = Jz / Jx
    return QAtlas.fetch(QAtlas.XXZ1D(; J=Jx, Δ=Δ), Energy{:per_site}(), Infinite())
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
    if !(isapprox(Jx, Jy; atol=1e-12) && isapprox(Jy, Jz; atol=1e-12))
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
