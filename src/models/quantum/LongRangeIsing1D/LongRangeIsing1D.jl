# ─────────────────────────────────────────────────────────────────────────────
# LongRangeIsing1D — 1D transverse-field Ising chain with power-law-decaying
# interactions.
#
# Hamiltonian:
#
#     H = -J Σ_{i<j} σᶻᵢ σᶻⱼ / |i-j|^α  -  h Σᵢ σˣᵢ,   J > 0, h > 0, α > 0.
#
# Two well-known closed-form limits:
#
#   - α → ∞  : nearest-neighbour TFIM   →  Δ = 2|h - J|  (Pfeuty 1970)
#   - α → 0  : Lipkin-Meshkov-Glick collective-spin model (S = N/2 sector)
#
# Finite α (between these limits) realises the long-range Ising phase
# diagram studied numerically by DMRG / VMC:
#
#   - L. Koffel, M. Lewenstein, L. Tagliacozzo,
#       *Phys. Rev. Lett.* **109**, 267203 (2012).
#   - Z.-X. Gong, M. Foss-Feig, S. Michalakis, A. V. Gorshkov,
#       *Phys. Rev. Lett.* **113**, 030601 (2014).
#
# Phase 1 of this file exposes only the α = Inf NN-TFIM limit, delegated
# to the existing exactly-solvable [`TFIM`](@ref) entry.  Any finite α
# raises `DomainError` pointing to the DMRG / VMC phase-2 follow-up.
#
# References:
#   - P. Pfeuty, *Annals Phys.* **57**, 79 (1970).
#   - L. Koffel, M. Lewenstein, L. Tagliacozzo,
#       *Phys. Rev. Lett.* **109**, 267203 (2012).
#   - Z.-X. Gong, M. Foss-Feig, S. Michalakis, A. V. Gorshkov,
#       *Phys. Rev. Lett.* **113**, 030601 (2014).
# ─────────────────────────────────────────────────────────────────────────────

"""
    LongRangeIsing1D(; J::Real = 1.0, h::Real = 1.0, α::Real = Inf)
        <: AbstractQAtlasModel

1D transverse-field Ising chain with power-law-decaying couplings:

    H = -J Σ_{i<j} σᶻᵢ σᶻⱼ / |i-j|^α  -  h Σᵢ σˣᵢ.

`J > 0` (ferromagnetic), `h > 0` (transverse field), `α > 0` (decay
exponent).  The default `α = Inf` sits at the nearest-neighbour TFIM
point, for which the gap is known in closed form.

Quantities registered (Phase 1):

| Quantity              | BC         | Method                                         |
| --------------------- | ---------- | ---------------------------------------------- |
| [`MassGap`](@ref)     | `Infinite` | delegated to [`TFIM`](@ref) at α = Inf         |

Any finite `α` raises `DomainError` — Phase 2 will add DMRG / VMC
support for the long-range phase diagram (Koffel-Lewenstein-Tagliacozzo
2012, Gong-Foss-Feig-Michalakis-Gorshkov 2014).

# References

- P. Pfeuty, *Annals Phys.* **57**, 79 (1970).
- L. Koffel, M. Lewenstein, L. Tagliacozzo,
  *Phys. Rev. Lett.* **109**, 267203 (2012).
- Z.-X. Gong, M. Foss-Feig, S. Michalakis, A. V. Gorshkov,
  *Phys. Rev. Lett.* **113**, 030601 (2014).
"""
struct LongRangeIsing1D <: AbstractQAtlasModel
    J::Float64
    h::Float64
    α::Float64
    function LongRangeIsing1D(J::Real, h::Real, α::Real)
        J > 0 || throw(DomainError(J, "LongRangeIsing1D requires J > 0; got J = $J."))
        h > 0 || throw(DomainError(h, "LongRangeIsing1D requires h > 0; got h = $h."))
        α > 0 || throw(DomainError(α, "LongRangeIsing1D requires α > 0; got α = $α."))
        return new(Float64(J), Float64(h), Float64(α))
    end
end
LongRangeIsing1D(; J::Real=1.0, h::Real=1.0, α::Real=Inf) = LongRangeIsing1D(J, h, α)

# ═══════════════════════════════════════════════════════════════════════════════
# Mass gap — α = Inf delegation to TFIM (nearest-neighbour limit)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::LongRangeIsing1D, ::MassGap, ::Infinite;
          J=m.J, h=m.h, α=m.α, kwargs...) -> Float64

Single-quasiparticle gap of the long-range TFIM at the **α = Inf**
nearest-neighbour limit.  Internally constructs `TFIM(; J=J, h=h)` and
forwards.  Closed-form result (Pfeuty 1970):

    Δ = 2 |h - J|.

Finite `α` raises `DomainError` — the long-range phase diagram
(Koffel-Lewenstein-Tagliacozzo 2012, Gong-Foss-Feig 2014) requires
numerical DMRG / VMC and is deferred to Phase 2.
"""
function fetch(
    m::LongRangeIsing1D,
    ::MassGap,
    ::Infinite;
    J::Real=m.J,
    h::Real=m.h,
    α::Real=m.α,
    kwargs...,
)
    J > 0 || throw(DomainError(J, "LongRangeIsing1D MassGap requires J > 0; got J = $J."))
    h > 0 || throw(DomainError(h, "LongRangeIsing1D MassGap requires h > 0; got h = $h."))
    if !isinf(α)
        throw(
            DomainError(
                α,
                "LongRangeIsing1D MassGap: closed-form supported only at α = Inf " *
                "(nearest-neighbour TFIM limit, delegated to TFIM). Finite α " *
                "requires numerical DMRG / VMC (Koffel-Lewenstein-Tagliacozzo 2012, " *
                "Gong-Foss-Feig 2014) and is deferred to Phase 2. Got α = $α.",
            ),
        )
    end
    return QAtlas.fetch(QAtlas.TFIM(; J=J, h=h), MassGap(), Infinite())
end
