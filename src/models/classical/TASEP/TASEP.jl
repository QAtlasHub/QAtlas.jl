# ─────────────────────────────────────────────────────────────────────────────
# TASEP — totally asymmetric simple exclusion process.
#
# Prototypical 1-D non-equilibrium lattice gas (Derrida-Evans-Hakim-Pasquier
# 1993; Derrida-Lebowitz 1998) in which particles on a 1-D lattice with PBC
# hop rightward at rate `p` if the next site is empty.  TASEP is the
# canonical KPZ-universality-class model for driven diffusive systems and
# the steady-state current
#
#     j(ρ) = p ρ (1 − ρ)
#
# is the prototypical non-equilibrium observable.  It is maximised at
# half-filling `ρ = 1/2` with `j_max = p/4`, and is exactly symmetric
# under the particle-hole exchange `ρ ↔ 1 − ρ`.
#
# Phase-1 entry registers only `SteadyStateCurrent`.  The full Derrida-
# Lebowitz cumulant generating function and finite-size corrections live
# in later phases (issue #241).
#
# References:
#   - M. Kardar, G. Parisi, Y.-C. Zhang,
#     Phys. Rev. Lett. 56, 889 (1986) — KPZ universality class.
#   - B. Derrida, M. R. Evans, V. Hakim, V. Pasquier,
#     J. Phys. A 26, 1493 (1993).
#   - B. Derrida, J. L. Lebowitz, Phys. Rev. Lett. 80, 209 (1998).
# ─────────────────────────────────────────────────────────────────────────────

"""
    TASEP(; p::Real = 1.0, ρ::Real = 0.5) <: AbstractQAtlasModel

Totally Asymmetric Simple Exclusion Process — 1-D KPZ-class non-equilibrium
lattice gas with rightward hopping rate `p > 0` at particle density
`0 ≤ ρ ≤ 1`.  The default `(p, ρ) = (1, 1/2)` is the maximal-current
half-filling point with `j_max = 1/4`.

Quantities registered (Phase 1):

| Quantity                           | BC         | Method                          |
| ---------------------------------- | ---------- | ------------------------------- |
| [`SteadyStateCurrent`](@ref)       | `Infinite` | analytic (`j = p ρ (1 − ρ)`)    |

# References

- M. Kardar, G. Parisi, Y.-C. Zhang,
  *Phys. Rev. Lett.* **56**, 889 (1986) — KPZ universality class.
- B. Derrida, M. R. Evans, V. Hakim, V. Pasquier,
  *J. Phys. A* **26**, 1493 (1993).
- B. Derrida, J. L. Lebowitz, *Phys. Rev. Lett.* **80**, 209 (1998).
"""
# CONVENTION
#   Hamiltonian: see file-header description above
#   Observable:  per src/core/quantities.jl (matches the dispatch tag)
#   Reference:   docs/src/conventions.md (project-wide convention policy)
#   STATUS:      backfilled by PR (audit gate); per-field domain content
#                left to a follow-up - see issue tracker for the model-specific
#                Hamiltonian sign / observable normalisation.

struct TASEP <: AbstractQAtlasModel
    p::Float64        # hopping rate
    ρ::Float64        # particle density
    function TASEP(p::Real, ρ::Real)
        p > 0 || throw(DomainError(p, "TASEP requires hopping rate p > 0; got p = $p."))
        (0 ≤ ρ ≤ 1) ||
            throw(DomainError(ρ, "TASEP requires density 0 ≤ ρ ≤ 1; got ρ = $ρ."))
        return new(Float64(p), Float64(ρ))
    end
end
TASEP(; p::Real=1.0, ρ::Real=0.5) = TASEP(p, ρ)

# ═══════════════════════════════════════════════════════════════════════════════
# Steady-state current
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::TASEP, ::SteadyStateCurrent, ::Infinite; p=m.p, ρ=m.ρ) -> Float64

Mean-field steady-state particle current of TASEP at hopping rate `p > 0`
and density `0 ≤ ρ ≤ 1` (Derrida-Lebowitz 1998):

    j(ρ) = p ρ (1 − ρ).

Special points:

- `ρ = 1/2`  →  `j_max = p/4`      (maximal current, half-filling)
- `ρ = 0` or `ρ = 1`  →  `j = 0`   (empty / fully packed lattice)
- `ρ ↔ 1 − ρ`                       (particle-hole symmetry)

`p ≤ 0` or `ρ ∉ [0, 1]` raises `DomainError`.

# References

- B. Derrida, J. L. Lebowitz, *Phys. Rev. Lett.* **80**, 209 (1998).
"""
function fetch(
    m::TASEP, ::SteadyStateCurrent, ::Infinite; p::Real=m.p, ρ::Real=m.ρ, kwargs...
)
    p > 0 || throw(DomainError(p, "TASEP SteadyStateCurrent requires p > 0; got p = $p."))
    (0 ≤ ρ ≤ 1) ||
        throw(DomainError(ρ, "TASEP SteadyStateCurrent requires 0 ≤ ρ ≤ 1; got ρ = $ρ."))
    return p * ρ * (1 - ρ)
end
