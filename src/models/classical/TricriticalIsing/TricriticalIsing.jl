# ─────────────────────────────────────────────────────────────────────────────
# TricriticalIsing — tricritical Ising CFT M(5, 4).
#
# The tricritical Ising model is the next-to-Ising member of the
# unitary minimal series M(p+1, p); it sits at the tricritical point
# of the Blume-Capel (vacancy-extended Ising) model and is also
# realised by the (p+1, p) = (5, 4) RSOS lattice models of
# Andrews-Baxter-Forrester (1984).
#
# Central charge:
#
#     c(5, 4) = 1 − 6 (5 − 4)² / (5 · 4) = 1 − 6/20 = 1 − 3/10 = 7/10.
#
# Independent Kac primaries (modulo Kac symmetry, total = 6),
# labelled by their literature (r̃, s̃) Kac indices.  Note QAtlas's
# MinimalModel API uses (r, s) = (s̃, r̃) so that r ∈ 1:(p_prime−1)
# and s ∈ 1:(p−1):
#
#     1   (identity)            h = 0       Kac (1,1)
#     ε   (energy)              h = 1/10    Kac (2,1)   API (1,2)
#     ε'  (vacancy)             h = 3/5     Kac (2,3)   API (3,2)
#     ε'' (irrelevant)          h = 3/2     Kac (1,3)   API (3,1)
#     σ   (spin)                h = 3/80    Kac (2,2)   API (2,2)
#     σ'  (subleading spin)     h = 7/16    Kac (1,2)   API (2,1)
#
# This Phase-1 entry delegates all three registered quantities
# (`CentralCharge`, `ConformalWeights`, `PrimaryFields`) to the
# already-implemented [`MinimalModel(5, 4)`](@ref).
#
# References:
#   - A. A. Belavin, A. M. Polyakov, A. B. Zamolodchikov,
#     Nucl. Phys. B 241, 333 (1984).
#   - D. Friedan, Z. Qiu, S. Shenker, Phys. Rev. Lett. 52, 1575 (1984).
#   - G. E. Andrews, R. J. Baxter, P. J. Forrester,
#     J. Stat. Phys. 35, 193 (1984).
# ─────────────────────────────────────────────────────────────────────────────

"""
    TricriticalIsing() <: AbstractQAtlasModel

Tricritical Ising CFT — the unitary Virasoro minimal model M(5, 4),
next-to-simplest unitary 2D CFT after Ising (M(4, 3)).  Lattice
realisations include the Blume-Capel model at its tricritical point
and the Andrews-Baxter-Forrester (1984) RSOS models at (p+1, p) = (5, 4).

Quantities registered (Phase 1):

| Quantity                       | BC         | Method                         |
| ------------------------------ | ---------- | ------------------------------ |
| [`CentralCharge`](@ref)        | `Infinite` | delegated to MinimalModel(5,4) |
| [`ConformalWeights`](@ref)     | `Infinite` | delegated to MinimalModel(5,4) |
| [`PrimaryFields`](@ref)        | `Infinite` | delegated to MinimalModel(5,4) |

# References

- A. A. Belavin, A. M. Polyakov, A. B. Zamolodchikov,
  *Nucl. Phys. B* **241**, 333 (1984).
- D. Friedan, Z. Qiu, S. Shenker, *Phys. Rev. Lett.* **52**, 1575 (1984).
"""
# CONVENTION
#   Hamiltonian: see file-header description above
#   Observable:  per src/core/quantities.jl (matches the dispatch tag)
#   Reference:   docs/src/conventions.md (project-wide convention policy)
#   STATUS:      backfilled by PR (audit gate); per-field domain content
#                left to a follow-up - see issue tracker for the model-specific
#                Hamiltonian sign / observable normalisation.

struct TricriticalIsing <: AbstractQAtlasModel end

# ═══════════════════════════════════════════════════════════════════════════════
# Central charge via MinimalModel(5, 4) delegation
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::TricriticalIsing, ::CentralCharge, ::Infinite; kwargs...)
        -> Rational{Int}

Central charge of the tricritical Ising CFT, delegated to
[`MinimalModel(5, 4)`](@ref):

    c = 7/10.

Returned as an exact `Rational{Int}` (`7//10`).
"""
function fetch(::TricriticalIsing, ::CentralCharge, ::Infinite; kwargs...)
    return QAtlas.fetch(QAtlas.MinimalModel(5, 4), CentralCharge())
end

# ═══════════════════════════════════════════════════════════════════════════════
# Conformal weights via MinimalModel(5, 4) Kac formula
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::TricriticalIsing, ::ConformalWeights, ::Infinite; r, s, kwargs...)
        -> Rational{Int}

Kac conformal weight h_{r,s} of the tricritical Ising CFT, delegated
to [`MinimalModel(5, 4)`](@ref).  Index range: 1 ≤ r ≤ p_prime − 1 = 3,
1 ≤ s ≤ p − 1 = 4, with

    h_{r,s} = ((p r − p_prime s)² − (p − p_prime)²) / (4 p p_prime)
            = ((5 r − 4 s)² − 1) / 80.

Famous primaries (literature (r̃, s̃) Kac labels mapped to the API
(r, s) = (s̃, r̃) convention):

| field         | lit. (r̃, s̃) | API (r, s) | h     |
| ------------- | ------------- | ---------- | ----- |
| 1 (identity)  | (1,1)         | (1,1)      | 0     |
| σ (spin)      | (2,2)         | (2,2)      | 3/80  |
| σ' (subld.)   | (1,2)         | (2,1)      | 7/16  |
| ε (energy)    | (2,1)         | (1,2)      | 1/10  |
| ε' (vacancy)  | (2,3)         | (3,2)      | 3/5   |
| ε'' (irrel.)  | (1,3)         | (3,1)      | 3/2   |
"""
function fetch(
    ::TricriticalIsing,
    ::ConformalWeights,
    ::Infinite;
    r::Integer=1,
    s::Integer=1,
    kwargs...,
)
    return QAtlas.fetch(QAtlas.MinimalModel(5, 4), ConformalWeights(); r=r, s=s)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Primary fields via MinimalModel(5, 4) delegation
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::TricriticalIsing, ::PrimaryFields, ::Infinite; kwargs...)
        -> Vector{NamedTuple}

Independent Kac primaries of M(5, 4) modulo the Kac symmetry
(r, s) ↔ (p_prime − r, p − s) = (4 − r, 5 − s) in the MinimalModel
API convention.  Total: 6 primaries.
"""
function fetch(::TricriticalIsing, ::PrimaryFields, ::Infinite; kwargs...)
    return QAtlas.fetch(QAtlas.MinimalModel(5, 4), PrimaryFields())
end
