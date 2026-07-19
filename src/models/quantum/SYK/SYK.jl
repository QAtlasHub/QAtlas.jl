# ─────────────────────────────────────────────────────────────────────────────
# SYK — Sachdev-Ye-Kitaev model (large-N maximally chaotic
# random-fermion quantum mechanics).
#
# Hamiltonian (q-body Majorana coupling, Maldacena-Stanford 2016
# convention):
#
#     H = (i^{q/2} / q!) Σ_{i_1<…<i_q} J_{i_1…i_q} ψ_{i_1} … ψ_{i_q},
#         ⟨J²⟩ = J² (q-1)! / N^{q-1},
#
# with `N` Majorana fermions `ψ_i` and all-to-all `q`-body Gaussian
# disorder.  In the large-`N`, low-energy (IR) limit the disorder-
# averaged theory becomes **emergently conformal**: the Schwarzian
# reparametrisation mode emerges as the soft sector and the Majorana
# fermion acquires a universal IR conformal dimension
#
#     Δ_ψ = 1 / q                  (Kitaev 2015; Maldacena-Stanford 2016)
#
# i.e. `⟨ψ(τ_1) ψ(τ_2)⟩ ∝ 1 / |τ_12|^{2 Δ_ψ}` at large `|τ_12|`.
# SYK is the prototype of holographic 1-D quantum mechanics dual to
# AdS₂ dilaton gravity / JT gravity, and saturates the Maldacena-
# Shenker-Stanford bound on chaos `λ_L = 2π / β`.
#
# Phase 1 of this file therefore exposes only the IR Majorana
# conformal dimension `Δ_ψ = 1 / q` via `ConformalWeights` at
# `field = :ψ`.  Composite-operator dimensions (e.g. the bilinear
# `O_n = ψ ∂_τ^{2n+1} ψ` tower), the Schwarzian zero-temperature
# entropy, the four-point function reproducing JT bulk
# diagrammatics, and the maximal Lyapunov exponent are deferred
# to Phase 2.
#
# References:
#   - S. Sachdev, J. Ye, [SachdevYe1993](@cite).
#   - A. Kitaev, "A simple model of quantum holography",
#     KITP talks (2015).
#   - J. Maldacena, D. Stanford, [MaldacenaStanford2016](@cite).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Fermion bilinears c†c
#   Observable:  Fermion (number n = c†c, bilinear ⟨c†_i c_j⟩); derived spin observables follow spin S = σ/2
#   Reference:   docs/src/conventions.md §Fermion convention

"""
    SYK(; q::Integer = 4) <: AbstractQAtlasModel

Sachdev-Ye-Kitaev model: `N` Majorana fermions with all-to-all
random `q`-body coupling.  The defining structural parameter for
the universal large-`N` IR is `q` (which must be even — Majorana
antisymmetry — and `≥ 2`).  The fermion number `N` and the
coupling scale `J` set the UV / finite-`N` data but do **not**
enter the leading large-`N` IR conformal dimension exposed here.

Default `q = 4` is the most-studied case: 4-body Majorana
coupling, exactly soluble Schwarzian-dressed two-point function,
and a maximal Lyapunov exponent saturating the Maldacena-
Shenker-Stanford bound.

Quantities registered (Phase 1):

| Quantity                        | BC         | Method                              |
| ------------------------------- | ---------- | ----------------------------------- |
| [`ConformalWeights`](@ref)      | `Infinite` | analytic large-N IR (`Δ_ψ = 1/q`)   |

Phase 2 will add composite-operator dimensions (bilinear tower),
the Schwarzian zero-T entropy, maximal Lyapunov saturation, and
the JT-gravity bulk dual.

# References

- S. Sachdev, J. Ye, *Phys. Rev. Lett.* **70**, 3339 (1993).
- A. Kitaev, KITP talks (2015) — "A simple model of quantum holography".
- J. Maldacena, D. Stanford, *Phys. Rev. D* **94**, 106002 (2016).
"""
struct SYK <: AbstractQAtlasModel
    q::Int
    function SYK(q::Integer)
        q ≥ 2 ||
            throw(DomainError(q, "SYK requires q-body coupling with q ≥ 2; got q = $q."))
        iseven(q) ||
            throw(DomainError(q, "SYK requires q even (Majorana indices); got q = $q."))
        return new(Int(q))
    end
end
SYK(; q::Integer=4) = SYK(q)

# ═══════════════════════════════════════════════════════════════════════════════
# IR Majorana conformal dimension — Δ_ψ = 1 / q
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::SYK, ::ConformalWeights, ::Infinite;
          q::Integer = m.q, field::Symbol = :ψ, kwargs...) -> Rational{Int}

Large-`N` IR (low-energy / emergent-conformal) Majorana fermion
conformal dimension of the Sachdev-Ye-Kitaev model:

    Δ_ψ = 1 / q.

This is the universal scaling of the disorder-averaged two-point
function

    ⟨ψ(τ_1) ψ(τ_2)⟩ ∝ 1 / |τ_12|^{2 Δ_ψ}

in the IR conformal window where the Schwarzian reparametrisation
soft mode dominates the dynamics (Kitaev 2015; Maldacena-Stanford
2016).  Finite-`N` and subleading-`1/(βJ)` Schwarzian corrections,
as well as composite-operator dimensions (the `O_n = ψ ∂_τ^{2n+1} ψ`
bilinear tower), are deferred to Phase 2.

# Arguments

- `q::Integer = m.q`: body-count of the Majorana coupling.  Must be
  even and `≥ 2`.
- `field::Symbol = :ψ`: which IR operator to query.  Phase 1
  supports only `:ψ` (the elementary Majorana fermion); any other
  symbol raises `DomainError`, deferring composite-operator
  dimensions to Phase 2.

# References

- S. Sachdev, J. Ye, *Phys. Rev. Lett.* **70**, 3339 (1993).
- A. Kitaev, KITP talks (2015).
- J. Maldacena, D. Stanford, *Phys. Rev. D* **94**, 106002 (2016).
"""
function fetch(
    m::SYK, ::ConformalWeights, ::Infinite; q::Integer=m.q, field::Symbol=:ψ, kwargs...
)
    q ≥ 2 || throw(DomainError(q, "SYK ConformalWeights requires q ≥ 2; got q = $q."))
    iseven(q) || throw(
        DomainError(
            q,
            "SYK ConformalWeights requires q even (q-body Majorana coupling); got q = $q.",
        ),
    )
    if field == :ψ
        return 1 // q
    else
        throw(
            DomainError(
                field,
                "SYK ConformalWeights: Phase 1 supports only field=:ψ (Majorana fermion IR dimension); composite-operator dimensions Phase 2. Got field=:$field.",
            ),
        )
    end
end
