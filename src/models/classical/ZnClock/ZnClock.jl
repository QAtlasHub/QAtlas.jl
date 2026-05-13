# ─────────────────────────────────────────────────────────────────────────────
# ZnClock — 2-D Z_n clock model (n-state generalisation of the Ising,
# 3-state Potts, and XY universality classes).
#
#     H = -J Σ_{<i,j>} cos(2π (σ_i − σ_j) / n),    σ_i ∈ {0, …, n − 1},  J > 0.
#
# Phase diagram (José-Kadanoff-Kirkpatrick-Nelson 1977;
# Elitzur-Pearson-Shigemitsu 1979):
#
#   - n = 2:       identical to the 2-D Ising model (c = 1/2).
#   - n = 3:       identical to the 2-D 3-state Potts model (c = 4/5).
#   - n = 4:       Ashkin-Teller line — continuous family of c = 1 CFTs
#                  parameterised by a marginal coupling.
#   - n ≥ 5:       an intermediate massless XY-like BKT phase appears between
#                  the low-T ordered and high-T disordered phases, bounded by
#                  two BKT transitions (José et al. 1977); the critical line
#                  carries a continuous family of c = 1 CFTs.
#   - n → ∞:       reduces to the classical XY model (BKT transition, c = 1).
#
# Phase-1 scope: closed-form `CentralCharge` for n = 2 and n = 3 via
# delegation to the already-implemented Virasoro minimal models
# `MinimalModel(4, 3)` (Ising) and `MinimalModel(6, 5)` (3-state Potts).
# n = 4 (Ashkin-Teller line) and n ≥ 5 (BKT line) require a coupling
# specification to pick a CFT out of the c = 1 family and are deferred to
# Phase 2; the corresponding branches throw a descriptive `DomainError`.
#
# References:
#   - J. V. José, L. P. Kadanoff, S. Kirkpatrick, D. R. Nelson,
#     Phys. Rev. B 16, 1217 (1977).
#   - S. Elitzur, R. B. Pearson, J. Shigemitsu, Phys. Rev. D 19, 3698 (1979).
# ─────────────────────────────────────────────────────────────────────────────

"""
    ZnClock(n::Integer) <: AbstractQAtlasModel
    ZnClock(; n::Integer=2)

2-D Z_n clock model with `n` clock states.

Phase-1 quantities:

| Quantity                | BC         | Method                                 |
| ----------------------- | ---------- | -------------------------------------- |
| [`CentralCharge`](@ref) | `Infinite` | delegated to `MinimalModel` (n = 2, 3) |

For n ≥ 4 the critical theory is a continuous family of c = 1 CFTs
(Ashkin-Teller for n = 4; intermediate BKT phase for n ≥ 5) whose
selection requires a coupling parameter; these branches throw a
`DomainError` and are deferred to Phase 2.

# References

- J. V. José, L. P. Kadanoff, S. Kirkpatrick, D. R. Nelson,
  *Phys. Rev. B* **16**, 1217 (1977).
- S. Elitzur, R. B. Pearson, J. Shigemitsu,
  *Phys. Rev. D* **19**, 3698 (1979).
"""
struct ZnClock <: AbstractQAtlasModel
    n::Int
    function ZnClock(n::Integer)
        n ≥ 2 || throw(DomainError(n, "ZnClock requires n ≥ 2; got n = $n."))
        return new(Int(n))
    end
end
ZnClock(; n::Integer=2) = ZnClock(n)

# ═══════════════════════════════════════════════════════════════════════════════
# Central charge via MinimalModel delegation (n = 2, 3 only)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::ZnClock, ::CentralCharge, ::Infinite; n::Integer=m.n, kwargs...)
        -> Rational{Int}

Central charge of the 2-D Z_n clock model, currently supported for
n = 2 and n = 3 via delegation:

  - n = 2: `MinimalModel(4, 3)` → c = 1/2 (Ising universality).
  - n = 3: `MinimalModel(6, 5)` → c = 4/5 (3-state Potts universality).

For n = 4 (Ashkin-Teller line) and n ≥ 5 (intermediate BKT phase) the
critical theory is a *continuous family* of c = 1 CFTs and selecting a
particular member requires a coupling; those branches throw a
`DomainError` and are deferred to Phase 2.

# References

- J. V. José, L. P. Kadanoff, S. Kirkpatrick, D. R. Nelson,
  *Phys. Rev. B* **16**, 1217 (1977).
- S. Elitzur, R. B. Pearson, J. Shigemitsu,
  *Phys. Rev. D* **19**, 3698 (1979).
"""
function fetch(m::ZnClock, ::CentralCharge, ::Infinite; n::Integer=m.n, kwargs...)
    n ≥ 2 || throw(DomainError(n, "ZnClock CentralCharge requires n ≥ 2; got n = $n."))
    if n == 2
        # Z_2 = Ising universality.
        return QAtlas.fetch(QAtlas.MinimalModel(4, 3), CentralCharge())
    elseif n == 3
        # Z_3 = 3-state Potts universality.
        return QAtlas.fetch(QAtlas.MinimalModel(6, 5), CentralCharge())
    elseif n == 4
        throw(
            DomainError(
                n,
                "ZnClock CentralCharge: n = 4 (Ashkin-Teller line, continuous family of " *
                "c = 1 CFTs) requires coupling specification (Kadanoff-Brown 1979). " *
                "Deferred to Phase 2.",
            ),
        )
    else  # n ≥ 5
        throw(
            DomainError(
                n,
                "ZnClock CentralCharge: n ≥ 5 (intermediate BKT phase, José-Kadanoff-" *
                "Kirkpatrick-Nelson 1977) has a critical line of c = 1 CFTs between " *
                "the low-T ordered and high-T disordered phases. Deferred to Phase 2.",
            ),
        )
    end
end
