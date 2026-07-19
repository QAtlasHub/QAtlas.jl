# ─────────────────────────────────────────────────────────────────────────────
# Compass1D — 1D alternating-bond compass chain.
#
# The 1D limit of the 2D Kugel–Khomskii / orbital compass model:
#
#   H = -J_x Σ_{i odd}  σ^x_i σ^x_{i+1}
#       -J_y Σ_{i even} σ^y_i σ^y_{i+1},        J_x, J_y > 0.
#
# Via Jordan–Wigner this maps to a dimerised Kitaev / p-wave wire with
# alternating couplings; the ground-state problem is exactly solvable.
# The bulk mass gap is
#
#       Δ = 2 |J_x − J_y|,
#
# which collapses to Δ = 0 at the symmetric point J_x = J_y — a
# first-order quantum phase transition between the X-bond-dominated and
# Y-bond-dominated dimerised ground states (Brzezicki–Dziarmaga–Oleś
# 2007).
#
# Phase 1 (this file) exposes a single closed-form fetch:
#
#   • `MassGap`, `Infinite`  →  2 · |J_x − J_y|.
#
# References:
#   • W. Brzezicki, J. Dziarmaga, A. M. Oleś, "Quantum phase transition
#     in the one-dimensional compass model", PRB 75, 134415 (2007).
#   • K. I. Kugel, D. I. Khomskii, "The Jahn–Teller effect and
#     magnetism: transition metal compounds", Sov. Phys. Usp. 25, 231
#     (1982) — original orbital compass / Kugel–Khomskii model.
#   • A. Kitaev, "Anyons in an exactly solved model and beyond", Annals
#     of [Kitaev2006](@cite) — JW-dual dimerised Kitaev chain.
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Pauli σ (this file)
#   Observable:  Spin S = σ/2  (QAtlas-wide spin convention; see docs/src/conventions.md)

"""
    Compass1D(; J_x::Real = 1.0, J_y::Real = 1.0) <: AbstractQAtlasModel

1D alternating-bond compass chain

    H = -J_x Σ_{i odd}  σ^x_i σ^x_{i+1}
        -J_y Σ_{i even} σ^y_i σ^y_{i+1},

with `J_x, J_y > 0`.  Dual to a dimerised Kitaev / p-wave wire via
Jordan–Wigner; exactly solvable.

# Phase 1 scope (this release)

Phase 1 exposes the closed-form bulk gap only:

- [`MassGap`](@ref) at `Infinite()` → `Δ = 2 |J_x − J_y|`.

The gap closes at the symmetric point `J_x = J_y`, which is a
**first-order** quantum phase transition between the X-bond-dimerised
and Y-bond-dimerised ground states (Brzezicki–Dziarmaga–Oleś 2007).

# References

- Brzezicki–Dziarmaga–Oleś, *PRB* **75**, 134415 (2007).
- Kugel–Khomskii, *Sov. Phys. Usp.* **25**, 231 (1982).
- Kitaev, *Annals of Physics* **321**, 2 (2006).
"""
struct Compass1D <: AbstractQAtlasModel
    J_x::Float64
    J_y::Float64
    function Compass1D(J_x::Real, J_y::Real)
        J_x > 0 || throw(DomainError(J_x, "Compass1D requires J_x > 0; got J_x = $J_x."))
        J_y > 0 || throw(DomainError(J_y, "Compass1D requires J_y > 0; got J_y = $J_y."))
        return new(Float64(J_x), Float64(J_y))
    end
end
Compass1D(; J_x::Real=1.0, J_y::Real=1.0) = Compass1D(J_x, J_y)

# ─── fetch methods ────────────────────────────────────────────────────

"""
    fetch(model::Compass1D, ::MassGap, ::Infinite) -> Float64

Bulk mass gap of the 1D alternating-bond compass chain:

    Δ = 2 |J_x − J_y|.

Closed-form result obtained from the Jordan–Wigner dual dimerised
Kitaev / p-wave wire (Brzezicki–Dziarmaga–Oleś 2007).  Δ vanishes at
the symmetric point `J_x = J_y`, which is a first-order quantum phase
transition between the X-bond-dimerised and Y-bond-dimerised ground
states.

# References

- Brzezicki–Dziarmaga–Oleś, *PRB* **75**, 134415 (2007).
"""
function fetch(model::Compass1D, ::MassGap, ::Infinite; kwargs...)
    return 2.0 * abs(model.J_x - model.J_y)
end
