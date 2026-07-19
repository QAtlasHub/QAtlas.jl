# ─────────────────────────────────────────────────────────────────────────────
# TodaLattice — 1-D Toda lattice (classical, integrable; Toda 1967).
#
# Hamiltonian (standard "Toda 1967" convention):
#
#     H = Σ_n [ p_n² / 2 + V(q_{n+1} - q_n) ],
#       V(r) = (a/b) e^{-b r} + a r - a/b      (a, b > 0),
#
# with `a = b = 1` reducing to the textbook
#
#     V(r) = e^{-r} + r - 1,   V(0) = 0,   V''(0) = 1.
#
# The model is completely integrable (Hénon 1974, Flaschka 1974) and
# admits a Lax-pair representation L̇ = [B, L] whose spectrum is
# conserved.  N-soliton solutions in closed form are due to Toda
# (1967); the quantum Toda chain admits a Bethe-ansatz spectrum
# (Gutzwiller 1981, Sklyanin 1985).
#
# This Phase-1 entry establishes the model type and exposes the
# **gapless acoustic phonon** scalar: small-amplitude oscillations
# around the trivial classical ground state q_n = const, p_n = 0
# linearise to ω(k) = 2 √(a b) |sin(k/2)|, so the bottom of the
# small-oscillation spectrum is `MassGap = 0`.  Soliton dispersions
# and the Gutzwiller quantum-Toda spectrum require new quantity types
# (position / momentum-parameter-dependent scalars) and are tracked
# as a follow-up phase.
#
# References:
#   - M. Toda, [Toda1967](@cite).
#   - H. Flaschka, [Flaschka1974](@cite) — Lax representation.
#   - M. Hénon, [Henon1974](@cite) — integrals of motion.
#   - M. C. Gutzwiller, [Gutzwiller1981](@cite) — quantum Toda.
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: see file-header description above
#   Observable:  per src/core/quantities.jl (matches the dispatch tag)
#   Reference:   docs/src/conventions.md (project-wide convention policy)
#   STATUS:      backfilled by PR (audit gate); per-field domain content
#                left to a follow-up - see issue tracker for the model-specific
#                Hamiltonian sign / observable normalisation.

"""
    TodaLattice(; a::Real = 1.0, b::Real = 1.0) <: AbstractQAtlasModel

Classical 1-D Toda lattice (Toda 1967) with standard nearest-neighbour
exponential potential `V(r) = (a/b) e^{-b r} + a r - a/b`.  The
defaults `a = b = 1` give the canonical textbook form
`V(r) = e^{-r} + r - 1`.

Quantities registered (Phase 1):

| Quantity                       | BC         | Method                       |
| ------------------------------ | ---------- | ---------------------------- |
| [`MassGap`](@ref)              | `Infinite` | linear-phonon (gapless)      |

Soliton energy-momentum relations and the Gutzwiller (1981) quantum-
Toda Bethe-ansatz spectrum are tracked as Phase 2: they need
momentum-parameter / quantum-number-indexed quantity types not yet
in QAtlas core.

# References

- M. Toda, *J. Phys. Soc. Jpn.* **22**, 431 (1967).
- H. Flaschka, *Phys. Rev. B* **9**, 1924 (1974).
- M. C. Gutzwiller, *Annals Phys.* **133**, 304 (1981).
"""
struct TodaLattice <: AbstractQAtlasModel
    a::Float64
    b::Float64
end
TodaLattice(; a::Real=1.0, b::Real=1.0) = TodaLattice(Float64(a), Float64(b))

# ═══════════════════════════════════════════════════════════════════════════════
# Linearised phonon mass gap
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::TodaLattice, ::MassGap, ::Infinite; kwargs...) -> Float64

Small-amplitude phonon gap of the classical Toda lattice around the
uniform-spacing ground state.  The harmonic expansion gives
`ω(k) = 2 √(a b) |sin(k/2)|`, an acoustic branch vanishing linearly
at `k = 0`.  The corresponding mass gap is identically zero for any
`(a, b) > 0`.

# References

- M. Toda, *J. Phys. Soc. Jpn.* **22**, 431 (1967).
- H. Flaschka, *Phys. Rev. B* **9**, 1924 (1974).
"""
function fetch(::TodaLattice, ::MassGap, ::Infinite; kwargs...)
    return 0.0
end
