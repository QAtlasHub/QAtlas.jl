# ─────────────────────────────────────────────────────────────────────────────
# SpinIce — classical spin ice on the pyrochlore lattice.
#
# Pauling 1935 originally derived the celebrated residual entropy of
# common water ice from the proton-disorder ice rule on a tetrahedral
# H₂O network.  The same combinatorics governs spin ice — Ising-like
# moments (Dy₂Ti₂O₇, Ho₂Ti₂O₇, …) on the pyrochlore lattice of
# corner-sharing tetrahedra — under the "2-in-2-out" ground-state
# constraint (Bramwell-Gingras 2001).
#
# Pauling derivation (per spin):
#   Each tetrahedron admits C(4,2) = 6 ice-rule configurations out of
#   the 16 a-priori Ising states; with two tetrahedra per spin in the
#   thermodynamic limit
#       W ≈ 2^N · (6/16)^{N/2}
#   so
#       S/N = log 2 + (1/2) log(3/8) = (1/2) log(3/2) ≈ 0.20273.
#
# More accurate values (Nagle 1966 series; Berg-Bohm-Kotrla numerical
# transfer-matrix) sit a few percent above the Pauling estimate; this
# file ships the textbook closed form.
#
# Coulomb-phase dipolar correlations and emergent magnetic monopoles
# (Castelnovo-Moessner-Sondhi 2008) require correlator and excitation
# infrastructure beyond bare residual entropy; they are tracked as
# follow-up scope (#257 phase 2).
#
# References:
#   - L. Pauling, [Pauling1935](@cite).
#   - S. T. Bramwell, M. J. P. Gingras, [BramwellGingras2001](@cite).
#   - C. Castelnovo, R. Moessner, S. L. Sondhi, [CastelnovoMoessnerSondhi2008](@cite).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: see file-header description above
#   Observable:  per src/core/quantities.jl (matches the dispatch tag)
#   Reference:   docs/src/conventions.md (project-wide convention policy)
#   STATUS:      backfilled by PR (audit gate); per-field domain content
#                left to a follow-up - see issue tracker for the model-specific
#                Hamiltonian sign / observable normalisation.

"""
    SpinIce() <: AbstractQAtlasModel

Classical spin ice on the pyrochlore lattice of corner-sharing
tetrahedra.  Ising-like local moments obey the "2-in-2-out" ice rule
in the ground-state manifold; the residual entropy follows the
Pauling 1935 closed form.

Quantities registered:

| Quantity                       | BC         | Method            |
| ------------------------------ | ---------- | ----------------- |
| [`ResidualEntropy`](@ref)      | `Infinite` | analytic (Pauling)|

The Coulomb-phase dipolar spin correlations and the
Castelnovo-Moessner-Sondhi (2008) emergent-monopole picture are
tracked as a follow-up phase (#257 phase 2) and not exposed here.

# References

- L. Pauling, *J. Am. Chem. Soc.* **57**, 2680 (1935).
- S. T. Bramwell, M. J. P. Gingras, *Science* **294**, 1495 (2001).
"""
struct SpinIce <: AbstractQAtlasModel end
# ═══════════════════════════════════════════════════════════════════════════════
# Pauling residual entropy
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::SpinIce, ::ResidualEntropy, ::Infinite; kwargs...) -> Float64

Zero-temperature residual entropy per spin of pyrochlore spin ice in
the Pauling 1935 closed form

    S/N = (1/2) log(3/2) ≈ 0.20273.

The estimate ignores ice-rule correlations between neighbouring
tetrahedra; subsequent series (Nagle 1966) and numerical transfer
matrix studies put the true value a few percent higher.  The textbook
closed form is shipped here because it is exact under the
mean-tetrahedron approximation and remains the canonical
constant-of-Nature for spin-ice phenomenology.

# References

- L. Pauling, *J. Am. Chem. Soc.* **57**, 2680 (1935).
"""
function fetch(::SpinIce, ::ResidualEntropy, ::Infinite; kwargs...)
    return 0.5 * log(3 / 2)
end
