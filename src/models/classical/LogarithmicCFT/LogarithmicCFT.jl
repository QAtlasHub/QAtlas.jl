# models/classical/LogarithmicCFT/LogarithmicCFT.jl
#
# Phase 1 implementation of the c = 0 logarithmic CFT (polymer / percolation
# universality class). Only the central charge is exposed; indecomposable
# representation structure and logarithmic operator content are deferred to
# Phase 2 (see issue #235).

# CONVENTION
#   Hamiltonian: see file-header description above
#   Observable:  per src/core/quantities.jl (matches the dispatch tag)
#   Reference:   docs/src/conventions.md (project-wide convention policy)
#   STATUS:      backfilled by PR (audit gate); per-field domain content
#                left to a follow-up - see issue tracker for the model-specific
#                Hamiltonian sign / observable normalisation.

"""
    LogarithmicCFT <: AbstractQAtlasModel

The c = 0 logarithmic conformal field theory describing the universality
class of self-avoiding polymers (Saleur 1992), critical percolation
(Cardy 2001), and dilute / dense polymer phases.

The theory is *logarithmic* because the dilatation operator `L_0` admits
indecomposable (non-diagonalisable) Jordan blocks — pairs of primary
fields share scaling dimensions but mix under conformal transformations,
producing logarithms in correlation functions. The total central charge
identically vanishes, yet the theory is non-trivial
(Pearce-Rasmussen-Zuber 2006; Vasseur-Jacobsen-Saleur 2011).

Phase 1 exposes only `CentralCharge = 0`. Logarithmic-operator structure,
indecomposable representations, and specific β-coupling parametrisations
are deferred to Phase 2.
"""
struct LogarithmicCFT <: AbstractQAtlasModel end
"""
    fetch(::LogarithmicCFT, ::CentralCharge, ::Infinite; kwargs...) -> Rational{Int}

Central charge of the c = 0 logarithmic CFT (polymer / percolation
universality):

    c = 0.

Despite the vanishing central charge the theory is non-trivial — the
dilatation operator `L_0` admits indecomposable Jordan blocks, mixing
pairs of primary fields under conformal transformations and producing
logarithms in correlation functions (Pearce-Rasmussen-Zuber 2006;
Vasseur-Jacobsen-Saleur 2011).

Phase 1 exposes only the central charge. Indecomposable representations
and logarithmic operator content require new quantity types and are
tracked as Phase 2.

# References

- H. Saleur, *Nucl. Phys. B* **382**, 486 (1992).
- J. Cardy, *J. Phys. A* **34**, 1419 (2001).
- P. A. Pearce, J. Rasmussen, J.-B. Zuber, *J. Stat. Mech.* P11017 (2006).
- R. Vasseur, J. L. Jacobsen, H. Saleur, *J. Stat. Mech.* L07001 (2011).
"""
function fetch(::LogarithmicCFT, ::CentralCharge, ::Infinite; kwargs...)
    return 0 // 1
end
