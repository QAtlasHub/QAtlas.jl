# core/universality.jl — the universality-class registry DESIGN.
#
# This is the *machinery* for the universality namespace (the parametric
# dispatch tag + the quantity tags it answers), kept in `core/` alongside the
# other registry types. The per-class DATA (critical exponents, central charge)
# and the universal-BEHAVIOUR implementations (entanglement scaling, Casimir
# correction) are registrations, and live under `universalities/` — class
# identity in `universalities/<Class>/`, cross-class universal behaviour in
# `universalities/behaviour/`.

# ─────────────────────────────────────────────────────────────────────────────
# Universality{C} — parametric type for universality classes
#
# Each universality class is identified by a Symbol parameter C
# (e.g., :Ising, :XY, :Heisenberg, :Potts3, :Percolation, :KPZ).
# The spatial dimension d is passed as a keyword argument to `fetch`.
# ─────────────────────────────────────────────────────────────────────────────

# `Universality{C}`, `CriticalExponents` and `GrowthExponents` are AbstractQAtlas's
# — see the note on the import list in QAtlas.jl. They were declared here as well,
# byte-for-byte, which made `QAtlas.Universality !== AbstractQAtlas.Universality` and
# cut this atlas out of every relation keyed on the base type. What stays in this
# file is the machinery that is genuinely QAtlas's: the `_is_universality` predicate
# and the per-class central charges.

# Namespace predicate: is `T` a Universality node?  Defined here (not in
# links.jl) so `register!` can derive `status=:universal` by construction — the
# universality `*_registry.jl` files run long before links.jl is included.
_is_universality(::Type{T}) where {T} = T <: Universality

raw"""
    _universality_central_charge(::Universality{C}) -> Real

Return the central charge `c` of the 1+1D CFT associated with
universality class `C`.  Used internally by the universal-behaviour
`CasimirEnergyCorrection` / Calabrese–Cardy entanglement fetches.

This generic method is the *contract* (and the error for classes without a
1+1D CFT). The per-class values are registered as class-identity data in
`universalities/<Class>/` (e.g. `:Ising` → `Ising2D/`, `:Potts3`/`:Potts4` →
`Potts/`, `:XY`/`:Heisenberg` → `ONModel/`):

| `C`           | `c`     | Reference                                            |
|---------------|---------|------------------------------------------------------|
| `:Ising`      | `1//2`  | BPZ minimal model M(3,4); Cardy 1986                 |
| `:Potts3`     | `4//5`  | M(5,6) minimal model; Dotsenko–Fateev 1984           |
| `:Potts4`     | `1//1`  | Free-boson radius limit; Saleur 1987                 |
| `:XY`         | `1//1`  | Compact free boson (1+1D Luttinger liquid)           |
| `:Heisenberg` | `1//1`  | SU(2)_1 WZW (Affleck–Haldane); 1+1D AFM chain        |

Other classes raise `ErrorException`.  Critical percolation is c = 0
(logarithmic, non-unitary), KPZ is non-equilibrium, and mean-field has no
1+1D CFT — all deliberately unsupported.
"""
function _universality_central_charge(::Universality{C}) where {C}
    return error(
        "QAtlas Universality{:$C}: no 1+1D CFT central charge is " *
        "registered for this universality class.  CasimirEnergyCorrection " *
        "requires `c` from a unitary 1+1D CFT; classes such as :KPZ " *
        "(non-equilibrium), :Percolation (non-unitary, c = 0 logarithmic), " *
        "and :MeanField (no 1+1D CFT) are not supported.  If `C` should " *
        "have a 1+1D CFT entry, add a method to " *
        "`QAtlas._universality_central_charge` in the class's " *
        "`universalities/<Class>/` file and document the source.",
    )
end
