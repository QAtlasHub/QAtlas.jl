# ─────────────────────────────────────────────────────────────────────────────
# bounds/Holographic — holographic / gravitational entropy bounds.
#
# Extracted from the former `Universality(:QuantumMechanics)` dumping ground:
# these are model-independent *bounds*, not a universality class.
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Bound{:Holographic}, ::BekensteinBound, ::Infinite; R, E)

Bekenstein 1981 universal upper bound on the entropy of a system of
radius `R` and total energy `E`,

    S ≤ 2π R E          (ℏ = c = k_B = 1).

Saturated (up to an O(1) factor) by black holes.  A `status=:bound`,
`direction=:upper` claim.
"""
function fetch(
    ::Bound{:Holographic}, ::BekensteinBound, ::Infinite; R::Real, E::Real, kwargs...
)
    (R ≥ 0 && E ≥ 0) ||
        throw(ArgumentError("BekensteinBound: R, E must be ≥ 0; got R=$(R), E=$(E)"))
    return 2π * R * E
end
