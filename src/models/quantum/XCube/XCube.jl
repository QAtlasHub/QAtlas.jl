# ─────────────────────────────────────────────────────────────────────────────
# XCube — fracton X-cube model (Vijay-Haah-Fu 2016).
#
# Hamiltonian (qubits on cube edges, both vertex and cube terms):
#
#     H = -Σ_v A_v^x - Σ_v A_v^y - Σ_v A_v^z - Σ_c B_c,
#
# where `A_v^μ` is the product of 4 σ^z spins around the cross at
# vertex `v` perpendicular to direction `μ`, and `B_c` is the
# product of 12 σ^x spins on the edges of cube `c`.  The model
# realises Type-I fracton topological order: anyon-like excitations
# are confined to live on 1-D / 2-D subsystems of the cubic lattice.
#
# On a closed periodic L_x × L_y × L_z cubic torus the ground-state
# degeneracy is (Vijay-Haah-Fu 2016; Slagle-Kim 2017)
#
#     log_2 GSD = 2 (L_x + L_y + L_z) − 3,
#
# i.e.
#
#     GSD = 2^{2 (L_x + L_y + L_z) − 3}.
#
# The subextensive (linear) scaling in `L` is the hallmark feature
# distinguishing fracton order from conventional topological order
# (e.g. ToricCode `GSD = 4^g`, independent of lattice size).
#
# This Phase-1 entry registers the closed-form `GroundStateDegeneracy`
# at PBC.  Subsystem-symmetry-protected excitation mobility (planar
# fractons, lineons), braiding / fusion rules, and the closely-related
# Haah's code (fully-mobile-free) are tracked as Phase 2.
#
# References:
#   - J. Haah, Phys. Rev. A 83, 042330 (2011).
#   - S. Vijay, J. Haah, L. Fu, Phys. Rev. B 94, 235157 (2016).
#   - K. Slagle, Y. B. Kim, Phys. Rev. B 96, 195139 (2017).
# ─────────────────────────────────────────────────────────────────────────────

"""
    XCube() <: AbstractQAtlasModel

Fracton X-cube model on the 3-D cubic lattice (Vijay-Haah-Fu 2016) —
prototype of Type-I fracton topological order with subextensive
ground-state degeneracy that scales linearly with each lattice
direction.

Quantities registered (Phase 1):

| Quantity                         | BC    | Method                                          |
| -------------------------------- | ----- | ----------------------------------------------- |
| [`GroundStateDegeneracy`](@ref)  | `PBC` | analytic (Vijay-Haah-Fu 2016 / Slagle-Kim 2017) |

# References

- J. Haah, *Phys. Rev. A* **83**, 042330 (2011).
- S. Vijay, J. Haah, L. Fu, *Phys. Rev. B* **94**, 235157 (2016).
- K. Slagle, Y. B. Kim, *Phys. Rev. B* **96**, 195139 (2017).
"""
struct XCube <: AbstractQAtlasModel end

# ═══════════════════════════════════════════════════════════════════════════════
# Ground-state degeneracy on a closed cubic torus
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::XCube, ::GroundStateDegeneracy, ::PBC;
          Lx::Int, Ly::Int, Lz::Int) -> BigInt

Subextensive ground-state degeneracy of the X-cube fracton model on
a closed L_x × L_y × L_z cubic torus (Vijay-Haah-Fu 2016 /
Slagle-Kim 2017):

    GSD(L_x, L_y, L_z) = 2^{2 (L_x + L_y + L_z) − 3}.

Returned as a `BigInt` because the exponent grows linearly in the
side lengths and can exceed `Int64` for modest L.  All `L_α ≥ 2` is
required (smaller tori fail the well-defined-stabiliser conditions).

# References

- S. Vijay, J. Haah, L. Fu, *Phys. Rev. B* **94**, 235157 (2016).
- K. Slagle, Y. B. Kim, *Phys. Rev. B* **96**, 195139 (2017).
"""
function fetch(
    ::XCube, ::GroundStateDegeneracy, ::PBC; Lx::Int, Ly::Int, Lz::Int, kwargs...
)
    (Lx ≥ 2 && Ly ≥ 2 && Lz ≥ 2) || throw(
        DomainError(
            (Lx, Ly, Lz),
            "XCube GSD: each L_α ≥ 2 is required for well-defined X-cube stabilisers; got (Lx, Ly, Lz) = ($Lx, $Ly, $Lz).",
        ),
    )
    exponent = 2 * (Lx + Ly + Lz) - 3
    return big(2)^exponent
end
