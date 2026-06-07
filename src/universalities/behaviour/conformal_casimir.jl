# universalities/behaviour/conformal_casimir.jl
#
# Universal BEHAVIOUR: the CFT Casimir / finite-size ground-state energy
# correction (Cardy 1986). Cross-class — it depends only on the central charge
# `c` of the class's 1+1D CFT (registered as class identity under
# `universalities/<Class>/`), not on the concrete model.
#
# For a 1+1D CFT with central charge c and CFT velocity v on size L:
#
#   E_0^PBC(L) = L · ε_∞          - π c v / (6  L) + O(L^{-2})
#   E_0^OBC(L) = L · ε_∞ + ε_surf - π c v / (24 L) + O(L^{-2})
#
# Only the universal 1/L correction is exposed. The PBC : OBC ratio of 4 is a
# kinematic consequence of the conformal cylinder ↔ strip map.
#
# References:
#   J. Cardy, Nucl. Phys. B 270, 186 (1986).
#   H. W. J. Blöte, J. L. Cardy, M. P. Nightingale, Phys. Rev. Lett. 56, 742 (1986).
#   I. Affleck, Phys. Rev. Lett. 56, 746 (1986).

raw"""
    fetch(::Universality{C}, ::CasimirEnergyCorrection, ::PBC; L, v) -> Real

Return the universal Cardy 1/L correction
``-\pi c v / (6 L)`` at periodic boundary conditions, where `c` is
the central charge of the 1+1D CFT for class `C` (see
[`_universality_central_charge`](@ref)) and `v` is the CFT velocity
supplied by the caller.

`L` and `v` must be positive.  The return type is `Rational` when both
`c` is rational and `v` is rational/integer; otherwise `Float64`.
"""
function fetch(
    u::Universality{C}, ::CasimirEnergyCorrection, ::PBC; L::Real, v::Real, kwargs...
) where {C}
    L > 0 || throw(ArgumentError("CasimirEnergyCorrection: L must be positive; got $L"))
    v > 0 || throw(ArgumentError("CasimirEnergyCorrection: v must be positive; got $v"))
    c = _universality_central_charge(u)
    return -π * c * v / (6 * L)
end

raw"""
    fetch(::Universality{C}, ::CasimirEnergyCorrection, ::OBC; L, v) -> Real

Return the universal Cardy 1/L correction
``-\pi c v / (24 L)`` at open boundary conditions.

The PBC : OBC ratio of the 1/L term is exactly 4, a kinematic
consequence of the conformal map between the cylinder (PBC) and the
strip (OBC).
"""
function fetch(
    u::Universality{C}, ::CasimirEnergyCorrection, ::OBC; L::Real, v::Real, kwargs...
) where {C}
    L > 0 || throw(ArgumentError("CasimirEnergyCorrection: L must be positive; got $L"))
    v > 0 || throw(ArgumentError("CasimirEnergyCorrection: v must be positive; got $v"))
    c = _universality_central_charge(u)
    return -π * c * v / (24 * L)
end
