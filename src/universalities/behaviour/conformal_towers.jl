# universalities/behaviour/conformal_towers.jl
#
# Universal BEHAVIOUR: Conformal tower of states excitation spectrum in 1+1D CFT.
#
# References:
#   J. Cardy, [Cardy1986](@cite).
#   H. W. J. Blöte, J. L. Cardy, M. P. Nightingale, [BloteCardyNightingale1986](@cite).

raw"""
    fetch(::Universality{C}, ::ConformalTower, bc::Union{PBC, OBC}; L::Real, v::Real) -> Vector{NamedTuple}

Return a sorted vector of NamedTuples representing the lowest-lying excitation energies
`E_n - E_0` relative to the ground state, their scaling dimensions `Δ_n` (or `h_n`),
and their degeneracies in the conformal tower of states:

    (energy = E_n - E_0, dimension = Δ_n, degeneracy = g_n)

For PBC (torus geometry):
    E_n - E_0 = (2π v / L) Δ_n

where `Δ_n = h_n + h̄_n` is the total scaling dimension. Degeneracy counts all states
at that level including all spin sectors (h_n - h̄_n). For example, at `Δ = 1` in the
Ising identity sector the two states `(h, h̄) = (1, 0)` and `(0, 1)` are each counted.

For OBC (strip geometry, fixed-free boundary conditions):
    E_n - E_0 = (π v / L) h_n

where `h_n` are the boundary operator scaling dimensions. The three lowest boundary
primary dimensions of the Ising model are `h ∈ {0, 1/16, 1/2}` corresponding to the
identity, spin (fixed-free BC), and energy boundary operators (Cardy 1986, Table 2).

# Arguments
- `L::Real`: system size (must be > 0).
- `v::Real`: CFT sound velocity (model-dependent; must be > 0).

# Returns
A `Vector{NamedTuple{(:energy, :dimension, :degeneracy), Tuple{Float64, Float64, Int}}}`
sorted by ascending `energy`.

# References
- J. Cardy, *Nucl. Phys. B* **270**, 186 (1986). — operator content and strip spectra.
- H. W. J. Blöte, J. L. Cardy, M. P. Nightingale, *Phys. Rev. Lett.* **56**, 742 (1986).
- I. Affleck, *Phys. Rev. Lett.* **56**, 746 (1986). — SU(2)_1 WZW spectrum.
"""
function fetch(
    u::Universality{C}, q::ConformalTower, bc::Union{PBC,OBC}; L::Real, v::Real, kwargs...
) where {C}
    L > 0 || throw(ArgumentError("ConformalTower: L must be positive; got $L"))
    v > 0 || throw(ArgumentError("ConformalTower: v must be positive; got $v"))
    return _conformal_tower(u, q, bc, L, v; kwargs...)
end

# Fallback for other boundary conditions or classes
function _conformal_tower(
    u::Universality{C}, ::ConformalTower, bc::BoundaryCondition, L::Real, v::Real; kwargs...
) where {C}
    return error(
        "QAtlas Universality{:$C}: ConformalTower is not implemented at $(typeof(bc)) boundary condition.",
    )
end

# Ising PBC (torus, M(4,3) c=1/2): lowest scaling dimensions Δ ∈ {0, 1/8, 1, 2}
# Primary operators: identity (Δ=0), spin σ (Δ=1/8=0.125), energy ε (Δ=1).
# At Δ=2: degeneracy=2 counts the two spin-±1 level-1 descendants of the identity,
# (h, h̄) = (1, 0) and (0, 1). The ε sector also contributes descendants at Δ=2
# (the level-1 ε-sector tower, Δ=1+1) but these are in a separate sector with
# different quantum numbers.
function _conformal_tower(
    ::Universality{:Ising}, ::ConformalTower, ::PBC, L::Real, v::Real; kwargs...
)
    scale = 2 * π * v / L
    return NamedTuple{(:energy, :dimension, :degeneracy),Tuple{Float64,Float64,Int}}[
        (energy=0.0, dimension=0.0, degeneracy=1),
        (energy=0.125 * scale, dimension=0.125, degeneracy=1),
        (energy=1.0 * scale, dimension=1.0, degeneracy=1),
        (energy=2.0 * scale, dimension=2.0, degeneracy=2),
    ]
end

# Ising OBC (strip, fixed-free boundary conditions): boundary h ∈ {0, 1/16, 1/2}
# These are the three lowest boundary primary dimensions from Cardy (1986) Table 2
# for the fixed-free (i-f) strip: identity (h=0), spin boundary (h=1/16), energy
# boundary (h=1/2). Note: free-free (f-f) BC gives h ∈ {0, 1/2, 2, ...} (no 1/16).
# The TFIM with OBC endpoints maps to fixed-free BC via the boundary-state formalism.
function _conformal_tower(
    ::Universality{:Ising}, ::ConformalTower, ::OBC, L::Real, v::Real; kwargs...
)
    scale = π * v / L
    return NamedTuple{(:energy, :dimension, :degeneracy),Tuple{Float64,Float64,Int}}[
        (energy=0.0, dimension=0.0, degeneracy=1),
        (energy=0.0625 * scale, dimension=0.0625, degeneracy=1),
        (energy=0.5 * scale, dimension=0.5, degeneracy=1),
    ]
end

# Heisenberg PBC (torus, SU(2)_1 WZW, c=1): primaries j=0 (Δ=0) and j=1/2 (Δ=1/4).
# SU(2)_k WZW primary dimensions: h_j = j(j+1)/(k+2). For k=1: h_{1/2} = (3/4)/3 = 1/4.
# j=1 is not a primary (allowed primaries satisfy j ≤ k/2 = 1/2 for k=1).
# Degeneracies: j=1/2 sector has (2j+1)^2 = 4 primary states (both chiral sectors);
# at Δ=1: 3 left-current modes × 3 right-current modes = 9 descendants of the j=0 vacuum.
function _conformal_tower(
    ::Universality{:Heisenberg}, ::ConformalTower, ::PBC, L::Real, v::Real; kwargs...
)
    scale = 2 * π * v / L
    return NamedTuple{(:energy, :dimension, :degeneracy),Tuple{Float64,Float64,Int}}[
        (energy=0.0, dimension=0.0, degeneracy=1),
        (energy=0.25 * scale, dimension=0.25, degeneracy=4),
        (energy=1.0 * scale, dimension=1.0, degeneracy=9),
    ]
end
