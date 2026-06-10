# universalities/behaviour/conformal_towers.jl
#
# Universal BEHAVIOUR: Conformal tower of states excitation spectrum in 1+1D CFT.
#
# References:
#   J. Cardy, Nucl. Phys. B 270, 186 (1986).
#   H. W. J. Blöte, J. L. Cardy, M. P. Nightingale, Phys. Rev. Lett. 56, 742 (1986).

raw"""
    fetch(::Universality{C}, ::ConformalTower, bc::Union{PBC, OBC}; L::Real, v::Real) -> Vector{NamedTuple}

Return a sorted vector of NamedTuples representing the lowest-lying excitation energies `E_n - E_0`,
their scaling dimensions `Δ_n` (or `h_n`), and their degeneracies in the conformal tower of states:

    (energy = E_n - E_0, dimension = Δ_n, degeneracy = g_n)

For PBC:
    E_n - E_0 = (2π v / L) Δ_n
For OBC:
    E_n - E_0 = (π v / L) h_n
"""
function fetch(
    u::Universality{C},
    q::ConformalTower,
    bc::Union{PBC, OBC};
    L::Real,
    v::Real,
    kwargs...,
) where {C}
    L > 0 || throw(ArgumentError("ConformalTower: L must be positive; got $L"))
    v > 0 || throw(ArgumentError("ConformalTower: v must be positive; got $v"))
    return _conformal_tower(u, q, bc, L, v; kwargs...)
end

# Fallback for other boundary conditions or classes
function _conformal_tower(
    u::Universality{C},
    ::ConformalTower,
    bc::BoundaryCondition,
    L::Real,
    v::Real;
    kwargs...,
) where {C}
    return error(
        "QAtlas Universality{:$C}: ConformalTower is not implemented at $(typeof(bc)) boundary condition.",
    )
end

# Ising PBC (bulk states): Δ ∈ {0, 0.125, 1.0, 2.0}
function _conformal_tower(
    ::Universality{:Ising},
    ::ConformalTower,
    ::PBC,
    L::Real,
    v::Real;
    kwargs...,
)
    scale = 2 * π * v / L
    return NamedTuple{(:energy, :dimension, :degeneracy),Tuple{Float64,Float64,Int}}[
        (energy=0.0, dimension=0.0, degeneracy=1),
        (energy=0.125 * scale, dimension=0.125, degeneracy=1),
        (energy=1.0 * scale, dimension=1.0, degeneracy=1),
        (energy=2.0 * scale, dimension=2.0, degeneracy=2),
    ]
end

# Ising OBC (boundary free-free states): h ∈ {0, 0.0625, 0.5}
function _conformal_tower(
    ::Universality{:Ising},
    ::ConformalTower,
    ::OBC,
    L::Real,
    v::Real;
    kwargs...,
)
    scale = π * v / L
    return NamedTuple{(:energy, :dimension, :degeneracy),Tuple{Float64,Float64,Int}}[
        (energy=0.0, dimension=0.0, degeneracy=1),
        (energy=0.0625 * scale, dimension=0.0625, degeneracy=1),
        (energy=0.5 * scale, dimension=0.5, degeneracy=1),
    ]
end

# Heisenberg PBC (bulk SU(2)_1 WZW states): Δ ∈ {0, 0.25, 1.0}
function _conformal_tower(
    ::Universality{:Heisenberg},
    ::ConformalTower,
    ::PBC,
    L::Real,
    v::Real;
    kwargs...,
)
    scale = 2 * π * v / L
    return NamedTuple{(:energy, :dimension, :degeneracy),Tuple{Float64,Float64,Int}}[
        (energy=0.0, dimension=0.0, degeneracy=1),
        (energy=0.25 * scale, dimension=0.25, degeneracy=4),
        (energy=1.0 * scale, dimension=1.0, degeneracy=9),
    ]
end
