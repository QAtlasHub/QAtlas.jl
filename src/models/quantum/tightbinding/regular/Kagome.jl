# ─────────────────────────────────────────────────────────────────────────────
# Kagome — nearest-neighbor tight-binding on the kagome lattice
#
# Hamiltonian (single-particle, spinless):
#   H = -t Σ_{⟨i,j⟩} (c†_i c_j + c†_j c_i)
#
# The kagome lattice has three sublattices (A, B, C). In the sublattice
# basis, the Bloch Hamiltonian is a 3×3 real symmetric matrix:
#
#   H(k) = -2t · [ 0              cos(θ₁/2)           cos(θ₂/2)        ]
#                [ cos(θ₁/2)      0                   cos((θ₂−θ₁)/2)   ]
#                [ cos(θ₂/2)      cos((θ₂−θ₁)/2)      0                ]
#
# with θ₁ = k·a₁, θ₂ = k·a₂ and Lattice2D's primitive basis
# a₁ = (1, 0), a₂ = (1/2, √3/2). Each off-diagonal element collects both
# the "same-cell" and "neighbouring-cell" hoppings of an A–B, A–C or B–C
# bond (two bonds per pair), yielding the factor 2.
#
# The key feature: one of the three bands is completely flat at E = +2t.
# The lower two bands touch it at the Γ-point (k = 0), where the spectrum
# degenerates to {−4t, +2t, +2t}.
#
# References:
#   I. Syôzi, "Statistics of Kagomé Lattice", Prog. Theor. Phys. 6, 306 (1951).
#   D. L. Bergman, C. Wu, L. Balents, "Band touching from real-space
#     topology in frustrated hopping models", Phys. Rev. B 78, 125104 (2008).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Fermion bilinears c†c
#   Observable:  Fermion (number n = c†c, bilinear ⟨c†_i c_j⟩); derived spin observables follow spin S = σ/2
#   Reference:   docs/src/conventions.md §Fermion convention

using LinearAlgebra: Symmetric, eigvals

# ═══════════════════════════════════════════════════════════════════════════════
# Dispatch tag
# ═══════════════════════════════════════════════════════════════════════════════

"""
    Kagome

Dispatch tag for the nearest-neighbor tight-binding model on the kagome
lattice (spinless, single-orbital, three sublattices).

Hamiltonian: H = -t Σ_{⟨i,j⟩} (c†_i c_j + h.c.)

The spectrum exhibits a dispersionless flat band at E = +2t touching the
lower dispersive bands at the Γ-point.
"""
struct Kagome <: AbstractQAtlasModel
    ;
    t::Float64;
    Lx::Int;
    Ly::Int;
end
Kagome(; t::Real=1.0, Lx::Integer=0, Ly::Integer=0) = Kagome(Float64(t), Int(Lx), Int(Ly))

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: closed-form Bloch spectrum for Lx × Ly kagome PBC
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::Kagome, ::TightBindingSpectrum; Lx, Ly, t=1.0) -> Vector{Float64}

Sorted single-particle spectrum of the nearest-neighbor tight-binding
Hamiltonian on an Lx × Ly kagome lattice with periodic boundary
conditions in both directions.

The 3×3 Bloch Hamiltonian is diagonalized numerically at each of the
`Lx·Ly` allowed momenta

    k_{mn} = (m/Lx) b₁ + (n/Ly) b₂,   m ∈ {0,…,Lx−1}, n ∈ {0,…,Ly−1}

with b₁, b₂ dual to `Lattice2D`'s kagome primitive basis. In the θ
parametrization (θ₁ = k·a₁ = 2π m/Lx, θ₂ = k·a₂ = 2π n/Ly),

    H(k) = -2t · [ 0              cos(θ₁/2)           cos(θ₂/2)        ]
                 [ cos(θ₁/2)      0                   cos((θ₂−θ₁)/2)   ]
                 [ cos(θ₂/2)      cos((θ₂−θ₁)/2)      0                ]

# Flat band

One band is exactly flat at `+2t`, so the returned spectrum contains at
least `Lx·Ly` eigenvalues equal to `2t`. At the Γ-point (k = 0) the
upper dispersive band touches the flat band, contributing one extra
`+2t` eigenvalue, giving a `Lx·Ly + 1`-fold degeneracy there. The
remaining `2·Lx·Ly − 1` eigenvalues lie in `[−4t, +2t)`.

# Arguments
- `Lx::Int`: number of unit cells in the first lattice direction
- `Ly::Int`: number of unit cells in the second lattice direction
- `t::Real`: nearest-neighbor hopping amplitude (default 1.0)

# Return
A sorted `Vector{Float64}` of length `3·Lx·Ly`.

# References
    I. Syôzi, Prog. Theor. Phys. 6, 306 (1951).
    D. L. Bergman et al., Phys. Rev. B 78, 125104 (2008).
"""
function fetch(
    m::Kagome, ::TightBindingSpectrum; Lx::Integer=m.Lx, Ly::Integer=m.Ly, t::Real=m.t
)
    Lx > 0 && Ly > 0 || error("Kagome TightBindingSpectrum: Lx, Ly must be positive")
    eigs = Float64[]
    sizehint!(eigs, 3 * Lx * Ly)
    for m in 0:(Lx - 1), n in 0:(Ly - 1)
        θ1 = 2π * m / Lx
        θ2 = 2π * n / Ly
        c1 = cos(θ1 / 2)
        c2 = cos(θ2 / 2)
        c12 = cos((θ2 - θ1) / 2)
        H = -2 * t * [0.0 c1 c2; c1 0.0 c12; c2 c12 0.0]
        append!(eigs, eigvals(Symmetric(H)))
    end
    return sort!(eigs)
end
