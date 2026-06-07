# ─────────────────────────────────────────────────────────────────────────────
# Triangular — nearest-neighbor tight-binding on the triangular lattice
#
# Hamiltonian (single-particle, spinless):
#   H = -t Σ_{⟨i,j⟩} (c†_i c_j + c†_j c_i)
#
# The triangular lattice has one sublattice per unit cell with six
# nearest neighbours per site. Since there is only one band, the Bloch
# "Hamiltonian" is a scalar:
#
#   E(k) = -t · Σ_{δ ∈ 6 NN} exp(i k·δ)
#        = -2t · [ cos(k·a₁) + cos(k·a₂) + cos(k·(a₂ − a₁)) ]
#
# with Lattice2D's primitive basis  a₁ = (1, 0),  a₂ = (1/2, √3/2).
#
# Frustration: unlike the (bipartite) square or honeycomb lattice, the
# triangular band is *not* symmetric about zero. With J ≡ t > 0,
#
#   minimum  E_min = -6t   at the Γ-point (unique),
#   maximum  E_max = +3t   at the K-points (m/Lx = n/Ly = 1/3 etc.),
#
# i.e. the band extends from -6t to +3t, giving the hallmark
# Van Hove singularity and the asymmetric density of states of the
# frustrated triangular hopping problem.
#
# References:
#   G. H. Wannier, "Antiferromagnetism. The Triangular Ising Net",
#     Phys. Rev. 79, 357 (1950).
#   T. Koretsune and M. Ogata, "Electronic structures of triangular
#     lattice models", J. Phys. Soc. Jpn. 76, 074706 (2007) — review
#     of the NN tight-binding spectrum and Van Hove singularity.
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Dispatch tag
# ═══════════════════════════════════════════════════════════════════════════════

# CONVENTION
#   Hamiltonian: Fermion bilinears c†c
#   Observable:  Fermion (number n = c†c, bilinear ⟨c†_i c_j⟩); derived spin observables follow spin S = σ/2
#   Reference:   docs/src/conventions.md §Fermion convention

"""
    Triangular

Dispatch tag for the nearest-neighbor tight-binding model on the
triangular lattice (spinless, single-orbital, one sublattice).

Hamiltonian: H = -t Σ_{⟨i,j⟩} (c†_i c_j + h.c.)

The single band ranges from `-6t` (at the Γ-point) to `+3t` (at the
K-points), reflecting geometric frustration: the absence of
particle-hole (chiral) symmetry on a non-bipartite lattice.
"""
struct Triangular <: AbstractQAtlasModel
    ;
    t::Float64;
    Lx::Int;
    Ly::Int;
end
function Triangular(; t::Real=1.0, Lx::Integer=0, Ly::Integer=0)
    return Triangular(Float64(t), Int(Lx), Int(Ly))
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: closed-form Bloch spectrum for Lx × Ly triangular PBC
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::Triangular, ::TightBindingSpectrum; Lx, Ly, t=1.0) -> Vector{Float64}

Sorted single-particle spectrum of the nearest-neighbor tight-binding
Hamiltonian on an `Lx × Ly` triangular lattice with periodic boundary
conditions in both directions.

Because the triangular lattice has a single sublattice per unit cell,
each allowed momentum contributes exactly one eigenvalue

    E_{mn} = -2t · [ cos(2π m/Lx) + cos(2π n/Ly) + cos(2π(n/Ly − m/Lx)) ]

for m ∈ {0,…,Lx−1}, n ∈ {0,…,Ly−1}. Note the lower cut-off is `-6t`
(at the Γ-point) but the upper cut-off is only `+3t` (at the K-points
`(1/3, 2/3)` and `(2/3, 1/3)` — allowed when `Lx` and `Ly` are
divisible by 3).

# Arguments
- `Lx::Int`: number of unit cells in the first lattice direction
- `Ly::Int`: number of unit cells in the second lattice direction
- `t::Real`: nearest-neighbor hopping amplitude (default 1.0; positive
  convention)

# Return
A sorted `Vector{Float64}` of length `Lx·Ly`.

# References
    G. H. Wannier, Phys. Rev. 79, 357 (1950).
"""
function fetch(
    m::Triangular, ::TightBindingSpectrum; Lx::Integer=m.Lx, Ly::Integer=m.Ly, t::Real=m.t
)
    Lx > 0 && Ly > 0 || error("Triangular TightBindingSpectrum: Lx, Ly must be positive")
    eigs = Vector{Float64}(undef, Lx * Ly)
    idx = 1
    for m in 0:(Lx - 1), n in 0:(Ly - 1)
        θ1 = 2π * m / Lx
        θ2 = 2π * n / Ly
        eigs[idx] = -2 * t * (cos(θ1) + cos(θ2) + cos(θ2 - θ1))
        idx += 1
    end
    return sort!(eigs)
end

# ─── Scalar invariants for verify() integration ────────────────────────────
#
# Sister to TightBindingChecksum + TightBindingMaxEnergy on Honeycomb.
# Quantities defined in Honeycomb.jl; forwarded here for Triangular.
#
# Note: Triangular is NOT bipartite. max(spectrum) saturates only at the
# K-point (θ₁=θ₂=2π/3), which is on the Lx × Ly k-grid only when both
# Lx, Ly are multiples of 3. So MaxEnergy is exposed for cross-check
# against real-space ED but has no Lx,Ly-independent literature value.
# Checksum tr(H²) = 2 t² · n_NN_bonds = 6 t² Lx Ly always holds.

"""
    fetch(::Triangular, ::TightBindingChecksum, ::Infinite; Lx, Ly, t=1.0) -> Float64

`Σ λᵢ² = tr(H²)` for the Lx × Ly triangular PBC TB spectrum. Each site
has 6 NN ⇒ `n_NN_bonds = 3 · Lx · Ly` and `tr(H²) = 6 · t² · Lx · Ly`.
"""
function fetch(
    m::Triangular,
    ::TightBindingChecksum,
    ::Infinite;
    Lx::Integer=m.Lx,
    Ly::Integer=m.Ly,
    t::Real=m.t,
)
    eigs = fetch(m, TightBindingSpectrum(); Lx=Lx, Ly=Ly, t=t)
    return sum(abs2, eigs)
end

"""
    fetch(::Triangular, ::TightBindingMaxEnergy, ::Infinite; Lx, Ly, t=1.0) -> Float64

`max(λᵢ)` for the Lx × Ly triangular PBC TB spectrum. Maximum reaches
`+3 t` only when the K-point is on the discrete grid (both Lx, Ly
multiples of 3); otherwise it is strictly less. Exposed for ED cross-
check rather than as a closed-form literature pin.
"""
function fetch(
    m::Triangular,
    ::TightBindingMaxEnergy,
    ::Infinite;
    Lx::Integer=m.Lx,
    Ly::Integer=m.Ly,
    t::Real=m.t,
)
    eigs = fetch(m, TightBindingSpectrum(); Lx=Lx, Ly=Ly, t=t)
    return maximum(eigs)
end
