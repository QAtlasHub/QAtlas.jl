# ─────────────────────────────────────────────────────────────────────────────
# Lieb — nearest-neighbor tight-binding on the Lieb lattice
#
# Hamiltonian (single-particle, spinless):
#   H = -t Σ_{⟨i,j⟩} (c†_i c_j + c†_j c_i)
#
# The Lieb lattice (line-centred square) has three sublattices per unit
# cell: A (corner) and B, C (right-edge and top-edge). Only A–B and A–C
# bonds exist — there is no direct B–C bond — so the lattice is
# bipartite between {A} and {B, C}.
#
# In the sublattice basis the Bloch Hamiltonian is:
#
#   H(k) = -2t · [  0              cos(θ₁/2)     cos(θ₂/2)  ]
#                [  cos(θ₁/2)      0             0          ]
#                [  cos(θ₂/2)      0             0          ]
#
# with θ₁ = k·a₁, θ₂ = k·a₂ and Lattice2D's primitive basis
# a₁ = (2, 0), a₂ = (0, 2). Each off-diagonal element collects the
# "same-cell" and "neighbouring-cell" A–B or A–C hopping (factor 2).
#
# The structure of H(k) admits an exact closed-form eigenvalue:
#
#   spec = { −E(k), 0, +E(k) },
#   E(k) = 2 |t| √(cos²(θ₁/2) + cos²(θ₂/2))
#
# so one band is dispersionless at E = 0 (the Lieb flat band, a
# consequence of the bipartite imbalance |A| − (|B| + |C|) per cell).
#
# At the M-point (θ₁ = θ₂ = π) both dispersive bands also collapse to
# zero, giving a three-fold band touching at E = 0.
#
# References:
#   E. H. Lieb, "Two Theorems on the Hubbard Model",
#     Phys. Rev. Lett. 62, 1201 (1989).
#   H. Tasaki, "From Nagaoka's Ferromagnetism to Flat-Band Ferromagnetism
#     and Beyond", Prog. Theor. Phys. 99, 489 (1998).
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Dispatch tag
# ═══════════════════════════════════════════════════════════════════════════════

"""
    Lieb

Dispatch tag for the nearest-neighbor tight-binding model on the Lieb
(line-centred square) lattice (spinless, single-orbital, three
sublattices). Hamiltonian:

    H = -t Σ_{⟨i,j⟩} (c†_i c_j + h.c.)

The spectrum contains a dispersionless flat band at `E = 0` for every
momentum, with an additional three-fold band touching at the M-point.
"""
struct Lieb <: AbstractQAtlasModel
    ;
    t::Float64;
    Lx::Int;
    Ly::Int;
end
Lieb(; t::Real=1.0, Lx::Integer=0, Ly::Integer=0) = Lieb(Float64(t), Int(Lx), Int(Ly))

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: closed-form Bloch spectrum for Lx × Ly Lieb PBC
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::Lieb, ::TightBindingSpectrum; Lx, Ly, t=1.0) -> Vector{Float64}

Sorted single-particle spectrum of the nearest-neighbor tight-binding
Hamiltonian on an `Lx × Ly` Lieb lattice with periodic boundary
conditions in both directions.

The closed-form eigenvalues of the 3×3 Bloch Hamiltonian are

    E_{mn,±} = ± 2t · √(cos²(π m/Lx) + cos²(π n/Ly)),
    E_{mn,0} = 0

for m ∈ {0,…,Lx−1}, n ∈ {0,…,Ly−1}. The `E = 0` contribution appears
for every (m, n) and forms the Lieb flat band.

# Flat band

For a generic `(Lx, Ly)` the spectrum contains exactly `Lx·Ly` zero
eigenvalues from the flat band. When both `Lx` and `Ly` are even the
M-point `(θ₁, θ₂) = (π, π)` becomes an allowed momentum; there the two
dispersive bands also collapse to zero, contributing two *extra* zero
modes (a three-fold band touching). The total number of `E = 0`
eigenvalues is therefore

    Lx·Ly + (2 if both Lx and Ly are even else 0).

# Arguments
- `Lx::Int`: number of unit cells in the first lattice direction
- `Ly::Int`: number of unit cells in the second lattice direction
- `t::Real`: nearest-neighbor hopping amplitude (default 1.0)

# Return
A sorted `Vector{Float64}` of length `3·Lx·Ly`.

# References
    E. H. Lieb, Phys. Rev. Lett. 62, 1201 (1989).
"""
function fetch(
    m::Lieb, ::TightBindingSpectrum; Lx::Integer=m.Lx, Ly::Integer=m.Ly, t::Real=m.t
)
    Lx > 0 && Ly > 0 || error("Lieb TightBindingSpectrum: Lx, Ly must be positive")
    eigs = Float64[]
    sizehint!(eigs, 3 * Lx * Ly)
    for m in 0:(Lx - 1), n in 0:(Ly - 1)
        θ1 = 2π * m / Lx
        θ2 = 2π * n / Ly
        c1 = cos(θ1 / 2)
        c2 = cos(θ2 / 2)
        energy = 2 * abs(t) * sqrt(c1 * c1 + c2 * c2)
        push!(eigs, -energy)
        push!(eigs, 0.0)
        push!(eigs, energy)
    end
    return sort!(eigs)
end

# ─── Scalar invariants for verify() integration ────────────────────────────
#
# Sister to TightBindingChecksum + TightBindingMaxEnergy on Honeycomb.
# Quantities defined in Honeycomb.jl; this file forwards them through
# Lieb's TightBindingSpectrum.

"""
    fetch(::Lieb, ::TightBindingChecksum, ::Infinite; Lx, Ly, t=1.0) -> Float64

`Σ λᵢ² = tr(H²)` for the Lx × Ly Lieb-lattice PBC TB spectrum. Lieb has
3 sites per unit cell (1 corner + 2 edge) with coordinations 4, 2, 2;
`n_NN_bonds = 4 · Lx · Ly` and `tr(H²) = 8 · t² · Lx · Ly`.
"""
function fetch(
    m::Lieb,
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
    fetch(::Lieb, ::TightBindingMaxEnergy, ::Infinite; Lx, Ly, t=1.0) -> Float64

`max(λᵢ)` for the Lx × Ly Lieb PBC TB spectrum. The dispersive bands
`±E(k) = ±2|t|·√(cos²(θ₁/2) + cos²(θ₂/2))` saturate at the Γ point
(θ₁=θ₂=0), giving `max(spectrum) = 2√2 · |t|`, independent of Lx, Ly.
"""
function fetch(
    m::Lieb,
    ::TightBindingMaxEnergy,
    ::Infinite;
    Lx::Integer=m.Lx,
    Ly::Integer=m.Ly,
    t::Real=m.t,
)
    eigs = fetch(m, TightBindingSpectrum(); Lx=Lx, Ly=Ly, t=t)
    return maximum(eigs)
end
