# ─────────────────────────────────────────────────────────────────────────────
# Honeycomb — nearest-neighbor tight-binding on the honeycomb lattice
#
# Also known as the graphene band structure (Wallace 1947).  The model
# is named after the lattice rather than the material so it composes
# cleanly with Lattice2D (which already uses the `Graphene` name for
# its topology type).  `const Graphene = Honeycomb` in
# `src/deprecate/legacy_honeycomb.jl` keeps pre-v0.13 callers working.
#
# Hamiltonian (single-particle, spinless):
#   H = -t Σ_{⟨i,j⟩ ∈ A-B} (c†_i c_j + c†_j c_i)
#
# The Bloch Hamiltonian on a bipartite (A, B) basis is:
#
#   H(k) = [  0      f(k)  ]
#          [ f(k)*   0     ]
#
# with f(k) = -t [exp(i k·δ₁) + exp(i k·δ₂) + exp(i k·δ₃)], where δᵢ are
# the three nearest-neighbor displacement vectors from an A site to its
# three B neighbors. The two bands are ±|f(k)|.
#
# Using the unit-cell basis (a₁, a₂) adopted by `Lattice2D`'s honeycomb
# topology (a₁ = (√3, 0), a₂ = (√3/2, 3/2), δ's derived accordingly):
#
#   |f(k)|² = t² [3 + 2cos(k·a₁) + 2cos(k·a₂) + 2cos(k·(a₂−a₁))]
#
# References:
#   P. R. Wallace, "The Band Theory of Graphite", Phys. Rev. 71, 622 (1947).
#   A. H. Castro Neto et al., "The electronic properties of graphene",
#     Rev. Mod. Phys. 81, 109 (2009).
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Dispatch tags
# ═══════════════════════════════════════════════════════════════════════════════

# CONVENTION
#   Hamiltonian: Fermion bilinears c†c
#   Observable:  Fermion (number n = c†c, bilinear ⟨c†_i c_j⟩); derived spin observables follow spin S = σ/2
#   Reference:   docs/src/conventions.md §Fermion convention

"""
    Honeycomb(; t::Real = 1.0, Lx::Int = 0, Ly::Int = 0) <: AbstractQAtlasModel

Nearest-neighbor tight-binding model on the honeycomb lattice
(spinless, single-orbital) — historically known as graphene.

Hamiltonian: `H = -t Σ_{⟨i,j⟩} (c†_i c_j + h.c.)`.

Fields:
- `t`  — nearest-neighbor hopping amplitude (default `1.0`).
- `Lx`, `Ly` — unit-cell counts in the two lattice directions.  `0` is
  a sentinel meaning "caller passes via kwargs at fetch time".

For backward compatibility with the pre-v0.13 call form
`fetch(Graphene(), TightBindingSpectrum(); Lx, Ly, t=…)`, the fetch
method still accepts `Lx`, `Ly`, `t` as kwargs (they override the
struct fields when supplied).  The legacy name `Graphene` is retained
as a type alias in `src/deprecate/legacy_honeycomb.jl`.
"""
struct Honeycomb <: AbstractQAtlasModel
    t::Float64
    Lx::Int
    Ly::Int
end
function Honeycomb(; t::Real=1.0, Lx::Integer=0, Ly::Integer=0)
    Honeycomb(Float64(t), Int(Lx), Int(Ly))
end

"""
    TightBindingSpectrum() <: AbstractQuantity

Single-particle Bloch spectrum of a tight-binding model.  Returned as
a sorted `Vector{Float64}` of length `n_orbitals · Lx · Ly`.
"""
struct TightBindingSpectrum <: AbstractQuantity end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: closed-form Bloch spectrum for Lx × Ly honeycomb PBC
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::Graphene, ::TightBindingSpectrum; Lx, Ly, t=1.0) -> Vector{Float64}

Sorted single-particle spectrum of the nearest-neighbor tight-binding
Hamiltonian on an Lx × Ly honeycomb lattice with periodic boundary
conditions in both directions.

The spectrum is obtained in closed form by diagonalizing the 2×2 Bloch
Hamiltonian at the `Lx·Ly` allowed momenta

    k_{mn} = (m/Lx) b₁ + (n/Ly) b₂,   m ∈ {0,…,Lx−1}, n ∈ {0,…,Ly−1}

where (b₁, b₂) is the reciprocal basis dual to the real-space
basis (a₁, a₂) used by `Lattice2D`'s `Honeycomb` topology.

Explicitly, using k·a₁ = 2π m/Lx and k·a₂ = 2π n/Ly,

    E_{mn,±} = ± t · √(3 + 2cos(2π m/Lx) + 2cos(2π n/Ly) + 2cos(2π(n/Ly − m/Lx)))

The chiral (sublattice) symmetry of the bipartite honeycomb lattice
guarantees that the spectrum is symmetric about zero. A small negative
floor is applied inside the square root to absorb rounding error at
Dirac points where the argument vanishes exactly.

# Arguments
- `Lx::Int`: number of unit cells in the first lattice direction
- `Ly::Int`: number of unit cells in the second lattice direction
- `t::Real`: nearest-neighbor hopping amplitude (default 1.0)

# Return
A sorted `Vector{Float64}` of length `2·Lx·Ly`.

# References
    P. R. Wallace, Phys. Rev. 71, 622 (1947).
    A. H. Castro Neto et al., Rev. Mod. Phys. 81, 109 (2009).
"""
function fetch(
    m::Honeycomb, ::TightBindingSpectrum; Lx::Integer=m.Lx, Ly::Integer=m.Ly, t::Real=m.t
)
    Lx > 0 && Ly > 0 || error(
        "Honeycomb TightBindingSpectrum: Lx and Ly must be positive. " *
        "Pass them in the struct (Honeycomb(; Lx, Ly)) or as kwargs.",
    )
    eigs = Float64[]
    sizehint!(eigs, 2 * Lx * Ly)
    for mi in 0:(Lx - 1), ni in 0:(Ly - 1)
        θ1 = 2π * mi / Lx
        θ2 = 2π * ni / Ly
        # |f(k)|² / t² = 3 + 2cos(θ1) + 2cos(θ2) + 2cos(θ2 - θ1)
        val = 3 + 2cos(θ1) + 2cos(θ2) + 2cos(θ2 - θ1)
        energy = abs(t) * sqrt(max(val, 0.0))
        push!(eigs, -energy)
        push!(eigs, energy)
    end
    return sort!(eigs)
end

# ─── Scalar invariants for verify() integration ────────────────────────────
#
# The full Bloch spectrum is a Vector{Float64}; verify() works on scalars.
# Exposing a chiral-symmetry-pinned scalar reduction lets cross-validate
# the closed-form spectrum against (a) the analytical tr(H²) identity and
# (b) the real-space ED ground truth (via Lattice2D.build_tight_binding,
# in test/util/tight_binding.jl), without leaving the verify() framework.

"""
    TightBindingChecksum() <: AbstractQuantity

Sum of squared single-particle eigenvalues `Σ λᵢ²` of the
tight-binding Bloch Hamiltonian, i.e. `tr(H²)`.  For a chiral
(bipartite) lattice this equals `2 t² · n_NN_bonds`, an analytical
identity independent of `Lx`, `Ly`, or k-grid details — a strong
sanity invariant on top of the spectrum closed form.
"""
struct TightBindingChecksum <: AbstractQuantity end

"""
    fetch(::Honeycomb, ::TightBindingChecksum, ::Infinite; Lx, Ly, t=1.0) -> Float64

`Σ λᵢ² = tr(H²)` for the Lx × Ly honeycomb PBC tight-binding spectrum.
Forwards through TightBindingSpectrum so it always agrees with the
spectrum closed form by construction; verify() cards then pin it
against the chiral-symmetry identity `2 t² · 3 · Lx · Ly` and against
the real-space ED ground truth.
"""
function fetch(
    m::Honeycomb,
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
    TightBindingMaxEnergy() <: AbstractQuantity

Largest single-particle eigenvalue of the tight-binding Bloch
Hamiltonian, i.e. `max(λᵢ)`.  For the chiral honeycomb this equals
`|t| · max_k |f(k)| = 3|t|` (saturated at the Γ point where the three
nearest-neighbor phases add coherently), independent of `Lx` and `Ly`
as long as the Γ point is on the discrete momentum grid (always true
for PBC with `Lx, Ly ≥ 1`).

Sister scalar to `TightBindingChecksum` (which captures `Σλᵢ² = tr(H²)`
and is invariant under chiral pair sign-flips).  The two together
discriminate any single-eigenvalue perturbation: a sign error on a
hopping displacement that flips a chiral pair changes neither the
checksum nor the bandwidth at Γ only if it doesn't reach the Γ point —
making them jointly sensitive to a broader class of typos.
"""
struct TightBindingMaxEnergy <: AbstractQuantity end

"""
    fetch(::Honeycomb, ::TightBindingMaxEnergy, ::Infinite; Lx, Ly, t=1.0) -> Float64

`max(λᵢ)` for the Lx × Ly honeycomb PBC tight-binding spectrum.
Forwards through `TightBindingSpectrum`.
"""
function fetch(
    m::Honeycomb,
    ::TightBindingMaxEnergy,
    ::Infinite;
    Lx::Integer=m.Lx,
    Ly::Integer=m.Ly,
    t::Real=m.t,
)
    eigs = fetch(m, TightBindingSpectrum(); Lx=Lx, Ly=Ly, t=t)
    return maximum(eigs)
end
