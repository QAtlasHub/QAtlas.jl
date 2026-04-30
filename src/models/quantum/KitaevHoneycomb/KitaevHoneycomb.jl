# ─────────────────────────────────────────────────────────────────────────────
# Kitaev honeycomb model — exact solutions
#
# Hamiltonian:
#   H = − Σ_{⟨ij⟩ ∈ x-bonds} Kₓ σˣᵢσˣⱼ
#       − Σ_{⟨ij⟩ ∈ y-bonds} Kᵧ σʸᵢσʸⱼ
#       − Σ_{⟨ij⟩ ∈ z-bonds} K_z σᶻᵢσᶻⱼ
#
# on the 2D honeycomb lattice (Lieb/Kitaev convention). The three bond
# types are assigned to the three nearest-neighbor directions emitted
# by Lattice2D's `Honeycomb` topology:
#
#   :type_1 ↔ z-bond (A→B same cell)
#   :type_2 ↔ x-bond (A→B upper-left cell, (+a₁, −a₂))
#   :type_3 ↔ y-bond (A→B upper cell, (−a₂))
#
# Kitaev (2006) showed that after a Jordan–Wigner–like four-Majorana
# mapping the gauge fields `u_{ij} = i bᵞᵢ bᵞⱼ` on each bond become
# conserved Z₂ variables, and the ground state sits in the flux-free
# sector (all plaquette fluxes `W_p = +1`, Lieb's theorem). In that
# sector the matter Majoranas decouple into a free hopping problem on
# the bipartite A/B lattice:
#
#   H_flux_free = (i/2) Σ_{⟨ij⟩, γ} Kᵞ c_i c_j,
#
# whose Bloch Hamiltonian on the two-sublattice basis is
#
#   H(k) = [  0     f(k)  ]
#          [ f(k)*   0    ],
#   f(k) = K_z + Kₓ exp(i θ₁) + Kᵧ exp(i θ₂)
#        (θ₁ = k·a₁, θ₂ = k·a₂; Lattice2D basis vectors).
#
# Bands ±|f(k)| fill the negative one at T=0. Per site (2 sites per
# unit cell):
#
#   |f|² = Kₓ² + Kᵧ² + K_z² + 2Kₓ Kᵧ cos θ₁
#                           + 2Kᵧ K_z cos θ₂
#                           + 2K_z Kₓ cos(θ₁ − θ₂).
#
# OBC of finite `Lx×Ly` uses the same flux-free-sector ansatz with
# `u_{ij} = +1`, giving a real `Lx·Ly × Lx·Ly` bipartite hopping
# matrix M whose singular values σ_k are the positive BdG single-
# particle energies. Ground state energy = −Σₖ σₖ; per site = energy
# divided by `2·Lx·Ly`.
#
# References
#   - A. Kitaev, "Anyons in an exactly solved model and beyond",
#     Ann. Phys. 321, 2 (2006).
#   - E. H. Lieb, PRL 73, 2158 (1994) — flux-free ground state.
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra: svdvals
using QuadGK: quadgk

"""
    KitaevHoneycomb(; Kx=1.0, Ky=1.0, Kz=1.0) <: AbstractQAtlasModel

Kitaev honeycomb model with coupling amplitudes on the three bond
types. `Kx = Ky = Kz` is the isotropic gapless point. See module header
for the full Hamiltonian and solution outline.
"""
struct KitaevHoneycomb <: AbstractQAtlasModel
    Kx::Float64
    Ky::Float64
    Kz::Float64
end
function KitaevHoneycomb(; Kx::Real=1.0, Ky::Real=1.0, Kz::Real=1.0)
    KitaevHoneycomb(Float64(Kx), Float64(Ky), Float64(Kz))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Internal: Bloch-form |f(k)|²
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _kitaev_fk_abs²(m, θ₁, θ₂) -> Float64

`|f(k)|² = Kₓ² + Kᵧ² + K_z² + 2Kₓ Kᵧ cos θ₁ + 2Kᵧ K_z cos θ₂ + 2K_z Kₓ cos(θ₁ − θ₂)`.
"""
@inline function _kitaev_fk_abs²(m::KitaevHoneycomb, θ₁::Real, θ₂::Real)
    return (
        m.Kx^2 +
        m.Ky^2 +
        m.Kz^2 +
        2 * m.Kx * m.Ky * cos(θ₁) +
        2 * m.Ky * m.Kz * cos(θ₂) +
        2 * m.Kz * m.Kx * cos(θ₁ - θ₂)
    )
end

@inline _kitaev_fk_abs(m::KitaevHoneycomb, θ₁::Real, θ₂::Real) = sqrt(
    max(_kitaev_fk_abs²(m, θ₁, θ₂), 0.0)
)

# ═══════════════════════════════════════════════════════════════════════════════
# Energy granularity convention: Kitaev is a 2D model — per-site is the
# native granularity at every BC.  The total-energy generic conversion in
# `src/core/quantities.jl` would need `_bc_size(bc, kwargs)` to read
# `Lx, Ly` from kwargs; we leave `Energy(:total)` unsupported for now and
# let it surface as the standard "no method registered" error.
# ═══════════════════════════════════════════════════════════════════════════════

native_energy_granularity(::KitaevHoneycomb, ::OBC) = :per_site
native_energy_granularity(::KitaevHoneycomb, ::PBC) = :per_site
native_energy_granularity(::KitaevHoneycomb, ::Infinite) = :per_site

# ═══════════════════════════════════════════════════════════════════════════════
# Energy — Infinite (per site)
#
# The unified GS + thermal `Energy{:per_site}, Infinite` fetch lives in
# `KitaevHoneycomb_thermal.jl` (sibling include).  At `beta = nothing`
# (default) it returns the T = 0 ground-state value
#
#     ε_gs = − 1/(8π²) ∫₀^{2π} dθ₁ ∫₀^{2π} dθ₂ |f(θ₁, θ₂)|
#
# computed by adaptive nested Gauss-Kronrod quadrature.  With a real
# `beta` it returns the matter-sector thermal value (see thermal file).
# Defining the method there keeps the GS path bit-for-bit identical
# while avoiding a method-overwrite warning on package load.
# ═══════════════════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════════════════
# Energy — PBC finite (per site)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _kitaev_pbc_sector_energy(m, Lx, Ly, νx, νy) -> Float64

Per-site ground state energy on the `Lx × Ly` torus **within a fixed
topological flux sector** `(νx, νy) ∈ {0, 1/2}²`. `νx = 1/2` flips
every bond that crosses the `m1 = Lx → m1 = 1` seam (anti-periodic
fermion boundary condition along the x-direction); `νy = 1/2` does the
same in the y-direction. The four sectors correspond to the four
inequivalent choices `(W_x, W_y) ∈ {±1}²` of the two non-contractible
Wilson-loop operators.

Built as a singular-value decomposition of the Majorana hopping matrix
`M` with entries `±K_γ`, where the sign flips whenever a bond of type
γ crosses one of the selected seams. Per-site energy = `−Σ σ_k / (2·Lx·Ly)`.
"""
function _kitaev_pbc_sector_energy(
    m::KitaevHoneycomb, Lx::Integer, Ly::Integer, νx::Int, νy::Int
)
    idx(m1, m2) = (m1 - 1) * Ly + m2
    N = Lx * Ly
    M = zeros(Float64, N, N)
    for m1 in 1:Lx, m2 in 1:Ly
        a = idx(m1, m2)
        # z-bond: (m1, m2) → (m1, m2), never crosses a seam.
        M[a, idx(m1, m2)] += m.Kz
        # x-bond: (m1, m2) → (m1 + 1 mod Lx, m2 − 1 mod Ly). Crosses the
        # x-seam if m1 == Lx, the y-seam if m2 == 1. A seam crossing
        # picks up the (−1)^ν sign for that direction.
        b_m1, crossed_x = m1 == Lx ? (1, true) : (m1 + 1, false)
        b_m2, crossed_y_xb = m2 == 1 ? (Ly, true) : (m2 - 1, false)
        s_x = (crossed_x && νx == 1) ? -1 : 1
        s_y_x = (crossed_y_xb && νy == 1) ? -1 : 1
        M[a, idx(b_m1, b_m2)] += s_x * s_y_x * m.Kx
        # y-bond: (m1, m2) → (m1, m2 − 1 mod Ly). Only the y-seam matters.
        c_m2, crossed_y = m2 == 1 ? (Ly, true) : (m2 - 1, false)
        s_y = (crossed_y && νy == 1) ? -1 : 1
        M[a, idx(m1, c_m2)] += s_y * m.Ky
    end
    return -sum(svdvals(M)) / (2 * Lx * Ly)
end

"""
    fetch(model::KitaevHoneycomb, ::Energy{:per_site}, bc::PBC; Lx, Ly) -> Float64

Per-site ground state energy on a `Lx × Ly` unit-cell torus (PBC in
both lattice directions) — enumerates all four topological flux
sectors and returns the minimum. Each sector corresponds to a choice
of fermion boundary conditions (`W_x, W_y ∈ {±1}`); Lieb's theorem
fixes plaquette fluxes at `+1`, so the spin-Hamiltonian ground state
is one of these four.

Bond connectivity matches `Lattice2D.build_lattice(Honeycomb, Lx, Ly;
boundary=PeriodicAxis())`. `bc.N` is ignored; pass `Lx`, `Ly` as kwargs.

For large `L` the four sectors converge to the same energy and
individual Bloch-sum terms dominate; for small `L` sector choice is
essential (e.g. `Lx = Ly = 2` gives distinct sector energies differing
by `~10%`).
"""
function fetch(model::KitaevHoneycomb, ::Energy{:per_site}, ::PBC; Lx::Integer, Ly::Integer)
    Lx > 0 && Ly > 0 ||
        error("KitaevHoneycomb PBC: Lx, Ly must be positive (got Lx=$Lx, Ly=$Ly).")
    return minimum(
        _kitaev_pbc_sector_energy(model, Lx, Ly, νx, νy) for νx in 0:1, νy in 0:1
    )
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy — OBC finite (per site)
#
# The unified GS + thermal `Energy{:per_site}, OBC` fetch lives in
# `KitaevHoneycomb_thermal.jl` (sibling include).  At `beta = nothing`
# (default) it returns the T = 0 ground-state value
# `−Σ σ_k / (2·Lx·Ly)`, where `σ_k` are the singular values of the
# bipartite hopping matrix `M` defined below.  With a real `beta` it
# returns the matter-sector thermal value (see thermal file).  Defining
# the method there keeps the GS path bit-for-bit identical while avoiding
# a method-overwrite warning on package load.
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _obc_hopping_matrix(model, Lx, Ly) -> Matrix{Float64}

`Lx·Ly × Lx·Ly` bipartite hopping matrix for the flux-free-sector
Majorana problem with open boundary conditions in both lattice
directions. `M[a, b] = K_γ` whenever the A-sublattice site of unit
cell `a` connects to the B-sublattice site of unit cell `b` via a
bond of type γ (matching Lattice2D's `Honeycomb` topology
conventions: type_1 ↔ K_z within the same cell, type_2 ↔ K_x crossing
`(+a₁, −a₂)`, type_3 ↔ K_y crossing `(−a₂)`). Bonds that would leave
the `Lx × Ly` strip are dropped.

Row indexing: `a = (m1 − 1) Ly + m2` for `(m1, m2) ∈ 1..Lx × 1..Ly`.
"""
function _obc_hopping_matrix(m::KitaevHoneycomb, Lx::Integer, Ly::Integer)
    idx(m1, m2) = (m1 - 1) * Ly + m2
    N = Lx * Ly
    M = zeros(Float64, N, N)
    for m1 in 1:Lx, m2 in 1:Ly
        a = idx(m1, m2)
        # z-bond: A(m1, m2) → B(m1, m2)
        M[a, idx(m1, m2)] += m.Kz
        # x-bond: A(m1, m2) → B(m1 + 1, m2 - 1)  iff in-range
        if (m1 + 1 ≤ Lx) && (m2 - 1 ≥ 1)
            M[a, idx(m1 + 1, m2 - 1)] += m.Kx
        end
        # y-bond: A(m1, m2) → B(m1, m2 - 1)      iff in-range
        if m2 - 1 ≥ 1
            M[a, idx(m1, m2 - 1)] += m.Ky
        end
    end
    return M
end

# ═══════════════════════════════════════════════════════════════════════════════
# MassGap — Infinite
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::KitaevHoneycomb, ::MassGap, ::Infinite) -> Float64

Single-Majorana gap in the thermodynamic limit.

    Δ = 2 · min_k |f(k)|.

In the A/B/C gapless phase (each |Kᵧ| ≤ sum of the other two), `f(k)`
has two linear (Dirac) zeros and `Δ = 0`. In the gapped A_z / A_x /
A_y phases (`|Kᵧ|` exceeds the sum of the other two), |f| is bounded
away from zero and `Δ = 2·( |K_γ_max| − |K_γ_other1| − |K_γ_other2| )`.
"""
function fetch(model::KitaevHoneycomb, ::MassGap, ::Infinite; kwargs...)
    ax, ay, az = abs(model.Kx), abs(model.Ky), abs(model.Kz)
    big = max(ax, ay, az)
    rest = ax + ay + az - big
    return 2 * max(big - rest, 0.0)
end
