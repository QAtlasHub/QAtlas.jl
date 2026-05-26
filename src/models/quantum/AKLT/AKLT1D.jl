# ─────────────────────────────────────────────────────────────────────────────
# AKLT1D — Spin-1 bilinear-biquadratic chain at the AKLT point
#
# Hamiltonian (J > 0 antiferromagnetic):
#
#   H = J Σᵢ [ Sᵢ · Sᵢ₊₁ + (1/3) (Sᵢ · Sᵢ₊₁)² ]
#
# where Sᵅ are 3×3 spin-1 generators (s = 1, s(s+1) = 2).  This is the
# Affleck–Kennedy–Lieb–Tasaki (1987/1988) point of the S=1 BLBQ chain;
# the Hamiltonian is, up to a constant, twice the projector onto the
# total-spin-2 subspace of each bond:
#
#   J [ Sᵢ·Sᵢ₊₁ + (1/3)(Sᵢ·Sᵢ₊₁)² ] = (2J) P₂(i, i+1) − (2J/3)
#
# so the unique frustration-free (Infinite, PBC) ground state — the
# Valence Bond Solid (VBS) — is the simultaneous null space of every
# bond projector, with energy density e₀ = −2J/3.  Under OBC the bulk
# remains in the VBS but two free spin-½ edge degrees of freedom give a
# 4-fold ground-state degeneracy.
#
# Closed-form values (AKLT 1988):
#   * GS energy density (Infinite):            e₀ = −2J/3
#   * Correlation length (Infinite):           ξ = 1/log 3 ≈ 0.910
#   * String order parameter (Kennedy–Tasaki): O_str = 4/9
#
# Numerical-exact (no closed form):
#   * Haldane gap (Infinite):  Δ ≈ 0.350 J  (García-Saez–Murg–Verstraete 2013, DMRG)
#
# References:
#   I. Affleck, T. Kennedy, E. H. Lieb, H. Tasaki,
#     "Valence bond ground states in isotropic quantum antiferromagnets",
#     Commun. Math. Phys. 115, 477 (1988).
#   T. Kennedy and H. Tasaki,
#     "Hidden Z₂ × Z₂ symmetry breaking in Haldane-gap antiferromagnets",
#     Phys. Rev. B 45, 304 (1992) — string order parameter.
#   A. García-Saez, V. Murg, and F. Verstraete,
#     "Spectral gap of the Affleck-Kennedy-Lieb-Tasaki Hamiltonian",
#     Phys. Rev. B 88, 245118 (2013); arXiv:1308.3631 — DMRG gap.
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Spin S (this file)
#   Observable:  Spin S         (QAtlas-wide spin convention; see docs/src/conventions.md)

using LinearAlgebra: I, Hermitian, eigvals, kron

"""
    AKLT1D(; J::Real = 1.0) <: AbstractQAtlasModel

Spin-1 bilinear-biquadratic chain at the AKLT point,

    H = J Σᵢ [ Sᵢ · Sᵢ₊₁ + (1/3) (Sᵢ · Sᵢ₊₁)² ],

with `J > 0` antiferromagnetic.  Same local Hilbert space as
[`S1Heisenberg1D`](@ref) but with the biquadratic coefficient tuned to
the special AKLT value `1/3` where the ground state is an exact
Valence Bond Solid (VBS).  The infinite-system observables exposed
through `fetch` are closed-form (energy density, correlation length,
string order parameter); the Haldane gap is the García-Saez–Murg–Verstraete (2013)
DMRG numerical-exact value.

The constructor requires `J > 0`: the AKLT bond-projector decomposition
`H = 2J Σ P₂ − (2J/3) N_bonds` and every analytic value registered in
`AKLT1D_registry.jl` (GS energy density, correlation length, string
order parameter, mass gap, β=∞ thermodynamic limits) assumes the
antiferromagnetic sign.  Passing `J ≤ 0` would silently return
sign-flipped or physically meaningless values, so the constructor
throws `ArgumentError` at construction time instead.
"""
struct AKLT1D <: AbstractQAtlasModel
    J::Float64
end
function AKLT1D(; J::Real=1.0)
    J > 0 || throw(ArgumentError(
        "AKLT1D requires J > 0 (antiferromagnetic); got J = $J.  " *
        "Every registered analytic observable assumes J > 0 — see the " *
        "AKLT1D module docstring."
    ))
    return AKLT1D(Float64(J))
end

# Reuse the spin-1 primitives + bond cap defined alongside `S1Heisenberg1D`
# so the AKLT chain inherits the same `_MAX_ED_SITES_S1 = 8` budget for
# 3^N dense ED.

"""
    _aklt_hamiltonian_matrix(model::AKLT1D, N::Int, ::OBC) -> Matrix{ComplexF64}

Assemble the dense `3^N × 3^N` OBC Hamiltonian

    H = J Σᵢ [ Sᵢ · Sᵢ₊₁ + (1/3) (Sᵢ · Sᵢ₊₁)² ]

via explicit tensor products built from the spin-1 primitives in
`HeisenbergS1.jl`.  Capped by `_MAX_ED_SITES_S1`.
"""
function _aklt_hamiltonian_matrix(model::AKLT1D, N::Int, ::OBC)
    N ≥ 2 || throw(ArgumentError("AKLT1D OBC chain needs N ≥ 2 (got N = $N)"))
    N ≤ _MAX_ED_SITES_S1 || throw(
        ArgumentError("spin-1 dense ED is capped at N ≤ $(_MAX_ED_SITES_S1) (got N = $N)"),
    )
    J = model.J
    D = 3^N
    # 9×9 single-bond block: S₁·S₂ + (1/3)(S₁·S₂)²
    SdotS = kron(_S1_x, _S1_x) + kron(_S1_y, _S1_y) + kron(_S1_z, _S1_z)
    bond = J * (SdotS + (1.0 / 3.0) * (SdotS * SdotS))
    H = zeros(ComplexF64, D, D)
    for i in 1:(N - 1)
        d_left = 3^(i - 1)
        d_right = 3^(N - i - 1)
        H .+= kron(
            Matrix{ComplexF64}(I, d_left, d_left),
            bond,
            Matrix{ComplexF64}(I, d_right, d_right),
        )
    end
    return H
end

# ═══════════════════════════════════════════════════════════════════════════════
# Infinite-system closed-form values
# ═══════════════════════════════════════════════════════════════════════════════

native_energy_granularity(::AKLT1D, ::Infinite) = :per_site

"""
    fetch(model::AKLT1D, ::GroundStateEnergyDensity, ::Infinite) -> Float64

Closed-form ground-state energy density of the AKLT chain in the
thermodynamic limit:

    e₀ = −(2/3) J

Derived analytically from the projector form `H = 2J Σ P₂(i,i+1) −
(2J/3) (N − 1)` (AKLT 1988): the VBS state is the exact null space of
every bond `P₂` projector, so per-bond energy is `−2J/3` and per-site
energy density (one bond per site in the bulk) is `−2J/3`.
"""
function fetch(model::AKLT1D, ::GroundStateEnergyDensity, ::Infinite; kwargs...)
    return -(2.0 / 3.0) * model.J
end

"""
    fetch(model::AKLT1D, ::Energy{:per_site}, ::Infinite) -> Float64

Per-site ground-state energy `e₀ = −2J/3` of the infinite AKLT chain.
Numerically identical to
`fetch(::AKLT1D, ::GroundStateEnergyDensity, ::Infinite)`; provided so
the BC-explicit `Energy(:per_site)` API resolves through the same
analytic path.
"""
function fetch(model::AKLT1D, ::Energy{:per_site}, ::Infinite; kwargs...)
    return -(2.0 / 3.0) * model.J
end

"""
    fetch(model::AKLT1D, ::CorrelationLength, ::Infinite) -> Float64

Closed-form bulk correlation length of the AKLT chain,

    ξ = 1 / log 3 ≈ 0.91024

(AKLT 1988).  Connected `⟨S^z_0 S^z_r⟩` decays as `(−1)^r (4/3) · 3^{−|r|}`
in the VBS state, giving `ξ = 1/log 3` independent of `J`.
"""
function fetch(model::AKLT1D, ::CorrelationLength, ::Infinite; kwargs...)
    return 1.0 / log(3.0)
end

"""
    fetch(model::AKLT1D, ::MassGap, ::Infinite) -> Float64

Haldane gap of the AKLT chain in the thermodynamic limit,

    Δ ≈ 0.350 J

(numerical-exact DMRG value; A. García-Saez, V. Murg, and F. Verstraete,
Phys. Rev. B **88**, 245118 (2013), arXiv:1308.3631).  No closed form is known; `reliability=:medium`
in the registry.
"""
function fetch(model::AKLT1D, ::MassGap, ::Infinite; kwargs...)
    return 0.350 * model.J
end

"""
    fetch(model::AKLT1D, ::StringOrderParameter, ::Infinite) -> Float64

Kennedy–Tasaki non-local (string) order parameter of the AKLT chain,

    O_str = 4/9

(closed form; AKLT 1988, Kennedy–Tasaki 1992).  This is the
infinite-distance limit of

    O_str(r) = −⟨ S^z_i exp[iπ Σ_{i<k<j} S^z_k] S^z_j ⟩

evaluated in the VBS ground state, and detects the hidden
`Z₂ × Z₂` symmetry breaking that defines the Haldane phase.  Independent
of `J`.
"""
function fetch(model::AKLT1D, ::StringOrderParameter, ::Infinite; kwargs...)
    return 4.0 / 9.0
end

# ═══════════════════════════════════════════════════════════════════════════════
# OBC dense-ED — full spectrum on N ≤ 8 sites
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::AKLT1D, ::ExactSpectrum, ::OBC; N::Int) -> Vector{Float64}

Sorted full eigenvalue spectrum of the OBC AKLT chain on `N` sites,
computed by dense ED on the `3^N`-dimensional Hilbert space.  Capped by
`_MAX_ED_SITES_S1` (so `N ≤ 8`, `3^8 = 6561`).

Under OBC the AKLT ground state is **4-fold degenerate** (S=½ edge
modes at each end, total spin `S_tot ∈ {0, 1}` from singlet ⊕ triplet
of the two edge ½-spins), and dense ED on `N ≤ 8` already exhibits
this degeneracy with `Δ_low_4 = E₃ − E₀` of order `10^{-13}` (only
floating-point noise).
"""
function fetch(model::AKLT1D, ::ExactSpectrum, bc::OBC; N::Int=bc.N, kwargs...)
    N > 0 || throw(ArgumentError("AKLT1D ExactSpectrum: N must be positive (got $N)"))
    H = _aklt_hamiltonian_matrix(model, N, OBC(N))
    return sort(real.(eigvals(Hermitian(H))))
end

# ═══════════════════════════════════════════════════════════════════════════════
# VBS ground-state spin correlations — exact closed form (AKLT 1988)
# ═══════════════════════════════════════════════════════════════════════════════
#
# The AKLT valence-bond-solid ground state has an exactly known,
# J-independent (the GS wavefunction is the same for every J > 0)
# equal-time two-point function:
#
#     ⟨Sᶻ₀ Sᶻ_r⟩ = (−1)^r · (4/3) · 3^{−|r|}     (r ≠ 0),
#     ⟨(Sᶻ)²⟩    = 2/3                            (r = 0, S = 1 on-site),
#
# equivalently ⟨S⃗₀·S⃗_r⟩ = (−1)^r · 4 · 3^{−|r|} by isotropy.  The decay
# ratio 3^{−1} per step is exactly the origin of ξ = 1/log 3.  ⟨Sᶻ⟩ = 0
# in the VBS, so this equal-time correlator is already the *connected*
# one.  Its lattice Fourier transform is the static structure factor
#
#     S_zz(q) = Σ_r e^{iqr} ⟨Sᶻ₀ Sᶻ_r⟩ = 2(1 − cos q) / (5 + 3 cos q)
#
# (Arovas–Auerbach–Haldane 1988): S_zz(0) = 0 (total-Sᶻ conservation),
# antiferromagnetic peak S_zz(π) = 2.
#
# References:
#   I. Affleck, T. Kennedy, E. H. Lieb, H. Tasaki, Commun. Math. Phys.
#     115, 477 (1988) — exact VBS two-point function.
#   D. P. Arovas, A. Auerbach, F. D. M. Haldane, Phys. Rev. Lett. 60,
#     531 (1988) — static structure factor of the AKLT chain.

"""
    fetch(model::AKLT1D, ::ZZCorrelation{:static}, ::Infinite; r::Integer) -> Float64

Exact equal-time spin-z two-point function of the AKLT VBS ground
state at separation `r`:

    ⟨Sᶻ₀ Sᶻ_r⟩ = (−1)^r · (4/3) · 3^{−|r|}   (r ≠ 0),
    ⟨(Sᶻ)²⟩    = 2/3                          (r = 0).

Closed form (AKLT 1988); `J`-independent (the VBS ground state does not
depend on `J > 0`).  Exponential decay with rate `log 3`, consistent
with [`fetch(::AKLT1D, ::CorrelationLength, ::Infinite)`](@ref)
`ξ = 1/log 3`.  Since `⟨Sᶻ⟩ = 0` in the VBS this equal-time value is
already the connected correlation.
"""
function fetch(::AKLT1D, ::ZZCorrelation{:static}, ::Infinite; r::Integer, kwargs...)
    r == 0 && return 2.0 / 3.0
    a = abs(r)
    return (iseven(a) ? 1.0 : -1.0) * (4.0 / 3.0) * 3.0^(-a)
end

"""
    fetch(model::AKLT1D, ::ZZStructureFactor, ::Infinite; q::Real) -> Float64

Exact static (equal-time) spin-z structure factor of the AKLT chain,

    S_zz(q) = Σ_r e^{iqr} ⟨Sᶻ₀ Sᶻ_r⟩ = 2 (1 − cos q) / (5 + 3 cos q)

(Arovas–Auerbach–Haldane 1988).  `S_zz(0) = 0` by total-`Sᶻ`
conservation; antiferromagnetic peak `S_zz(π) = 2`.  `J`-independent.
This is the lattice-Fourier sum of the closed-form
[`ZZCorrelation`](@ref) `(−1)^r (4/3) 3^{−|r|}`.
"""
function fetch(::AKLT1D, ::ZZStructureFactor, ::Infinite; q::Real, kwargs...)
    c = cos(q)
    return 2.0 * (1.0 - c) / (5.0 + 3.0 * c)
end
