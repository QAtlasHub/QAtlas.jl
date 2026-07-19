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
# References (DOIs given for direct lookup via doiget / DOI resolver):
#   I. Affleck, T. Kennedy, E. H. Lieb, H. Tasaki,
#     "Valence bond ground states in isotropic quantum antiferromagnets",
#     [AKLT1988](@cite).  doi:10.1007/BF01218021
#   T. Kennedy and H. Tasaki,
#     "Hidden Z₂ × Z₂ symmetry breaking in Haldane-gap antiferromagnets",
#     [KennedyTasaki1992](@cite) — string order parameter.
#     doi:10.1103/PhysRevB.45.304
#   D. P. Arovas, A. Auerbach, F. D. M. Haldane,
#     "Extended Heisenberg models of antiferromagnetism: Analogies to the
#     fractional quantum Hall effect", [Arovas1988](@cite) —
#     static structure factor of the AKLT chain.
#     doi:10.1103/PhysRevLett.60.531
#   A. García-Saez, V. Murg, and F. Verstraete,
#     "Spectral gap of the Affleck-Kennedy-Lieb-Tasaki Hamiltonian",
#     [GarciaSaez2013](@cite); arXiv:1308.3631 — DMRG gap.
#     doi:10.1103/PhysRevB.88.245118
#
# Finite-temperature context (no analytic reduction known at the AKLT
# point; the model is not Bethe-ansatz integrable, unlike the nearby
# Takhtajan-Babujian θ = π/4 and Lai-Sutherland θ = -π/4 BLBQ points):
#   S. Lou, T. Xiang, Z. Su,
#     "Thermodynamics of the Bilinear-Biquadratic Spin-One Heisenberg
#     Chain", [Lou2000](@cite) — TMRG numerical reference
#     for BLBQ thermodynamic curves (includes the AKLT point).
#     doi:10.1103/PhysRevLett.85.2380
#   A. Lohmann, H.-J. Schmidt, J. Richter,
#     "Tenth-order high-temperature expansion for the susceptibility and
#     the specific heat of spin-s Heisenberg models with arbitrary
#     exchange patterns", [Lohmann2014](@cite) — HTSE
#     methodology baseline (bilinear Heisenberg only; biquadratic
#     extension to AKLT is not in the literature; see QAtlas issue #506).
#     doi:10.1103/PhysRevB.89.014415
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
    J > 0 || throw(
        ArgumentError(
            "AKLT1D requires J > 0 (antiferromagnetic); got J = $J.  " *
            "Every registered analytic observable assumes J > 0 — see the " *
            "AKLT1D module docstring.",
        ),
    )
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
Phys. Rev. B **88**, 245118 (2013), arXiv:1308.3631, doi:10.1103/PhysRevB.88.245118).
No closed form is known; `reliability=:medium` in the registry.
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
#     doi:10.1007/BF01218021
#   D. P. Arovas, A. Auerbach, F. D. M. Haldane, Phys. Rev. Lett. 60,
#     531 (1988) — static structure factor of the AKLT chain.
#     doi:10.1103/PhysRevLett.60.531

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

# ═══════════════════════════════════════════════════════════════════════════════
# Finite-temperature thermodynamics — β = ∞ (T = 0) limits only
# ═══════════════════════════════════════════════════════════════════════════════
#
# The AKLT model is NOT Bethe-ansatz integrable (the integrable BLBQ
# points are Takhtajan-Babujian θ = π/4 and Lai-Sutherland θ = -π/4;
# AKLT sits at θ = arctan(1/3) ≈ 18.4° between them).  As a result, no
# closed-form analytic, transfer-matrix, or integral-equation reduction
# is known for the AKLT chain finite-T thermodynamic functions (free
# energy, entropy, specific heat, susceptibility) at arbitrary β.
#
# The published HTSE machinery (Lohmann-Schmidt-Richter 2014,
# doi:10.1103/PhysRevB.89.014415) handles only bilinear Heisenberg
# `J s_i · s_j`; the biquadratic `(s_i · s_j)²` term required by AKLT is
# not supported.  Lou-Xiang-Su 2000 (doi:10.1103/PhysRevLett.85.2380)
# give TMRG **numerical** thermodynamic curves for the BLBQ chain
# including the AKLT point, but no analytic form.
#
# Per the "src/ holds only analytic / TM / integral-equation reductions"
# policy, the only finite-T values shipped here are the **β → ∞ ground-
# state limits**, all of which are closed-form from the AKLT bond-
# projector decomposition (AKLT 1988):
#
#   H = 2J Σ_bonds P₂(i,j) − (2J/3) · N_bonds
#
# The VBS state annihilates every P₂(i,j) projector, so the GS energy is
# exactly  −(2J/3) · N_bonds  for any (OBC, PBC, Infinite) BC.
#
#   OBC (N ≥ 2, N_bonds = N − 1):
#     • f(∞) = −(2J/3) · (N − 1) / N        [exact]
#     • s(∞) = log 4 / N                    [4-fold edge-mode degeneracy:
#                                            two free spin-½ at the open
#                                            ends give a singlet ⊕ triplet
#                                            ground manifold; AKLT 1988]
#     • c(∞) = 0                            [Var(H) = 0 on any GS manifold]
#     • χ(∞) = ∞                            [Curie tail from edge spins;
#                                            ⟨(S^z_tot)²⟩_GS = 1/2 over
#                                            the 4-state manifold, so
#                                            χ ≈ β/(2N) → ∞.  Throws
#                                            DomainError; the divergence
#                                            coefficient 1/(2N) is exact
#                                            but not registered as a
#                                            separate observable here.]
#
#   PBC (N ≥ 2, N_bonds = N):  VBS is the unique frustration-free GS.
#     • f(∞) = −2J/3                        [exact, N-independent]
#     • s(∞) = 0, c(∞) = 0, χ(∞) = 0        [Haldane gap suppression]
#
#   Infinite:  same bulk as PBC at N → ∞.
#     • f(∞) = −2J/3, s(∞) = 0, c(∞) = 0, χ(∞) = 0
#
# Finite β > 0 (β ≠ ∞) throws DomainError everywhere.  Future work to
# extend this with a biquadratic-aware HTSE truncation (or a Padé-
# accelerated form) is tracked in QAtlas issue #506.

# ── OBC: β = ∞ limits ─────────────────────────────────────────────────────────

const _AKLT_NO_FINITE_BETA_MSG =
    "AKLT1D finite-β thermodynamics has no known " *
    "analytic / TM / integral-equation reduction (the model is not Bethe-ansatz " *
    "integrable).  Only β = Inf is supported here.  See AKLT1D module docstring " *
    "and QAtlas issue #506 for the planned HTSE extension."

"""
    fetch(model::AKLT1D, ::FreeEnergy, bc::OBC; beta::Real) -> Float64

Per-site Helmholtz free energy at `β = ∞` for the AKLT OBC chain on `N`
sites:

    f_OBC(∞) = −(2J/3) · (N − 1) / N

Exact closed form from the bond-projector decomposition (AKLT 1988,
doi:10.1007/BF01218021): the VBS ground state annihilates every bond
projector, so `E_GS_OBC = −(2J/3)(N − 1)` for any of the 4 degenerate
ground states.  Only `beta = Inf` is supported; finite-β throws
`DomainError`.
"""
function fetch(model::AKLT1D, ::FreeEnergy, bc::OBC; beta::Real, kwargs...)
    bc.N ≥ 2 || throw(ArgumentError("AKLT1D OBC FreeEnergy needs N ≥ 2 (got N = $(bc.N))"))
    isinf(beta) && beta > 0 || throw(DomainError(beta, _AKLT_NO_FINITE_BETA_MSG))
    return -(2.0 / 3.0) * model.J * (bc.N - 1) / bc.N
end

"""
    fetch(model::AKLT1D, ::ThermalEntropy, bc::OBC; beta::Real) -> Float64

Per-site Gibbs entropy at `β = ∞` for the AKLT OBC chain:

    s_OBC(∞) = log(4) / N

The OBC ground state is 4-fold degenerate (singlet ⊕ triplet of the two
free spin-½ edge modes; AKLT 1988, doi:10.1007/BF01218021), so the
β → ∞ Gibbs entropy of the equal-weight ground-manifold mixed state is
`log 4` total, hence `log 4 / N` per site.
"""
function fetch(::AKLT1D, ::ThermalEntropy, bc::OBC; beta::Real, kwargs...)
    bc.N ≥ 2 ||
        throw(ArgumentError("AKLT1D OBC ThermalEntropy needs N ≥ 2 (got N = $(bc.N))"))
    isinf(beta) && beta > 0 || throw(DomainError(beta, _AKLT_NO_FINITE_BETA_MSG))
    return log(4.0) / bc.N
end

"""
    fetch(model::AKLT1D, ::SpecificHeat, bc::OBC; beta::Real) -> Float64

Per-site heat capacity at `β = ∞` for the AKLT OBC chain:

    c_OBC(∞) = 0

The 4-fold degenerate ground manifold has zero energy variance (all four
states share the same energy `−(2J/3)(N − 1)`), so the variance form
`c = β² Var(H) / N` evaluates to 0 identically at `β = ∞`.
"""
function fetch(::AKLT1D, ::SpecificHeat, bc::OBC; beta::Real, kwargs...)
    bc.N ≥ 2 ||
        throw(ArgumentError("AKLT1D OBC SpecificHeat needs N ≥ 2 (got N = $(bc.N))"))
    isinf(beta) && beta > 0 || throw(DomainError(beta, _AKLT_NO_FINITE_BETA_MSG))
    return 0.0
end

"""
    fetch(model::AKLT1D, ::SusceptibilityZZ, bc::OBC; beta::Real)

Always throws `DomainError`.  The AKLT OBC longitudinal susceptibility
diverges as `β → ∞` (Curie tail from the two free spin-½ edge modes):
averaging `(S^z_tot)²` over the 4-fold ground manifold gives `1/2`, so

    χ_OBC(β → ∞) ≈ β / (2N) → ∞.

The divergence coefficient `1/(2N)` is exact (AKLT 1988,
doi:10.1007/BF01218021) but is not registered as a separate observable
here.  Finite β has no analytic reduction.  Use the PBC or Infinite
fetches for finite gap-suppressed values, or HeisenbergS1.jl's
ED-based path for finite-N numerical references.
"""
function fetch(::AKLT1D, ::SusceptibilityZZ, ::OBC; beta::Real, kwargs...)
    return throw(
        DomainError(
            beta,
            "AKLT1D OBC SusceptibilityZZ diverges at β = Inf (edge-mode Curie tail; " *
            "coefficient 1/(2N) is exact but not a registered observable).  " *
            "Finite β has no analytic reduction.  See AKLT1D module docstring.",
        ),
    )
end

# ── PBC: β = ∞ limits ─────────────────────────────────────────────────────────

"""
    fetch(model::AKLT1D, ::FreeEnergy, bc::PBC; beta::Real) -> Float64

Per-site Helmholtz free energy at `β = ∞` for the AKLT PBC chain:

    f_PBC(∞) = −2J / 3

Exact, `N`-independent for any `N ≥ 2`.  PBC has `N` bonds and the VBS
state is the unique frustration-free ground state annihilating every
bond projector (AKLT 1988, doi:10.1007/BF01218021); the GS energy is
`E_GS = −(2J/3) · N`.
"""
function fetch(model::AKLT1D, ::FreeEnergy, bc::PBC; beta::Real, kwargs...)
    bc.N ≥ 2 || throw(ArgumentError("AKLT1D PBC FreeEnergy needs N ≥ 2 (got N = $(bc.N))"))
    isinf(beta) && beta > 0 || throw(DomainError(beta, _AKLT_NO_FINITE_BETA_MSG))
    return -(2.0 / 3.0) * model.J
end

"""
    fetch(model::AKLT1D, ::ThermalEntropy, bc::PBC; beta::Real) -> Float64

Per-site Gibbs entropy at `β = ∞` for the AKLT PBC chain: `0`.  The PBC
VBS ground state is unique and gapped (Haldane gap, García-Saez–Murg–
Verstraete 2013, doi:10.1103/PhysRevB.88.245118), so the β → ∞ Gibbs
state is the pure GS with zero entropy.
"""
function fetch(::AKLT1D, ::ThermalEntropy, bc::PBC; beta::Real, kwargs...)
    bc.N ≥ 2 ||
        throw(ArgumentError("AKLT1D PBC ThermalEntropy needs N ≥ 2 (got N = $(bc.N))"))
    isinf(beta) && beta > 0 || throw(DomainError(beta, _AKLT_NO_FINITE_BETA_MSG))
    return 0.0
end

"""
    fetch(model::AKLT1D, ::SpecificHeat, bc::PBC; beta::Real) -> Float64

Per-site heat capacity at `β = ∞` for the AKLT PBC chain: `0`.  Pure GS
implies zero energy variance.
"""
function fetch(::AKLT1D, ::SpecificHeat, bc::PBC; beta::Real, kwargs...)
    bc.N ≥ 2 ||
        throw(ArgumentError("AKLT1D PBC SpecificHeat needs N ≥ 2 (got N = $(bc.N))"))
    isinf(beta) && beta > 0 || throw(DomainError(beta, _AKLT_NO_FINITE_BETA_MSG))
    return 0.0
end

"""
    fetch(model::AKLT1D, ::SusceptibilityZZ, bc::PBC; beta::Real) -> Float64

Per-site uniform longitudinal susceptibility at `β = ∞` for the AKLT
PBC chain: `0`.  The unique VBS ground state has `⟨S^z_i⟩ = 0` and the
Haldane-gapped excitation spectrum gives exponential suppression
`χ ~ exp(−β Δ)` (Δ ≈ 0.350 J, García-Saez–Murg–Verstraete 2013,
doi:10.1103/PhysRevB.88.245118), so the strict β → ∞ limit is `0`.
"""
function fetch(::AKLT1D, ::SusceptibilityZZ, bc::PBC; beta::Real, kwargs...)
    bc.N ≥ 2 ||
        throw(ArgumentError("AKLT1D PBC SusceptibilityZZ needs N ≥ 2 (got N = $(bc.N))"))
    isinf(beta) && beta > 0 || throw(DomainError(beta, _AKLT_NO_FINITE_BETA_MSG))
    return 0.0
end

# ── Infinite: β = ∞ limits ────────────────────────────────────────────────────

"""
    fetch(model::AKLT1D, ::FreeEnergy, ::Infinite; beta::Real) -> Float64

Per-site Helmholtz free energy at `β = ∞` for the infinite AKLT chain:

    f(∞) = −2J / 3

Same closed form as `GroundStateEnergyDensity` at `Infinite()`; exposed
through `FreeEnergy` for callers that route generic thermodynamic
observables.
"""
function fetch(
    model::AKLT1D,
    ::FreeEnergy,
    ::Infinite;
    scheme::Symbol=:canonical,
    beta::Real,
    kwargs...,
)
    scheme === :canonical ||
        return _aklt_thermo_infinite_scheme(model, FreeEnergy(), Val(scheme); beta=beta)
    isinf(beta) && beta > 0 || throw(DomainError(beta, _AKLT_NO_FINITE_BETA_MSG))
    return -(2.0 / 3.0) * model.J
end

"""
    fetch(model::AKLT1D, ::ThermalEntropy, ::Infinite; beta::Real) -> Float64

Per-site Gibbs entropy at `β = ∞` for the infinite AKLT chain: `0`.
Bulk Haldane phase has a unique GS modulo finite edge-mode corrections
(which vanish per-site as `N → ∞`).
"""
function fetch(
    model::AKLT1D,
    ::ThermalEntropy,
    ::Infinite;
    scheme::Symbol=:canonical,
    beta::Real,
    kwargs...,
)
    scheme === :canonical ||
        return _aklt_thermo_infinite_scheme(model, ThermalEntropy(), Val(scheme); beta=beta)
    isinf(beta) && beta > 0 || throw(DomainError(beta, _AKLT_NO_FINITE_BETA_MSG))
    return 0.0
end

"""
    fetch(model::AKLT1D, ::SpecificHeat, ::Infinite; beta::Real) -> Float64

Per-site heat capacity at `β = ∞` for the infinite AKLT chain: `0`.
"""
function fetch(
    model::AKLT1D,
    ::SpecificHeat,
    ::Infinite;
    scheme::Symbol=:canonical,
    beta::Real,
    kwargs...,
)
    scheme === :canonical ||
        return _aklt_thermo_infinite_scheme(model, SpecificHeat(), Val(scheme); beta=beta)
    isinf(beta) && beta > 0 || throw(DomainError(beta, _AKLT_NO_FINITE_BETA_MSG))
    return 0.0
end

"""
    fetch(model::AKLT1D, ::SusceptibilityZZ, ::Infinite; beta::Real) -> Float64

Per-site uniform longitudinal susceptibility at `β = ∞` for the infinite
AKLT chain: `0`.  Haldane-gapped bulk → exponential suppression
`χ ~ exp(−β Δ)`, strict limit is `0`.
"""
function fetch(::AKLT1D, ::SusceptibilityZZ, ::Infinite; beta::Real, kwargs...)
    isinf(beta) && beta > 0 || throw(DomainError(beta, _AKLT_NO_FINITE_BETA_MSG))
    return 0.0
end
