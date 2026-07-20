# core/quantities.jl — concrete quantity struct library.
#
# Every physical observable that `fetch` can return is represented by a
# concrete subtype of `AbstractQuantity`.  Compared with the older
# `Quantity{:foo}` phantom-type pattern this gains:
#
#   * static dispatch (compiler sees the type, not a Symbol)
#   * compile-time argument checks (e.g. `RenyiEntropy(-1)` is rejected
#     by the inner constructor)
#   * unambiguous names — axis-indexed for tensor quantities, entropy
#     flavour spelled out, real-space / Fourier-space correlators kept
#     as separate types
#
# The legacy symbol dispatch still works through the `Quantity{S}()` shim
# in `core/type.jl` + canonicalize aliases in `core/alias.jl`.  That path
# is routed through `_symbol_to_quantity` in `deprecate/` (Milestone 1).

# ─── Quantity taxonomy: abstract family layer (#690) ────────────────────
#
# Families that are physically the same object up to a component / index
# (χ_xx/χ_yy/χ_zz, m_x/m_y/m_z, …) share an intermediate abstract supertype
# between the concrete leaf and `AbstractQuantity`.  Purely additive: every
# leaf keeps its name and `<: AbstractQuantity` still holds transitively, so
# all existing dispatch is unchanged.  The layer is what lets
#   * identities be declared once against a family instead of per leaf
#     (`@identity … family=AbstractSusceptibility`, core/identity.jl), and
#   * the atlas graph group a family as one node cluster instead of N
#     disconnected leaves.
# The component that the leaf's *name* encodes is recovered by the
# [`component`](@ref) trait.

# The 8 abstract quantity families now live in AbstractQAtlas, imported in src/QAtlas.jl (#734).

"""
    component(q) -> Union{Symbol,Nothing}
    component(::Type{<:AbstractQuantity}) -> Union{Symbol,Nothing}

The component / index that a family leaf's *type name* encodes: the spin
axis of a magnetization (`:x`/`:y`/`:z`), the diagonal axis pair of a
susceptibility / correlator / structure factor (`:xx`/`:yy`/`:zz`), or the
excitation channel of a gap (`:mass`/`:charge`/`:spin`).  `nothing` for
quantities that carry no component (the default), including the
site-resolved `…Local` magnetizations (whose extra site argument makes them
a different fetch shape).

Identities that hold per-component (e.g. the static FDT
`χ_αα = β·Var(M_α)/N`, or the SU(2) isotropy `χ_xx = χ_yy = χ_zz`) pair
family members by matching `component` — see `core/identity.jl`.
"""
component(q::AbstractQuantity) = component(typeof(q))
component(::Type{<:AbstractQuantity}) = nothing

# ─── Scalar thermodynamics ──────────────────────────────────────────────

# `Energy{G}`, its outer constructors, and the `native_energy_granularity` trait
# stub now live in AbstractQAtlas (#734) -- imported in src/QAtlas.jl.  The
# BC-granularity `fetch` routing methods below dispatch on the imported symbols.

"""
    FidelitySusceptibility() <: AbstractQuantity

Fidelity susceptibility `χ_F(λ) = −∂²⟨ψ(λ)|ψ(λ + δλ)⟩/∂δλ²`.
"""
struct FidelitySusceptibility <: AbstractQuantity end

# `PartitionFunction`, `CriticalTemperature`, `SpontaneousMagnetization`
# are currently defined in src/models/classical/IsingSquare/IsingSquare.jl
# as bare `struct X end` tags.  They will be migrated to subtype
# `AbstractQuantity` in the IsingSquare refactor commit (M1.7).

# ─── Entropies (explicit variants; see user-requested naming) ──────────

"""
    VonNeumannEntropy{M}() <: AbstractQuantity
    VonNeumannEntropy()                   # M = :equilibrium (default)
    VonNeumannEntropy(:equilibrium)       # explicit equilibrium S(ℓ)
    VonNeumannEntropy(:quench)            # post-quench S(ℓ, t)

Von Neumann entanglement entropy of a reduced density matrix,
`S_vN = −Tr ρ_A log ρ_A`.  The mode parameter `M::Symbol` is a
phantom type that splits the dispatch into:

- `:equilibrium` — equilibrium / thermal value
  `S_vN = -Tr ρ_A log ρ_A` of the GS or thermal reduced density
  matrix on the first `ℓ` sites (kwarg `ℓ`).  This is the original
  meaning; the no-argument constructor `VonNeumannEntropy()` keeps
  back-compatibility by routing here.

- `:quench` — time-evolved entanglement entropy
  `S_vN(ℓ, t) = -Tr ρ_A(t) log ρ_A(t)` after a sudden quench from
  the ground state of an `initial::AbstractQAtlasModel` (`H_0`) to
  the post-quench Hamiltonian (the `model` argument to `fetch`).
  Calabrese–Cardy quasi-particle picture (J. Stat. Mech. P04010
  (2005)): linear growth `S(ℓ, t) ≈ (c/3) v_E t` for `t ≪ ℓ/(2 v_E)`,
  saturating at `(c/3) log ℓ + const` for `t ≫ ℓ/(2 v_E)`.

See `docs/src/calc/tfim-quench-entanglement.md` for the
free-fermion derivation in the TFIM.
"""
struct VonNeumannEntropy{M} <: AbstractEntanglementMeasure
    function VonNeumannEntropy{M}() where {M}
        M isa Symbol || error("VonNeumannEntropy mode must be a Symbol, got $(typeof(M))")
        M in (:equilibrium, :quench) ||
            error("unknown VonNeumannEntropy mode :$M; expected :equilibrium or :quench")
        return new{M}()
    end
end
VonNeumannEntropy() = VonNeumannEntropy{:equilibrium}()
VonNeumannEntropy(m::Symbol) = VonNeumannEntropy{m}()

"""
    RenyiEntropy(α) <: AbstractQuantity

Rényi entropy of order `α`, `S_α = (1 − α)⁻¹ log Tr ρ_A^α`.

- `α = 1` recovers [`VonNeumannEntropy`](@ref) (implementations may
  dispatch accordingly).
- `α = 2` is the second Rényi entropy, frequently measured
  experimentally.
- `α > 0`, `α ≠ 1` are the supported generic cases.

The inner constructor rejects `α ≤ 0` and `α = 1` (use
`VonNeumannEntropy()` explicitly) — this is intentional, to force the
call site to be explicit about which entropy it wants.
"""
struct RenyiEntropy <: AbstractEntanglementMeasure
    α::Float64
    function RenyiEntropy(α::Real)
        α > 0 || throw(ArgumentError("RenyiEntropy: α must be positive; got $α"))
        α == 1 && throw(
            ArgumentError(
                "RenyiEntropy(1) is ambiguous; use VonNeumannEntropy() explicitly."
            ),
        )
        return new(Float64(α))
    end
end

"""
    ResidualEntropy() <: AbstractQuantity

Zero-temperature configurational (residual) entropy density.  Real-
valued, non-negative; non-zero in the presence of macroscopic
ground-state degeneracy (e.g. ice rule, frustrated Ising AFM,
Pauling-1935-style models).  Distinct from [`ThermalEntropy`](@ref):
`ThermalEntropy` is a finite-temperature thermodynamic quantity, while
`ResidualEntropy` is the lim_{T -> 0} S(T) / N residual term.
Zero-temperature ground-state entropy per site,

    S_residual / (N k_B) = lim_{T → 0⁺} S(T) / N,

i.e. the entropy density of the (possibly degenerate) ground-state
manifold.  Non-zero for frustrated classical models with extensive
ground-state degeneracy — e.g. the antiferromagnetic Ising model on
the triangular lattice (Wannier 1950, ≈ 0.3230659669) and on the
kagome lattice (Houtappel 1950).  Defined as a separate quantity
from [`ThermalEntropy`](@ref) to keep the zero-temperature limit
explicit at the dispatch level (avoiding β → ∞ extrapolations of a
finite-T fetch).
"""
struct ResidualEntropy <: AbstractThermalPotential end

# ─── Magnetizations (axis explicit) ─────────────────────────────────────

# The global bulk magnetizations adopt AbstractQAtlas's axis-parametric
# `Magnetization{A}` (#734).  The fused `MagnetizationX/Y/Z` names stay as
# deprecated aliases so existing fetch methods, registry rows, and downstream
# code keep working; new code should use `Magnetization(:x)` / `(:y)` / `(:z)`.

"""
    const MagnetizationX = Magnetization{:x}   # deprecated alias

Bulk-averaged `⟨σˣ⟩` in Pauli convention (= 2 ⟨Sˣ⟩ in spin-1/2 units).
For a spin-1/2 chain `H = -J ΣSᶻSᶻ - h ΣSˣ` this is the transverse
magnetization.  **Deprecated** — use `Magnetization(:x)`.
"""
const MagnetizationX = Magnetization{:x}

"""
    const MagnetizationY = Magnetization{:y}   # deprecated alias

Bulk-averaged `⟨σʸ⟩`.  **Deprecated** — use `Magnetization(:y)`.
"""
const MagnetizationY = Magnetization{:y}

"""
    const MagnetizationZ = Magnetization{:z}   # deprecated alias

Bulk-averaged `⟨σᶻ⟩`.  For Z₂-symmetric phases on an infinite system this is
the order parameter at low temperature.  **Deprecated** — use `Magnetization(:z)`.
"""
const MagnetizationZ = Magnetization{:z}

"""
    MagnetizationXLocal{M}() <: AbstractQuantity
    MagnetizationXLocal()                       # M = :equilibrium (default)
    MagnetizationXLocal(:equilibrium)           # explicit equilibrium ⟨σˣ_i⟩_β
    MagnetizationXLocal(:quench)                # post-quench ⟨σˣ_i⟩(t)

Site-resolved `⟨σˣ_i⟩` quantity.  The mode parameter `M::Symbol` is a
phantom type that splits the dispatch into:

- `:equilibrium` — site-resolved thermal expectation
  `[⟨σˣ_i⟩_β for i = 1:N]` (Vector{Float64}).  This is the original
  meaning; the no-argument constructor `MagnetizationXLocal()` keeps
  back-compatibility by routing here.

- `:quench` — time-evolved local transverse magnetisation
  `⟨σˣ_i⟩(t) = ⟨ψ_0|e^{iH_f t} σˣ_i e^{-iH_f t}|ψ_0⟩` after a sudden
  quench from the ground state of an `initial::AbstractQAtlasModel`
  (`H_0`) to the post-quench Hamiltonian (the `model` argument to
  `fetch`).  Returns a single `Float64` for one `(i, t)` pair.

See `docs/src/calc/tfim-sigma-x-quench.md` for the closed-form
derivation in the TFIM (Calabrese–Essler–Fagotti, J. Stat. Mech.
P07016 (2012); Barouch–McCoy–Dresden, PRA **2** (1970)).
"""
struct MagnetizationXLocal{M} <: AbstractMagnetization
    function MagnetizationXLocal{M}() where {M}
        M isa Symbol ||
            error("MagnetizationXLocal mode must be a Symbol, got \$(typeof(M))")
        M in (:equilibrium, :quench) ||
            error("unknown MagnetizationXLocal mode :\$M; expected :equilibrium or :quench")
        return new{M}()
    end
end
MagnetizationXLocal() = MagnetizationXLocal{:equilibrium}()
MagnetizationXLocal(m::Symbol) = MagnetizationXLocal{m}()

"""
    MagnetizationYLocal() <: AbstractQuantity

Site-resolved `⟨σʸ_i⟩` vector of length `N_bulk`.  Identically zero
for any real Hermitian Hamiltonian (parity / time-reversal); a model
that returns it explicitly does so as an exact baseline against
random-sample estimators that fluctuate around zero.
"""
struct MagnetizationYLocal <: AbstractMagnetization end

"""
    MagnetizationZLocal() <: AbstractQuantity

Site-resolved `⟨σᶻ_i⟩` vector of length `N_bulk`.
"""
struct MagnetizationZLocal <: AbstractMagnetization end

"""
    EnergyLocal() <: AbstractQuantity

Bond-resolved energy density vector, length `N_bulk − 1` for a bond
Hamiltonian `Σ_b h_b`.
"""
struct EnergyLocal <: AbstractQuantity end

# ─── Susceptibilities (axis pair) ────────────────────────────────────────

# The static axis-pair susceptibilities adopt AbstractQAtlas's index-parametric
# `Susceptibility{I}` (#734).  The fused `SusceptibilityXX/YY/ZZ` names stay as
# deprecated aliases so existing fetch methods, registry rows, and downstream
# code keep working; new code should use `Susceptibility(:x, :x)` etc.

"""
    const SusceptibilityXX = Susceptibility{(:x, :x)}   # deprecated alias

Static transverse susceptibility,
`χ_xx(β) = β · (⟨M_x²⟩ − ⟨M_x⟩²) / N`.  **Deprecated** — use `Susceptibility(:x, :x)`.
"""
const SusceptibilityXX = Susceptibility{(:x, :x)}

"""
    const SusceptibilityYY = Susceptibility{(:y, :y)}   # deprecated alias

Analogue for the y-axis.  **Deprecated** — use `Susceptibility(:y, :y)`.
"""
const SusceptibilityYY = Susceptibility{(:y, :y)}

"""
    const SusceptibilityZZ = Susceptibility{(:z, :z)}   # deprecated alias

Uniform longitudinal susceptibility,
`χ_zz(β) = β · (⟨M_z²⟩ − ⟨M_z⟩²) / N`.  **Deprecated** — use `Susceptibility(:z, :z)`.
"""
const SusceptibilityZZ = Susceptibility{(:z, :z)}

# ─── Real-space two-point correlators ───────────────────────────────────
#
# `XXCorrelation` / `YYCorrelation` / `ZZCorrelation` all carry a `mode`
# field so the same type dispatches static / dynamic / light-cone / …
# variants.  A model may implement only a subset of modes; `fetch`
# methods should error explicitly for unsupported modes.

"""
    ZZCorrelation{M}() <: AbstractQuantity
    ZZCorrelation(; mode::Symbol = :static)

Real-space 2-point correlator `⟨σᶻ_i σᶻ_j⟩`.  The mode `M::Symbol` is
a phantom type parameter so dispatch can specialise on it.

Supported `mode` values (by convention; individual models need only
implement the ones they support):

- `:static` — equal-time, thermal or zero-temperature value
- `:connected` — `⟨σᶻ_i σᶻ_j⟩ − ⟨σᶻ_i⟩⟨σᶻ_j⟩`
- `:dynamic` — retarded real-time correlator `⟨σᶻ_i(t) σᶻ_j(0)⟩`
- `:lightcone` — space-time spreading `⟨σᶻ_i(t) σᶻ_j(0)⟩` as a
  matrix over (site, time)

The companion type for Fourier-space structure factors is
[`ZZStructureFactor`](@ref), kept separate because it carries (q, ω)
arguments instead of (i, j, t).
"""
struct ZZCorrelation{M} <: AbstractTwoPointCorrelation end
ZZCorrelation(; mode::Symbol=:static) = ZZCorrelation{mode}()

"""
    XXCorrelation{M}() <: AbstractQuantity
    XXCorrelation(; mode::Symbol = :static)

Real-space 2-point `⟨σˣ_i σˣ_j⟩` correlator.  See
[`ZZCorrelation`](@ref) for the `mode` semantics.
"""
struct XXCorrelation{M} <: AbstractTwoPointCorrelation end
XXCorrelation(; mode::Symbol=:static) = XXCorrelation{mode}()

"""
    YYCorrelation{M}() <: AbstractQuantity
    YYCorrelation(; mode::Symbol = :static)

Real-space 2-point `⟨σʸ_i σʸ_j⟩` correlator.
"""
struct YYCorrelation{M} <: AbstractTwoPointCorrelation end
YYCorrelation(; mode::Symbol=:static) = YYCorrelation{mode}()

# ─── Fourier-space structure factors (q, ω) ────────────────────────────

"""
    ZZStructureFactor() <: AbstractQuantity

Fourier-space structure factor
`S_zz(q, ω) = ∫ dt e^{iωt} (1/N) Σ_{ij} e^{iq·(i-j)} ⟨σᶻ_i(t)σᶻ_j(0)⟩`
(or its static limit, depending on the model's fetch signature).

Kept as a separate type from [`ZZCorrelation`](@ref) because the
argument domain is (q, ω) instead of (i, j, t) and because existing
users already expect a dedicated `StructureFactor` quantity.
"""
struct ZZStructureFactor <: AbstractStructureFactor end

"""
    XXStructureFactor() <: AbstractQuantity

Fourier-space equivalent of [`XXCorrelation`](@ref).
"""
struct XXStructureFactor <: AbstractStructureFactor end

"""
    YYStructureFactor() <: AbstractQuantity

Fourier-space equivalent of [`YYCorrelation`](@ref).
"""
struct YYStructureFactor <: AbstractStructureFactor end

# ─── Universality / lattice spectra / advanced ─────────────────────────

"""
    ConformalWeights() <: AbstractQuantity

Primary scaling dimension `h` of a 2D rational CFT.  For Virasoro
[`MinimalModel`](@ref) this is the Kac-table entry `h_{r,s}`; for
[`WZWSU2`](@ref) it is the SU(2)-spin label `h_j = j(j+1)/(k+2)`.

Concrete model fetch methods take additional keyword arguments
identifying the primary (`r`, `s` for `MinimalModel`; `j` for
`WZWSU2`) and return an exact `Rational{Int}`.
"""
struct ConformalWeights <: AbstractQuantity end

"""
    PrimaryFields() <: AbstractQuantity

Full list of primary fields of a 2D rational CFT.  For
[`MinimalModel`](@ref) the result is a `Vector{NamedTuple{(:r, :s, :h)}}`
of length `(p - 1)(p_prime - 1) / 2`, with one entry per Kac-symmetry
orbit.

Future CFT classes may return different NamedTuple schemas (e.g.
`(j, h)` for WZW). The return type is therefore a
`Vector{<:NamedTuple}` whose schema depends on the model.
"""
struct PrimaryFields <: AbstractQuantity end

"""
    FractalDimension() <: AbstractQuantity

Hausdorff dimension `d_H` of the random geometric set associated with
a model — e.g. the SLE_κ curve's `d_H(κ) = min(2, 1 + κ/8)`
(Beffara 2008).  Real-valued, dimensionless, capped at the ambient
space dimension.
"""
struct FractalDimension <: AbstractQuantity end

"""
    StringOrderParameter() <: AbstractQuantity

Kennedy-Tasaki non-local (string) order parameter

    O_str = lim_{|i-j| -> infty} -<S^z_i exp[i pi sum_{i<k<j} S^z_k] S^z_j>

for S=1 chains.  Detects the hidden Z_2 x Z_2 symmetry breaking that
defines the Haldane phase (T. Kennedy and H. Tasaki, Phys. Rev. B **45**,
304 (1992)).  At the AKLT point the closed-form value is O_str = 4/9
(AKLT 1988), making it the canonical analytic test bed for any
implementation that aims to detect topologically non-trivial gapped
phases of integer-spin chains.
"""
struct StringOrderParameter <: AbstractQuantity end

"""
    LuttingerParameter() <: AbstractQuantity

Luttinger liquid parameter `K`.  Meaningful for critical 1D models
with U(1) symmetry (e.g. XXZ in the critical regime `|Δ| < 1`).
"""
struct LuttingerParameter <: AbstractQuantity end

"""
    FermiVelocity() <: AbstractQuantity

Fermi velocity `v_F = ∂ε/∂k |_{k_F}`.  Meaningful for non-interacting
/ mean-field fermionic band structures (tight-binding lattices,
Bogoliubov-de Gennes diagonalisations).  In QAtlas this is the type
returned by models like `Honeycomb` (at the Dirac cones), the
other tight-binding lattices, and the TFIM Majorana mode at the
critical field.
"""
struct FermiVelocity <: AbstractVelocity end

"""
    LuttingerVelocity() <: AbstractQuantity

Luttinger-liquid / bosonisation velocity `u` (a.k.a. `v_{LL}`) of the
low-energy linear-dispersion mode in a 1D critical interacting system.
Used by models like [`XXZ1D`](@ref) in the Luttinger regime
`|Δ| < 1`, the Heisenberg chain at the SU(2) point, and any other
bosonised 1D critical theory.

For a free-fermion model this coincides with [`FermiVelocity`](@ref);
for interacting systems `u` includes the Luttinger renormalisation.
"""
struct LuttingerVelocity <: AbstractVelocity end

"""
    const SpinWaveVelocity = LuttingerVelocity

Spin-chain community alias for [`LuttingerVelocity`](@ref).  The "spin
wave velocity" (e.g. in the Haldane / Affleck literature on the AFM
Heisenberg chain) is the same quantity as the Luttinger velocity once
bosonised; both dispatch through the same fetch method via the type
identity.
"""
const SpinWaveVelocity = LuttingerVelocity

"""
    LiebRobinsonVelocity() <: AbstractQuantity

Lieb-Robinson velocity `v_LR` setting the linear light cone for
information propagation in a local lattice quantum system: for any
local operators `A_x`, `B_y` separated by `|x - y|`,

    || [A_x(t), B_y(0)] || <= C * exp(-mu * (|x - y| - v_LR * t)).

For free-fermion-mappable spin chains (TFIM, XY/XYh1D, the XX limit
of XXZ) the bound is saturated and `v_LR` equals twice the maximum
single-particle group velocity. Reference: Lieb-Robinson, *Commun.
Math. Phys.* **28**, 251 (1972); Hastings-Koma, *Commun. Math. Phys.*
**265**, 781 (2006). Tracking: issue #579 inequality framework.
"""
struct LiebRobinsonVelocity <: AbstractVelocity end

"""
    EntanglementGrowthSlope() <: AbstractQuantity

Linear-growth slope of the half-system entanglement entropy after a
global quench from a thermal-like initial state. Calabrese-Cardy 2005
predicts that, for `t < L / (2 v)`,

    dS_A / dt = (π c v) / (3 β_eff),

where `c` is the central charge of the critical post-quench
Hamiltonian, `v` is the Lieb-Robinson velocity of correlation
spreading, and `β_eff` is the effective inverse temperature of the
generalised-Gibbs steady state set by the initial state. This struct
is the type tag; concrete dispatches live at the universality layer
and on model files. Reference: Calabrese-Cardy 2005, *J. Stat. Mech.*
P04010. Tracking: #580 quench-dynamics phase.
"""
struct EntanglementGrowthSlope <: AbstractQuantity end

"""
    EntanglementSaturationDensity() <: AbstractQuantity

Per-unit-length saturation value of post-quench entanglement entropy
in the Calabrese-Cardy 2005 picture: in the long-time regime
`t > L / (2 v)`, the half-system entropy saturates at

    S_A(infty) / L = π c / (6 beta_eff),

where `c` is the central charge of the post-quench critical
Hamiltonian and `beta_eff` is the effective inverse temperature of the
generalised-Gibbs steady state. Universal in (c, beta_eff). Partner to
[`EntanglementGrowthSlope`](@ref) (which gives the dS/dt of the
linear regime). Reference: Calabrese-Cardy J. Stat. Mech. P04010
(2005). Tracking: #580.
"""
struct EntanglementSaturationDensity <: AbstractQuantity end

"""
    ThermalEnergyDensity <: AbstractQuantity

Leading low-temperature thermal energy density of a (1+1)D
conformal field theory above its ground state,

    e(T) - e_0 = pi c T^2 / 6 = pi c / (6 beta^2),

where  is the central charge and  (Affleck 1986;
Bloete-Cardy-Nightingale 1986). This is the universal counterpart of
[`ConformalCasimirEnergy`](@ref): the same `c` prefactor that
controls the Casimir term in finite size also fixes the leading
thermal-excitation density in finite temperature, via modular
invariance.
"""
struct ThermalEnergyDensity <: AbstractQuantity end

"""
    CFTThermalEntropyDensity <: AbstractQuantity

Leading low-temperature thermal entropy density (entropy per unit
length) of a (1+1)D conformal field theory,

    s(T) = pi c T / 3 = pi c / (3 beta),

with  the central charge. This is the temperature derivative of
the universal CFT free energy density  (Bloete-Cardy-
Nightingale 1986), and the operational complement of
[`ThermalEnergyDensity`](@ref).
"""
struct CFTThermalEntropyDensity <: AbstractQuantity end

"""
    WignerSemicircleMoment <: AbstractQuantity

Moments of the Wigner semicircle distribution

    rho(x) = (1 / (2 pi)) sqrt(4 - x^2),   x in [-2, 2],

the universal large-N eigenvalue density of Gaussian random matrix
ensembles (GOE / GUE / GSE) under Wigner-Mehta normalisation.

The even moments are Catalan numbers,

    m_{2k} = C_k = (2k)! / (k! (k+1)!),

and the odd moments vanish by symmetry. These are the universal large-N
free-probability moments underlying RMT spectral statistics; they also
count rooted plane trees / non-crossing pair partitions.

Reference: E. P. Wigner, *Ann. Math.* **62**, 548 (1955); M. L. Mehta
*Random Matrices* (1991).
"""
struct WignerSemicircleMoment <: AbstractQuantity end

"""
    CardyEntropy() <: AbstractQuantity

Asymptotic high-energy entropy (log of the density of states) of a
1+1D CFT, given by the Cardy 1986 formula

    S_Cardy(E) = 2 π sqrt(c E / 6),

where `c` is the central charge of the CFT and `E` is the excitation
energy (in units where the cylinder circumference is 1). This counts
the number of CFT states at fixed energy `E` and underlies, e.g., the
Cardy-Verlinde / black-hole-entropy correspondences. Tracking: #580.
"""
struct CardyEntropy <: AbstractQuantity end

"""
    ConformalCasimirEnergy() <: AbstractQuantity

Universal Casimir (ground-state) energy of a 1+1D CFT on a cylinder
of circumference `L` (PBC). Cardy 1986 / Blote-Cardy-Nightingale 1986
/ Affleck 1986 showed it is determined entirely by the central charge:

    E_0(L) = -π c / (6 L).

This is the strict thermodynamic-limit subtraction `lim_{L->∞}
(E_GS(L) - L * e_∞) * L` that extracts the universal finite-size
correction.  Sign convention follows the original PRL: `E_0 < 0` for
unitary CFTs with `c > 0`. Tracking: #580.
"""
struct ConformalCasimirEnergy <: AbstractQuantity end

"""
    LogarithmicNegativity() <: AbstractQuantity

Logarithmic negativity `E = log Tr |ρ^{T_B}|` measuring mixed-state
entanglement between two subsystems. For two adjacent intervals on
an infinite 1+1D-CFT chain at T = 0, the universal closed form
(Calabrese-Cardy-Tonni 2012) is

    E(ℓ_A, ℓ_B) = (c/4) log[ℓ_A · ℓ_B / (ℓ_A + ℓ_B)],

i.e., the same geometric-mean log of the mutual-information universal
formula with the prefactor c/3 replaced by c/4. Tracking: #580.
"""
struct LogarithmicNegativity <: AbstractEntanglementMeasure end

"""
    BoundaryEntropy() <: AbstractQuantity

Affleck-Ludwig universal (non-integer) boundary entropy `log g`
of a conformal boundary state in a 1+1D rational CFT, given by

    g_a = S_{0a} / sqrt(S_{00})

for the Cardy boundary state |a⟩ corresponding to primary `a`, where
`S_{ab}` is the modular S-matrix. The quantity `log g` is non-negative
under unitary RG and decreases monotonically (g-theorem). The
universal "ground-state degeneracy" interpretation goes back to
Affleck-Ludwig 1991. Tracking: #580.
"""
struct BoundaryEntropy <: AbstractQuantity end

"""
    PageEntropy() <: AbstractQuantity

Page average entropy of a subsystem for a Haar-random pure state in
`H_A ⊗ H_B`. For `dim(H_A) = m`, `dim(H_B) = n` with `m ≤ n` (else
swap by purity symmetry), Page 1993 found

    <S_A> = sum_{k=n+1}^{m·n} 1/k - (m-1)/(2n).

For `m = n` this gives `<S_A> ≈ log m - 1/2` (close to maximal but
reduced by 1/2); for `m << n` it gives `<S_A> ≈ log m - m/(2n)`.
This is the famous Page curve in dimension space underlying e.g. the
information-paradox / Page-time analysis of evaporating black holes.

Reference: D. N. Page, *Phys. Rev. Lett.* **71**, 1291 (1993),
DOI 10.1103/PhysRevLett.71.1291. Tracking: #580.
"""
struct PageEntropy <: AbstractEntanglementMeasure end

"""
    E8Spectrum() <: AbstractQuantity

Zamolodchikov E8 mass spectrum (8 stable particles).  Concrete
implementation lives in `src/universalities/E8.jl`; the type is defined
here so `src/core/alias.jl` can reference it without circular loads.
"""
struct E8Spectrum <: AbstractQuantity end

"""
    EdgeModeEnergy() <: AbstractQuantity

Energy of the lowest-lying boundary mode on an open chain.  In a
topological 1D superconductor (Kitaev 2001) the OBC chain hosts two
Majorana zero modes at the chain ends; their hybridization energy
decays exponentially with chain length,

```math
E_\\text{edge}(N) \\sim e^{-N/\\xi},
```

where `ξ` is the bulk correlation length.  In the trivial phase the
OBC lowest single-particle excitation is set by the bulk gap.

`EdgeModeEnergy` is the smallest positive BdG eigenvalue of the OBC
chain — the same quantity as [`MassGap`](@ref) at `OBC`, exposed under
a name that makes the boundary-mode interpretation explicit at the
call site.

Currently used by [`Kitaev1D`](@ref).
"""
struct EdgeModeEnergy <: AbstractQuantity end
# ─── Quench dynamics: Loschmidt echo / DQPT rate function ──────────────

"""
    LoschmidtEcho{M}() <: AbstractQuantity
    LoschmidtEcho(:amplitude)
    LoschmidtEcho(:rate)
    LoschmidtRateFunction()        # alias for LoschmidtEcho{:rate}

Loschmidt-echo family for sudden-quench dynamics.  After preparing
`|ψ_0⟩` as the ground state of an "initial" model `H_0` and quenching to
the "final" model `H_f` (passed as the first positional argument to
`fetch`), the Loschmidt amplitude is

    G(t) = ⟨ψ_0 | e^{-i H_f t} | ψ_0⟩,

with the Loschmidt echo `L(t) = |G(t)|² ∈ [0, 1]` and the rate function

    λ(t) = -log L(t) / N         (finite N)
    λ(t) = -lim_{N→∞} log L(t)/N (thermodynamic limit / Infinite)

Non-analytic cusps in `λ(t)` are dynamical quantum phase transitions
(DQPT).  See Heyl, Polkovnikov, Kehrein, [Heyl2013](@cite) and the
review Heyl, [Heyl2018](@cite).

The mode `M::Symbol ∈ (:amplitude, :rate)` is a phantom type parameter
so that `:amplitude` (returns `L(t)`) and `:rate` (returns `λ(t)`)
dispatch separately.  The convenience alias
[`LoschmidtRateFunction`](@ref) is the only flavour defined for
`Infinite`, since the `:amplitude` itself is identically zero in the
thermodynamic limit (extensive cumulants).

The pre-quench Hamiltonian is passed via the `initial` keyword on
`fetch`, e.g.

    fetch(TFIM(J=1.0, h=0.5), LoschmidtRateFunction(), Infinite();
          initial=TFIM(J=1.0, h=2.0), t=1.0)
"""
struct LoschmidtEcho{M} <: AbstractQuantity
    function LoschmidtEcho{M}() where {M}
        M isa Symbol || error("LoschmidtEcho mode must be a Symbol, got $(typeof(M))")
        M in (:amplitude, :rate) ||
            error("unknown LoschmidtEcho mode :$M; expected :amplitude or :rate")
        return new{M}()
    end
end
LoschmidtEcho(m::Symbol) = LoschmidtEcho{m}()
LoschmidtEcho(; mode::Symbol=:rate) = LoschmidtEcho{mode}()

"""
    const LoschmidtRateFunction = LoschmidtEcho{:rate}

Convenience alias for the rate-function flavour
`λ(t) = -log L(t)/N`.  See [`LoschmidtEcho`](@ref).
"""
const LoschmidtRateFunction = LoschmidtEcho{:rate}

# Other spectrum / universality tag types (`TightBindingSpectrum`,
# `ExactSpectrum`, `GroundStateEnergyDensity`, `CriticalExponents`,
# `GrowthExponents`) are currently defined in their respective model /
# universality source files as bare `struct X end`.  Later commits
# (M1.6-M1.8) subtype them to `AbstractQuantity` in place.

# ─── Energy granularity routing (depends on BoundaryCondition / _bc_size) ───
#
# `Energy()` is dispatch-routed to the model's native granularity through the
# `native_energy_granularity` trait.  The non-native granularity is provided
# by a generic conversion fallback that uses `_bc_size`.  Models on lattices
# whose system size is not encoded in `bc.N` (Kitaev's `Lx, Ly` kwargs) can
# define `fetch(::Model, ::Energy{:total}, bc; ...)` directly to bypass the
# fallback.

function fetch(
    model::AbstractQAtlasModel, ::Energy{:natural}, bc::BoundaryCondition; kwargs...
)
    g = native_energy_granularity(model, bc)
    return fetch(model, Energy{g}(), bc; kwargs...)
end

function fetch(
    model::AbstractQAtlasModel, ::Energy{:per_site}, bc::Union{OBC,PBC}; kwargs...
)
    g = native_energy_granularity(model, bc)
    g === :per_site && error(
        "QAtlas Energy(:per_site): $(typeof(model)) declares native :per_site at " *
        "$(typeof(bc)) but no direct method is registered.  Implement " *
        "`fetch(::$(typeof(model)), ::Energy{:per_site}, ::$(typeof(bc)); ...)` " *
        "to register it.",
    )
    return fetch(model, Energy{:total}(), bc; kwargs...) / _bc_size(bc, kwargs)
end

function fetch(model::AbstractQAtlasModel, ::Energy{:total}, bc::Union{OBC,PBC}; kwargs...)
    g = native_energy_granularity(model, bc)
    g === :total && error(
        "QAtlas Energy(:total): $(typeof(model)) declares native :total at " *
        "$(typeof(bc)) but no direct method is registered.  Implement " *
        "`fetch(::$(typeof(model)), ::Energy{:total}, ::$(typeof(bc)); ...)` " *
        "to register it.",
    )
    return fetch(model, Energy{:per_site}(), bc; kwargs...) * _bc_size(bc, kwargs)
end

# ─── Charge / spin gaps (correlated electron systems) ──────────────────

"""
    ChargeGap() <: AbstractQuantity

Charge (Mott) gap of an electron system,

    Δ_c = E₀(N+1) + E₀(N-1) - 2 E₀(N),

i.e. the energy cost of adding a particle plus the cost of removing
one, equivalent to the gap between the half-filled ground state and
the lowest charged excitation.  Strictly positive in a Mott insulator
and exactly zero in a metal / superconductor.

Implemented analytically for [`Hubbard1D`](@ref) at half filling via
the Lieb–Wu (1968) closed-form integral.
"""
struct ChargeGap <: AbstractGap end

"""
    SpinGap() <: AbstractQuantity

Spin gap of an electron system,

    Δ_s = E₀(S^z = 1) - E₀(S^z = 0),

i.e. the lowest excitation energy at fixed total particle number that
flips one spin.  Zero whenever the spinon branch is gapless (e.g. the
half-filled 1D Hubbard chain — rigorous Lieb–Wu result), positive in a
spin-gapped phase (Haldane chain, BCS superconductor, …).
"""
struct SpinGap <: AbstractGap end
# ─── Quench / nonequilibrium long-time ensembles ────────────────────────

"""
    GGEValue{Q<:AbstractQuantity}(inner) <: AbstractQuantity

Wrapper quantity carrying an underlying observable `inner::Q` whose
*generalised Gibbs ensemble* (GGE) stationary value is to be computed —
i.e. the `t → ∞` long-time average that an integrable (free-fermion)
quench reaches.

For an integrable system the ordinary (canonical) Gibbs ensemble does
not describe the long-time relaxed state: every mode-occupation
`n_k = ⟨c_k† c_k⟩` is a separate conserved quantity, so the diagonal
ensemble is a *generalised* Gibbs ensemble fixed by the full
distribution `{n_k}`.  See Rigol et al. [Rigol2007](@cite) for the
foundational argument and Calabrese, Essler, Fagotti J. Stat. Mech.
(2012) P07016 / P07022 for the TFIM-specific closed-form expressions.

`fetch(model_f, ::GGEValue{Q}, bc; initial::ModelType, kwargs...)`
returns the GGE expectation of the `Q` observable in the post-quench
Hamiltonian `model_f`, with the conserved mode occupations frozen by
the initial-state (`initial`) Bogoliubov rotation.

# Construction

```julia
GGEValue(Energy())                  # ⟨H_f⟩ stationary value
GGEValue(MagnetizationX())          # ⟨σˣ⟩ stationary value
```

# Fetch signature (TFIM)

```julia
fetch(TFIM(h = h_f), GGEValue(Energy()), Infinite();
      initial = TFIM(h = h_0)) -> Float64
```

A no-quench limit `h_0 = h_f` reduces to the static ground-state value
of the inner observable.
"""
struct GGEValue{Q<:AbstractQuantity} <: AbstractQuantity
    inner::Q
end
# core/quantities.jl — concrete quantity struct library.

# ─── Topological order (introduced for ToricCode, Kitaev 2003) ─────────

"""
    GroundStateDegeneracy() <: AbstractQuantity

Dimension of the ground-state subspace as an `Int`.  In topologically
ordered phases this is a robust, lattice-independent invariant
determined by the ambient surface (e.g. `4^g` on a closed orientable
genus-`g` surface for the toric code) and is set by the kwarg `genus`
on the `fetch` call.  Trivially `1` for any gapped, symmetry-unbroken
phase.
"""
struct GroundStateDegeneracy <: AbstractQuantity end

"""
    AnyonStatistics() <: AbstractQuantity

Topological data of the anyon content of a 2D topologically ordered
phase.  The dispatched `fetch` method takes a `type::Symbol` kwarg
selecting an anyon (or a mutual-braiding row) and returns a
`NamedTuple` whose shape depends on the requested row — typically
`(label, statistics, self_phase, quantum_dim, fusion)` for individual
anyons, and `(label, mutual_phase, anyons)` for two-anyon braids.

For Abelian theories like the toric code's Z₂ order all quantum
dimensions are 1; for non-Abelian theories the schema is unchanged but
`quantum_dim` becomes irrational and additional fields (e.g. `F`, `R`
matrices) may be added by the implementing method.

Has no boundary-condition argument: anyon statistics are purely
topological invariants.
"""
struct AnyonStatistics <: AbstractQuantity end

# ─── Random Matrix Theory level statistics (introduced for RMT/Poisson universality classes) ───

"""
    WignerSurmise() <: AbstractQuantity

Wigner surmise nearest-neighbour level-spacing distribution `P_β(s)`
for the three Wigner-Dyson ensembles (β ∈ {1, 2, 4}: GOE, GUE, GSE).
The surmise is the exact `N = 2` Gaussian-ensemble spacing
distribution; it is also a celebrated, accurate approximation to the
bulk `N → ∞` spacing distribution. Returns the value `P_β(s)` at
the requested `s` (the universality fetch carries `β`).
"""
struct WignerSurmise <: AbstractQuantity end

"""
    TracyWidom() <: AbstractQuantity

Tracy-Widom largest-eigenvalue cumulative distribution `F_β(x)` for
the three Wigner-Dyson ensembles (β ∈ {1, 2, 4}). Returns the value
`F_β(x) = P[ξ_β ≤ x]` at the requested `x`.

QAtlas Phase 1 evaluates `F_β` from a precomputed table compiled
from Bornemann, *On the numerical evaluation of Fredholm
determinants*, Math. Comp. **79**, 871 (2010), Table 1, with
monotone linear interpolation on the table support and
Tracy-Widom 1994/1996 tail asymptotics outside it. A direct
Painlevé-II integrator is deferred to Phase 2.
"""
struct TracyWidom <: AbstractQuantity end

"""
    MeanRatio() <: AbstractQuantity

Mean of the consecutive level-spacing ratio
`r_n = min(s_n, s_{n+1}) / max(s_n, s_{n+1})`,
introduced by Oganesyan-Huse (2007) and tabulated for the
Wigner-Dyson and Poisson ensembles by Atas-Bogomolny-Giraud-Roux,
Phys. Rev. Lett. **110**, 084101 (2013):

| Ensemble  | ⟨r⟩                  |
|-----------|----------------------|
| Poisson   | 2 log 2 − 1 ≈ 0.3863 |
| GOE (β=1) | 0.5307               |
| GUE (β=2) | 0.5996               |
| GSE (β=4) | 0.6744               |
"""
struct MeanRatio <: AbstractQuantity end

# ─── CFT finite-size scaling (introduced for Universality{:Ising,...} classes) ───

"""
    CasimirEnergyCorrection() <: AbstractQuantity

Universal `1/L` finite-size correction to the ground-state energy of a
1+1D conformal field theory.

For a critical 1+1D system with central charge `c` and CFT velocity
`v` on a system of size `L`:

- Periodic boundary (PBC):
  ``E_0(L) = L\\,\\varepsilon_\\infty - \\dfrac{\\pi c v}{6 L} + O(L^{-2})``
- Open boundary (OBC):
  ``E_0(L) = L\\,\\varepsilon_\\infty + \\varepsilon_{\\mathrm{surf}} - \\dfrac{\\pi c v}{24 L} + O(L^{-2})``

This quantity returns *only* the universal ``1/L`` correction term
(``-\\pi c v/(6 L)`` at PBC, ``-\\pi c v/(24 L)`` at OBC), not the
extensive ``L \\varepsilon_\\infty`` piece nor the OBC surface term
``\\varepsilon_{\\mathrm{surf}}``.  The PBC-to-OBC ratio is exactly 4,
independent of the universality class.

The CFT velocity `v` is model-dependent (e.g. ``v = 2J`` for the TFIM
at the critical point, ``v = (\\pi/2) J`` for the AFM Heisenberg chain,
``v = v_F`` for the XXZ Luttinger liquid) and is supplied by the caller
as a kwarg.  The central charge `c` is read from the universality
class via the same data the `Universality{C}` entry exposes for
[`CriticalExponents`](@ref).

# References
- J. Cardy, *Nucl. Phys. B* **270**, 186 (1986).
- H. W. J. Blöte, J. L. Cardy, M. P. Nightingale, *Phys. Rev. Lett.*
  **56**, 742 (1986).
- I. Affleck, *Phys. Rev. Lett.* **56**, 746 (1986).

!!! note "Phase 2 (TODO)"
    The conformal *tower of states* --- primary scaling dimensions
    ``(h, \\bar h)`` and the
    ``E_n - E_0 = (2\\pi v/L)(h_n + \\bar h_n)`` excitation pattern ---
    is tracked separately as future work (Phase 2 of issue #150) and
    will be exposed via a `ConformalTower` quantity once implemented.
"""
struct CasimirEnergyCorrection <: AbstractQuantity end

@doc raw"""
    ConformalTower() <: AbstractQuantity

The conformal tower of states excitation spectrum in 1+1D conformal field theories.
At boundary condition `bc::Union{PBC, OBC}`, returns a sorted vector of NamedTuples
representing the lowest-lying excitation energies `E_n - E_0` relative to the
ground state, their scaling dimensions `Δ_n` (or `h_n`), and their degeneracies:

    (energy = E_n - E_0, dimension = Δ_n, degeneracy = g_n)

For periodic boundary conditions (PBC), the excitation energies scale as:
    E_n - E_0 = (2π v / L) Δ_n
where `Δ_n = h_n + \bar{h}_n` is the scaling dimension.

For open boundary conditions (OBC), the excitation energies scale as:
    E_n - E_0 = (π v / L) h_n
where `h_n` is the boundary scaling dimension.

Keyword arguments:
- `L::Real`: system size.
- `v::Real`: CFT sound velocity.
"""
struct ConformalTower <: AbstractQuantity end

"""
    SteadyStateCurrent() <: AbstractQuantity

Steady-state mass / particle current `j` in a 1D non-equilibrium lattice
gas (e.g. ASEP / TASEP, Derrida-Lebowitz 1998).  For TASEP at hopping
rate `p` and density `ρ`,

    j(ρ) = p ρ (1 − ρ)              (TASEP mean-field steady state)

— the canonical KPZ-class non-equilibrium observable.
"""
struct SteadyStateCurrent <: AbstractQuantity end

"""
    ChiralCondensate() <: AbstractQuantity

Vacuum expectation value `⟨ψ̄ψ⟩` of a fermion bilinear, signalling
spontaneous (anomalous) chiral-symmetry breaking.  The massless
Schwinger model is the canonical 1+1-D example: even though the
classical Lagrangian is chirally symmetric, the anomaly forces a
non-zero condensate

    ⟨ψ̄ψ⟩ = − exp(γ_E) · e / (2π^{3/2}),    m_γ = e/√π.

(Schwinger 1962; Coleman-Jackiw-Susskind 1975.)
"""
struct ChiralCondensate <: AbstractQuantity end
# ─── RMT spectral form factor (introduced for Universality{:RMT}, issue #243) ─

"""
    SpectralFormFactor() <: AbstractQuantity

Disorder-averaged spectral form factor
`K(t) = ⟨|Σ_n e^{−iE_n t}|²⟩ / Z²`
— the canonical late-time quantum-chaology diagnostic
(Mehta 2004 §16; Cotler et al. 2017).

For GUE random-matrix-theory eigenvalues in the large-`N`
thermodynamic limit, with rescaled time `τ = t / N`, the
disorder-averaged SFF has the universal closed form

* `K(τ) = (τ/(2π)) − (τ/(4π)) log|1 − τ/(2π)|`   for `τ ≤ 2π`
* `K(τ) = 1`                                       for `τ ≥ 2π`

so that `K` exhibits a linear ramp `K(τ) ≈ τ/π` for small `τ` and
saturates to the universal plateau `K(τ→∞) = 1` for `τ` beyond
the Heisenberg time `τ_H = 2π`.

QAtlas Phase 1 (issue #243) exposes only the late-time plateau
`τ ≥ τ_H` for the GUE ensemble; the ramp regime and the GOE/GSE
sigma-model closed forms (Mehta 2004 §16) are deferred to Phase 2.

# References
- M. L. Mehta, *Random Matrices*, 3rd ed., Elsevier (2004), §16.
- E. Brézin, S. Hikami, *Phys. Rev. E* **55**, 4067 (1997).
- J. S. Cotler, G. Gur-Ari, M. Hanada, J. Polchinski, P. Saad,
  S. H. Shenker, D. Stanford, A. Streicher, M. Tezuka,
  *JHEP* **05**, 118 (2017), arXiv:1611.04650 — ramp-plateau picture.
"""
struct SpectralFormFactor <: AbstractQuantity end

# ─── Registered-status worked example (status axis, v0.24) ──────────────
# A one-sided BOUND quantity (Lieb-Robinson velocity cone) carried in with
# the new `status` registry axis — exercising the :bound status,
# verify_bound, and the atlas status rendering.  (The :approx status and
# verify_approx ship too; their worked example is deferred to the
# definition-list redesign, which expresses approximations as non-canonical
# definitions of an existing quantity rather than a new quantity.)

"""
    LiebRobinsonBound() <: AbstractQuantity

Lieb-Robinson velocity `v_LR` — the slope of the causal cone bounding the
spread of operator commutators,

    ‖[A_x(t), B_y(0)]‖ ≤ C · exp(−μ (|x − y| − v_LR · t)).

`fetch` returns `v_LR` itself (the bounding cone slope). This is a
one-sided `:bound`: any genuinely measured information velocity stays
`≤ v_LR`, and for free-fermion models the bound is *saturated* by the
maximum group velocity `max_k |dΛ/dk|`. Registered with `status=:bound`.

(Lieb & Robinson 1972; Hastings & Koma 2006.)
"""
struct LiebRobinsonBound <: AbstractQuantity end

"""
    CHSHBound() <: AbstractQuantity

The CHSH (Bell-inequality) correlator bound — the maximum of
`S = E(a,b) + E(a,b′) + E(a′,b) − E(a′,b′)` admissible in a given physical
theory.  A `status=:bound` quantity with the historical name (like
[`LiebRobinsonBound`](@ref)); fetched against a [`Bound`](@ref) domain
(not a model), with a `scheme=` selector choosing the theory regime
(`:classical` → 2, `:quantum` → 2√2, `:no_signalling` → 4).
"""
struct CHSHBound <: AbstractQuantity end

"""
    MerminGHZBound() <: AbstractQuantity

The Mermin 3-party Bell-type bound — the maximum of the Mermin operator
`|⟨M₃⟩|` admissible in a given theory.  A `status=:bound` quantity
(Mermin 1990); fetched against a [`Bound`](@ref) domain with `scheme=`
choosing the theory regime (`:classical` → 2 local-realistic, `:quantum`
→ 4 quantum, saturated by the GHZ state).
"""
struct MerminGHZBound <: AbstractQuantity end

"""
    ChaosBound() <: AbstractQuantity

The Maldacena–Shenker–Stanford bound on quantum chaos — an upper bound on the
Lyapunov exponent `λ_L` of out-of-time-order correlators (`λ_L ≤ 2π/β`).  A
`status=:bound` quantity; fetched against a [`Bound`](@ref) domain
(`Bound(:Dynamics)`).
"""
struct ChaosBound <: AbstractQuantity end

"""
    ScramblingTime() <: AbstractQuantity

The fast-scrambling time `t_* = (β/2π) log N` (Sekino–Susskind 2008) — the
conjectured *lower* bound on the time for a thermal system of `N` degrees of
freedom to scramble local information into global entanglement; saturated by
black holes (the fastest scramblers).  A `status=:bound`, `direction=:lower`
quantity; fetched against a [`Bound`](@ref) domain (`Bound(:Dynamics)`).
"""
struct ScramblingTime <: AbstractQuantity end

"""
    BekensteinBound() <: AbstractQuantity

The Bekenstein universal entropy bound — an upper bound on the entropy of a
bounded system (`S ≤ 2π R E`).  A `status=:bound` quantity; fetched against a
[`Bound`](@ref) domain (`Bound(:Holographic)`).
"""
struct BekensteinBound <: AbstractQuantity end

"""
    QuantumSpeedLimit() <: AbstractQuantity

The quantum speed limit — a *lower* bound on the time to evolve a state to an
orthogonal one (Margolus–Levitin `τ ≥ π/(2E)`).  A `status=:bound`,
`direction=:lower` quantity; fetched against a [`Bound`](@ref) domain
(`Bound(:Dynamics)`).
"""
struct QuantumSpeedLimit <: AbstractQuantity end

"""
    OptimalCloningFidelity() <: AbstractQuantity

The optimal universal quantum cloning fidelity — an upper bound on the
single-copy fidelity of a `1 → 2` qubit cloner (Bužek–Hillery `F ≤ 5/6`).  A
`status=:bound`, `direction=:upper` quantity; fetched against a [`Bound`](@ref)
domain (`Bound(:QuantumInformation)`).
"""
struct OptimalCloningFidelity <: AbstractQuantity end

"""
    BB84KeyRate() <: AbstractQuantity

The BB84 asymptotic secret-key rate `R(e) = 1 − 2 H₂(e)` (Shor–Preskill 2000),
with `H₂` the binary entropy and `e` the qubit error rate (QBER).  A provably
achievable rate — a *lower* bound on the extractable secret-key fraction;
positive for `e < 11%`.  A `status=:bound`, `direction=:lower` quantity; fetched
against a [`Bound`](@ref) domain (`Bound(:QuantumInformation)`).
"""
struct BB84KeyRate <: AbstractQuantity end

"""
    Polarization() <: AbstractQuantity

The bulk polarization density (or order parameter) per site. For the
classical 2D six-vertex model, it corresponds to the spontaneous polarization
(in the ferroelectric phase Δ > 1) or the spontaneous staggered polarization
(in the antiferroelectric phase Δ < -1).
"""
struct Polarization <: AbstractQuantity end

@doc raw"""
    SphereFreeEnergy() <: AbstractQuantity

Universal sphere free energy $F = -\ln |Z(S^3)|$ of a 2+1D conformal field theory.
Acts as a measure of the degrees of freedom in 2+1D (the $F$-theorem).
"""
struct SphereFreeEnergy <: AbstractQuantity end

@doc raw"""
    CornerEntanglementCoefficient() <: AbstractQuantity

Universal corner coefficient in the bipartite entanglement entropy of a 2+1D CFT
with a boundary corner of angle $\theta$:

    S(ρ_A) = a |∂A| - c(\theta) \ln(L/\epsilon) + o(\ln(L/\epsilon)).

For a nearly smooth boundary $\theta \to \pi$, $c(\theta) \approx \sigma (\pi - \theta)^2$
where $\sigma = \frac{\pi^2}{24} C_T$. If no angle `theta` is provided, the fetch
method returns the smooth-limit prefactor $\sigma$.
"""
struct CornerEntanglementCoefficient <: AbstractQuantity end

"""
    DynamicLocalization() <: AbstractQuantity

Cycle-averaged effective-hopping renormalization of a tight-binding band driven
by a spatially-uniform monochromatic ac electric field (Peierls coupling).  In
units `e = ℏ = a = 1`, a field `E(τ) = E₀ cos(ωτ)` gives the dimensionless drive
`K = E₀/ω`, and the hopping is renormalized by the exact, nonperturbative Bessel
factor

    t_eff / t = J₀(K).

The band collapses — "dynamic localization" — at the zeros of `J₀` (first at
`K = 2.404826…`), where a static tilt drives no current despite `E₀ ≠ 0`.  This
is the hallmark exact nonlinear (all-orders-in-field) response of the ac-driven
free-fermion chain (Dunlap–Kenkre 1986; Holthaus–Hone 1996); the full harmonic
content of the current is the Bessel spectrum `Jₙ(K)`
(see `driven_band_harmonic_weights`).
"""
struct DynamicLocalization <: AbstractQuantity end

"""
    HighHarmonicAmplitude() <: AbstractQuantity

Peak amplitude of the `harmonic`-th harmonic (frequency `n ω`) of the intraband
current of a tight-binding band driven by a monochromatic ac field — the
exact, all-orders-in-field higher-order response (high-harmonic generation).

For drive `K = E₀/ω` the n-th harmonic amplitude of the current, maximized over
crystal momentum, is the Bessel envelope

    A₀(K) = 2t |J₀(K)|,      Aₙ(K) = 4t |Jₙ(K)|   (n ≥ 1),

so `harmonic = 1` is the linear response, `harmonic ≥ 2` the genuinely nonlinear
higher harmonics.  For small `K`, `Aₙ ∝ Kⁿ` — the n-th harmonic is the order-n
(χ⁽ⁿ⁾) response, whose leading coefficient is `nonlinear_susceptibility`.  The
`n = 0` value is the dynamic-localization envelope (`2·|`[`DynamicLocalization`](@ref)`|`).
(Dunlap–Kenkre 1986; Holthaus–Hone 1996.)
"""
struct HighHarmonicAmplitude <: AbstractQuantity end

# ─── component trait: concrete methods (#690) ────────────────────────────
# The component a leaf's type name encodes; `nothing` (the AbstractQuantity
# default above) everywhere else.  The `…Local` magnetizations deliberately
# define NO component: their extra site argument makes them a different
# fetch shape, so component-paired identity generation must not pick them up.

component(::Type{MagnetizationX}) = :x
component(::Type{MagnetizationY}) = :y
component(::Type{MagnetizationZ}) = :z
component(::Type{SusceptibilityXX}) = :xx
component(::Type{SusceptibilityYY}) = :yy
component(::Type{SusceptibilityZZ}) = :zz
component(::Type{<:XXCorrelation}) = :xx
component(::Type{<:YYCorrelation}) = :yy
component(::Type{<:ZZCorrelation}) = :zz
component(::Type{XXStructureFactor}) = :xx
component(::Type{YYStructureFactor}) = :yy
component(::Type{ZZStructureFactor}) = :zz
component(::Type{MassGap}) = :mass
component(::Type{ChargeGap}) = :charge
component(::Type{SpinGap}) = :spin
