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

# ─── Scalar thermodynamics ──────────────────────────────────────────────

"""
    Energy{G}() <: AbstractQuantity
    Energy()                 # G = :natural — model-and-BC-natural granularity
    Energy(:total)           # explicit ⟨H⟩
    Energy(:per_site)        # explicit ⟨H⟩ / N

Ground-state / thermal energy expectation.  The type parameter `G` makes
the granularity (total vs per-site) a dispatch axis instead of a hidden
docstring contract.

`Energy()` resolves to the model's native granularity via the
[`native_energy_granularity`](@ref) trait — keeping every existing
`fetch(model, Energy(), bc; ...)` call site working unchanged.  Use the
explicit constructors when the caller needs a specific granularity (e.g.
the thermodynamic-identity harness comparing `f + T·s` against per-site
`ε`).

The non-native granularity is provided automatically by a generic
conversion fallback for 1D BCs (`OBC` / `PBC`) that uses
[`_bc_size`](@ref).  Models on lattices whose size is not captured by
`bc.N` (e.g. 2D Kitaev with `Lx, Ly` kwargs) currently support only
their declared native granularity.
"""
struct Energy{G} <: AbstractQuantity
    function Energy{G}() where {G}
        G isa Symbol || error("Energy granularity must be a Symbol, got $(typeof(G))")
        G in (:natural, :total, :per_site) ||
            error("unknown Energy granularity :$G; expected :natural, :total, or :per_site")
        return new{G}()
    end
end
Energy() = Energy{:natural}()
Energy(g::Symbol) = Energy{g}()

"""
    native_energy_granularity(model, bc) -> :total | :per_site

Trait declaring which granularity the given `model` returns natively for
[`Energy`](@ref) at boundary condition `bc`.  Used by the `Energy()`
(`:natural`) router and by the generic conversion fallbacks.

Every model that supports `Energy` must add a method per supported BC,
e.g.

```julia
QAtlas.native_energy_granularity(::TFIM, ::OBC) = :total
QAtlas.native_energy_granularity(::TFIM, ::Infinite) = :per_site
```

A missing method is caught at the call site as a `MethodError`, which
is intentional: it forces new models to declare the convention rather
than silently inheriting an unrelated default.
"""
function native_energy_granularity end

"""
    FreeEnergy() <: AbstractQuantity

Helmholtz free energy per site, `f = -β⁻¹ log Z / N`.
"""
struct FreeEnergy <: AbstractQuantity end

"""
    SpecificHeat() <: AbstractQuantity

Specific heat per site, `c_v(β) = β² (⟨H²⟩ − ⟨H⟩²) / N`.
"""
struct SpecificHeat <: AbstractQuantity end

"""
    MassGap() <: AbstractQuantity

Energy gap between the ground state and the first excited state.
"""
struct MassGap <: AbstractQuantity end

"""
    FidelitySusceptibility() <: AbstractQuantity

Fidelity susceptibility `χ_F(λ) = −∂²⟨ψ(λ)|ψ(λ + δλ)⟩/∂δλ²`.
"""
struct FidelitySusceptibility <: AbstractQuantity end

"""
    LoschmidtEcho{M}() <: AbstractQuantity
    LoschmidtEcho(; mode::Symbol = :rate)

Quench observable derived from the Loschmidt amplitude
`L(t) = ⟨ψ₀ | e^{-i H_f t} | ψ₀⟩` of an initial state `|ψ₀⟩` evolving
under a final Hamiltonian `H_f`.  The phantom-type parameter `M`
selects which derived quantity is returned:

- `:probability`  →  `|L(t)|²` (the Loschmidt echo proper, in `[0, 1]`)
- `:rate`          →  Loschmidt rate function
                      `λ(t) = -lim_{N→∞} (1/N) log |L(t)|²`,
                      the standard order parameter for dynamical
                      quantum phase transitions in the thermodynamic
                      limit (Heyl, Polkovnikov & Kehrein 2013)

A model-specific `fetch` method is required to provide the initial
Hamiltonian (typically as the kwarg `initial::Model`) and the time
`t::Real`; a Δ-mismatch / DomainError convention is left to the model
author.

This type is shared infrastructure for issues #143 (TFIM Loschmidt),
#148 (XX free-fermion quench), and any future quench observable that
factorises through the Loschmidt amplitude.

# References
- M. Heyl, A. Polkovnikov, S. Kehrein, *Phys. Rev. Lett.* 110, 135704
  (2013).
- F.H.L. Essler, M. Fagotti, *J. Stat. Mech.* (2016) 064002.
"""
struct LoschmidtEcho{M} <: AbstractQuantity
    function LoschmidtEcho{M}() where {M}
        M isa Symbol ||
            throw(ArgumentError("LoschmidtEcho mode must be a Symbol, got $(typeof(M))"))
        M in (:probability, :rate) || throw(
            ArgumentError("unknown LoschmidtEcho mode :$M; expected :probability or :rate"),
        )
        return new{M}()
    end
end
LoschmidtEcho(; mode::Symbol=:rate) = LoschmidtEcho{mode}()

# `PartitionFunction`, `CriticalTemperature`, `SpontaneousMagnetization`
# are currently defined in src/models/classical/IsingSquare/IsingSquare.jl
# as bare `struct X end` tags.  They will be migrated to subtype
# `AbstractQuantity` in the IsingSquare refactor commit (M1.7).

# ─── Entropies (explicit variants; see user-requested naming) ──────────

"""
    ThermalEntropy() <: AbstractQuantity

Thermal / thermodynamic entropy per site, `s(β) = −∂f/∂T` where `f` is the
free energy per site.  Real-valued, non-negative, monotone in `T`.
"""
struct ThermalEntropy <: AbstractQuantity end

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
struct VonNeumannEntropy{M} <: AbstractQuantity
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
struct RenyiEntropy <: AbstractQuantity
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
struct ResidualEntropy <: AbstractQuantity end

# ─── Magnetizations (axis explicit) ─────────────────────────────────────

"""
    MagnetizationX() <: AbstractQuantity

Bulk-averaged `⟨σˣ⟩` in Pauli convention (= 2 ⟨Sˣ⟩ in spin-1/2 units).
For a spin-1/2 chain `H = -J ΣSᶻSᶻ - h ΣSˣ` this is the transverse
magnetization; the axis-explicit name avoids the "transverse" /
"longitudinal" ambiguity that depends on the model's Hamiltonian
choice.
"""
struct MagnetizationX <: AbstractQuantity end

"""
    MagnetizationY() <: AbstractQuantity

Bulk-averaged `⟨σʸ⟩`.
"""
struct MagnetizationY <: AbstractQuantity end

"""
    MagnetizationZ() <: AbstractQuantity

Bulk-averaged `⟨σᶻ⟩`.  For Z₂-symmetric phases on an infinite system
this is the order parameter at low temperature; finite-system fetch
methods may return the absolute value / the ordered-phase limit as
documented.
"""
struct MagnetizationZ <: AbstractQuantity end

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
struct MagnetizationXLocal{M} <: AbstractQuantity
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
struct MagnetizationYLocal <: AbstractQuantity end

"""
    MagnetizationZLocal() <: AbstractQuantity

Site-resolved `⟨σᶻ_i⟩` vector of length `N_bulk`.
"""
struct MagnetizationZLocal <: AbstractQuantity end

"""
    EnergyLocal() <: AbstractQuantity

Bond-resolved energy density vector, length `N_bulk − 1` for a bond
Hamiltonian `Σ_b h_b`.
"""
struct EnergyLocal <: AbstractQuantity end

# ─── Susceptibilities (axis pair) ────────────────────────────────────────

"""
    SusceptibilityXX() <: AbstractQuantity

Static transverse susceptibility,
`χ_xx(β) = β · (⟨M_x²⟩ − ⟨M_x⟩²) / N`.
"""
struct SusceptibilityXX <: AbstractQuantity end

"""
    SusceptibilityYY() <: AbstractQuantity

Analogue for the y-axis.
"""
struct SusceptibilityYY <: AbstractQuantity end

"""
    SusceptibilityZZ() <: AbstractQuantity

Uniform longitudinal susceptibility,
`χ_zz(β) = β · (⟨M_z²⟩ − ⟨M_z⟩²) / N`.
"""
struct SusceptibilityZZ <: AbstractQuantity end

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
struct ZZCorrelation{M} <: AbstractQuantity end
ZZCorrelation(; mode::Symbol=:static) = ZZCorrelation{mode}()

"""
    XXCorrelation{M}() <: AbstractQuantity
    XXCorrelation(; mode::Symbol = :static)

Real-space 2-point `⟨σˣ_i σˣ_j⟩` correlator.  See
[`ZZCorrelation`](@ref) for the `mode` semantics.
"""
struct XXCorrelation{M} <: AbstractQuantity end
XXCorrelation(; mode::Symbol=:static) = XXCorrelation{mode}()

"""
    YYCorrelation{M}() <: AbstractQuantity
    YYCorrelation(; mode::Symbol = :static)

Real-space 2-point `⟨σʸ_i σʸ_j⟩` correlator.
"""
struct YYCorrelation{M} <: AbstractQuantity end
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
struct ZZStructureFactor <: AbstractQuantity end

"""
    XXStructureFactor() <: AbstractQuantity

Fourier-space equivalent of [`XXCorrelation`](@ref).
"""
struct XXStructureFactor <: AbstractQuantity end

"""
    YYStructureFactor() <: AbstractQuantity

Fourier-space equivalent of [`YYCorrelation`](@ref).
"""
struct YYStructureFactor <: AbstractQuantity end

# ─── Universality / lattice spectra / advanced ─────────────────────────

"""
    CentralCharge() <: AbstractQuantity

Central charge `c` of the emergent CFT.  For 1D critical systems
extracted from the Calabrese–Cardy entanglement formula; universality
pages return literature values.
"""
struct CentralCharge <: AbstractQuantity end

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
    CorrelationLength() <: AbstractQuantity

Two-point correlation length `ξ` controlling the exponential decay of
connected equal-time correlators in a gapped phase,

    ⟨σ_α(0) σ_α(r)⟩_c ~ e^{-r/ξ}    (r → ∞).

For a critical system `ξ = ∞`; implementations return `Inf` in that
case.  At `T = 0` and 1D free-fermion models like TFIM, `ξ` is set by
the inverse mass gap (`ξ = 1/(2|h - J|)`).
"""
struct CorrelationLength <: AbstractQuantity end

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
returned by models like [`Honeycomb`](@ref) (at the Dirac cones), the
other tight-binding lattices, and the TFIM Majorana mode at the
critical field.
"""
struct FermiVelocity <: AbstractQuantity end

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
struct LuttingerVelocity <: AbstractQuantity end

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
    E8Spectrum() <: AbstractQuantity

Zamolodchikov E8 mass spectrum (8 stable particles).  Concrete
implementation lives in `src/universalities/E8.jl`; the type is defined
here so `src/core/alias.jl` can reference it without circular loads.
"""
struct E8Spectrum <: AbstractQuantity end

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
(DQPT).  See Heyl, Polkovnikov, Kehrein, PRL 110, 135704 (2013) and the
review Heyl, Rep. Prog. Phys. 81, 054001 (2018).

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
LoschmidtEcho() = LoschmidtEcho{:rate}()
LoschmidtEcho(m::Symbol) = LoschmidtEcho{m}()

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

Implemented analytically for [](@ref) at half filling via
the Lieb–Wu (1968) closed-form integral.
"""
struct ChargeGap <: AbstractQuantity end

"""
    SpinGap() <: AbstractQuantity

Spin gap of an electron system,

    Δ_s = E₀(S^z = 1) - E₀(S^z = 0),

i.e. the lowest excitation energy at fixed total particle number that
flips one spin.  Zero whenever the spinon branch is gapless (e.g. the
half-filled 1D Hubbard chain — rigorous Lieb–Wu result), positive in a
spin-gapped phase (Haldane chain, BCS superconductor, …).
"""
struct SpinGap <: AbstractQuantity end
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
distribution `{n_k}`.  See Rigol et al. PRL 98, 050405 (2007) for the
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
