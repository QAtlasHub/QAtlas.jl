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
    VonNeumannEntropy() <: AbstractQuantity

Von Neumann entanglement entropy of a reduced density matrix:
`S_vN = −Tr ρ_A log ρ_A`.  Requires a subsystem specification through the
model's fetch kwargs (e.g. `ℓ`, the subsystem length).
"""
struct VonNeumannEntropy <: AbstractQuantity end

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
    MagnetizationXLocal() <: AbstractQuantity

Site-resolved `⟨σˣ_i⟩` vector of length `N_bulk`.
"""
struct MagnetizationXLocal <: AbstractQuantity end

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
