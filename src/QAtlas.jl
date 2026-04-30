module QAtlas

export AbstractModel, Model
export BoundaryCondition, Infinite, PBC, OBC
export AbstractQuantity, Quantity
export fetch

# --- Classical Models ---
export IsingSquare, PartitionFunction, CriticalTemperature, SpontaneousMagnetization

# --- Quantum Models ---
export TFIM                                             # v0.13 concrete struct
export E8                                               # v0.13 concrete struct
export XXZ1D                                            # v0.13 new model
export KitaevHoneycomb                                  # spin-½ Kitaev honeycomb
export TightBindingSpectrum
# NOTE: `Honeycomb`, `Kagome`, `Lieb`, `Triangular` are NOT exported —
# they all conflict with Lattice2D's topology types of the same name.
# Access them as `QAtlas.Honeycomb()` / `QAtlas.Kagome()` / etc. in code
# that also uses `Lattice2D`.  `Graphene` *is* exported as the
# backward-compat top-level alias for `Honeycomb` (see src/deprecate/)
# since the name does not collide with anything in Lattice2D.
export Heisenberg1D, ExactSpectrum, GroundStateEnergyDensity
export S1Heisenberg1D                                    # spin-1 (Haldane chain)

# --- Core Implementation ---
include("core/alias.jl")
include("core/type.jl")
include("core/quantities.jl")
include("core/registry.jl")
include("core/pfaffian.jl")
include("core/dense_ed.jl")

# --- Implementation registry public API ---
export Implementation, implementation_status, implementation_status_markdown

# --- Quantity struct exports (new, axis-explicit naming) ---
export Energy, FreeEnergy, SpecificHeat, MassGap, FidelitySusceptibility
export ThermalEntropy, VonNeumannEntropy, RenyiEntropy
export MagnetizationX, MagnetizationY, MagnetizationZ
export MagnetizationXLocal, MagnetizationYLocal, MagnetizationZLocal, EnergyLocal
export SusceptibilityXX, SusceptibilityYY, SusceptibilityZZ
export XXCorrelation, YYCorrelation, ZZCorrelation
export XXStructureFactor, YYStructureFactor, ZZStructureFactor
export CentralCharge, LuttingerParameter, CorrelationLength
export FermiVelocity, LuttingerVelocity, SpinWaveVelocity
export E8Spectrum

# --- TFIM Infinite dynamic helpers ---
export tfim_quasiparticle_dispersion, tfim_two_spinon_dos

# --- Universality Classes ---
export Universality, CriticalExponents, GrowthExponents
export Ising2D, KPZ1D, MeanField  # backward-compatible aliases
include("universalities/Universality.jl")
include("universalities/E8.jl")
include("universalities/MeanField.jl")
include("universalities/Ising2D.jl")
include("universalities/KPZ.jl")
include("universalities/Percolation.jl")
include("universalities/Potts.jl")
include("universalities/ONModel.jl")

# --- Models ---
# Layout: `<class>/<Model>/<Model>.jl` (with optional sibling axis files like
# `TFIM/TFIM_thermal.jl`).  `tightbinding/<lattice-class>/` groups multiple
# tight-binding Hamiltonians by lattice type (regular = Bloch-diagonalisable;
# future: quasicrystalline, fractal, disordered).
include("models/classical/IsingSquare/IsingSquare.jl")
include("models/quantum/tightbinding/regular/Honeycomb.jl")
include("models/quantum/tightbinding/regular/Kagome.jl")
include("models/quantum/tightbinding/regular/Lieb.jl")
include("models/quantum/tightbinding/regular/Triangular.jl")
include("models/quantum/TFIM/TFIM.jl")
include("models/quantum/TFIM/TFIM_dynamics.jl")
include("models/quantum/TFIM/TFIM_xx_static.jl")
include("models/quantum/TFIM/TFIM_thermal.jl")
include("models/quantum/TFIM/TFIM_pbc_thermal.jl")
include("models/quantum/TFIM/TFIM_zaxis.jl")
include("models/quantum/TFIM/TFIM_local.jl")
include("models/quantum/TFIM/TFIM_entanglement.jl")
include("models/quantum/TFIM/TFIM_cft_entanglement.jl")
include("models/quantum/TFIM/TFIM_infinite_dynamics.jl")
include("models/quantum/TFIM/TFIM_registry.jl")  # populates REGISTRY for TFIM
include("models/quantum/Heisenberg/Heisenberg.jl")
include("models/quantum/Heisenberg/HeisenbergS1.jl")
include("models/quantum/Heisenberg/HeisenbergS1_observables.jl")
include("models/quantum/Heisenberg/HeisenbergS1_registry.jl")
include("models/quantum/KitaevHoneycomb/KitaevHoneycomb.jl")
include("models/quantum/KitaevHoneycomb/KitaevHoneycomb_thermal.jl")
include("models/quantum/KitaevHoneycomb/KitaevHoneycomb_registry.jl")
include("models/quantum/XXZ/XXZ.jl")
include("models/quantum/XXZ/XXZ_thermal.jl")
include("models/quantum/XXZ/XXZ_registry.jl")  # populates REGISTRY for XXZ1D
include("models/quantum/Heisenberg/Heisenberg_registry.jl")  # populates REGISTRY for Heisenberg1D

# --- Deprecation shims (legacy API) ---
# Loaded last so they can route into any already-registered concrete
# `fetch` method.  See src/deprecate/README.md.
include("deprecate/legacy_fetch.jl")
include("deprecate/legacy_tfim.jl")
include("deprecate/legacy_e8.jl")
include("deprecate/legacy_honeycomb.jl")
export Graphene                                         # backward-compat alias
include("deprecate/legacy_xxz.jl")

end # module QAtlas
