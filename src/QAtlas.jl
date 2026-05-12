module QAtlas

export AbstractModel, Model
export BoundaryCondition, Infinite, PBC, OBC
export AbstractQuantity, Quantity
export fetch

# --- Classical Models ---
export IsingSquare, PartitionFunction, CriticalTemperature, SpontaneousMagnetization
export SixVertex
export IsingTriangular

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
export Hubbard1D                                         # Lieb-Wu Bethe ansatz half-filling
export MajumdarGhosh                                     # spin-1/2 J1-J2 chain at MG point
export S1Heisenberg1D                                    # spin-1 (Haldane chain)
export AKLT1D                                            # spin-1 BLBQ at AKLT point

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
export Energy, FreeEnergy, SpecificHeat, MassGap, FidelitySusceptibility, LoschmidtEcho
export ThermalEntropy, VonNeumannEntropy, RenyiEntropy
export Energy, FreeEnergy, SpecificHeat, MassGap, FidelitySusceptibility
export ChargeGap, SpinGap                                # Hubbard / correlated-electron gaps
export ThermalEntropy, VonNeumannEntropy, RenyiEntropy
export ThermalEntropy, VonNeumannEntropy, RenyiEntropy, ResidualEntropy
export MagnetizationX, MagnetizationY, MagnetizationZ
export MagnetizationXLocal, MagnetizationYLocal, MagnetizationZLocal, EnergyLocal
export SusceptibilityXX, SusceptibilityYY, SusceptibilityZZ
export XXCorrelation, YYCorrelation, ZZCorrelation
export XXStructureFactor, YYStructureFactor, ZZStructureFactor
export CentralCharge, LuttingerParameter, CorrelationLength
export ConformalWeights, PrimaryFields
export StringOrderParameter
export FermiVelocity, LuttingerVelocity, SpinWaveVelocity
export E8Spectrum
export LoschmidtEcho, LoschmidtRateFunction
export GGEValue                                          # quench long-time wrapper

# --- TFIM Infinite dynamic helpers ---
export tfim_quasiparticle_dispersion, tfim_two_spinon_dos

# --- Heisenberg1D Infinite spinon kinematics (issue #154 phase 1) ---
export heisenberg_spinon_dispersion,
    heisenberg_two_spinon_lower_edge, heisenberg_two_spinon_upper_edge

# --- Universality Classes ---
export Universality, CriticalExponents, GrowthExponents
export Ising2D, KPZ1D, MeanField  # backward-compatible aliases
export MinimalModel, WZWSU2       # 2D rational-CFT dispatch tags
include("universalities/Universality.jl")
include("universalities/E8.jl")
include("universalities/MeanField.jl")
include("universalities/Ising2D.jl")
include("universalities/KPZ.jl")
include("universalities/Percolation.jl")
include("universalities/Potts.jl")
include("universalities/ONModel.jl")
include("universalities/MinimalModel.jl")
include("universalities/WZW.jl")
include("universalities/CardyEntanglement.jl")

# --- Models ---
# Layout: `<class>/<Model>/<Model>.jl` (with optional sibling axis files like
# `TFIM/TFIM_thermal.jl`).  `tightbinding/<lattice-class>/` groups multiple
# tight-binding Hamiltonians by lattice type (regular = Bloch-diagonalisable;
# future: quasicrystalline, fractal, disordered).
include("models/classical/IsingSquare/IsingSquare.jl")
include("models/classical/IsingSquare/IsingSquare_thermal.jl")
include("models/classical/IsingSquare/IsingSquare_registry.jl")
include("models/classical/SixVertex/SixVertex.jl")
include("models/classical/SixVertex/SixVertex_registry.jl")
include("models/classical/IsingTriangular/IsingTriangular.jl")
include("models/classical/IsingTriangular/IsingTriangular_registry.jl")
include("models/quantum/tightbinding/regular/Honeycomb.jl")
include("models/quantum/tightbinding/regular/Kagome.jl")
include("models/quantum/tightbinding/regular/Lieb.jl")
include("models/quantum/tightbinding/regular/Triangular.jl")
include("models/quantum/TFIM/TFIM.jl")
include("models/quantum/TFIM/TFIM_dynamics.jl")
include("models/quantum/TFIM/TFIM_xx_static.jl")
include("models/quantum/TFIM/TFIM_yy.jl")
include("models/quantum/TFIM/TFIM_xx_yy_structure_factor.jl")
include("models/quantum/TFIM/TFIM_thermal.jl")
include("models/quantum/TFIM/TFIM_pbc_thermal.jl")
include("models/quantum/TFIM/TFIM_zaxis.jl")
include("models/quantum/TFIM/TFIM_local.jl")
include("models/quantum/TFIM/TFIM_sigma_x_quench.jl")
include("models/quantum/TFIM/TFIM_entanglement.jl")
include("models/quantum/TFIM/TFIM_cft_entanglement.jl")
include("models/quantum/TFIM/TFIM_infinite_dynamics.jl")
include("models/quantum/TFIM/TFIM_loschmidt.jl")
include("models/quantum/TFIM/TFIM_gge.jl")
include("models/quantum/TFIM/TFIM_registry.jl")  # populates REGISTRY for TFIM
include("models/quantum/Heisenberg/Heisenberg.jl")
include("models/quantum/Heisenberg/Heisenberg_spinon.jl")
include("models/quantum/Heisenberg/HeisenbergS1.jl")
include("models/quantum/Heisenberg/HeisenbergS1_observables.jl")
include("models/quantum/Heisenberg/HeisenbergS1_registry.jl")
include("models/quantum/AKLT/AKLT1D.jl")
include("models/quantum/AKLT/AKLT1D_registry.jl")
include("models/quantum/KitaevHoneycomb/KitaevHoneycomb.jl")
include("models/quantum/KitaevHoneycomb/KitaevHoneycomb_thermal.jl")
include("models/quantum/KitaevHoneycomb/KitaevHoneycomb_registry.jl")
include("models/quantum/XXZ/XXZ.jl")
include("models/quantum/XXZ/XXZ_bethe.jl")     # Yang-Yang single integral, used by XXZ.jl dispatch
include("models/quantum/XXZ/XXZ_thermal.jl")
include("models/quantum/XXZ/XXZ_xx_quench.jl")  # XX (Δ=0) Loschmidt rate at Infinite (issue #148)
include("models/quantum/XXZ/XXZ_registry.jl")  # populates REGISTRY for XXZ1D
include("models/quantum/Heisenberg/Heisenberg_registry.jl")  # populates REGISTRY for Heisenberg1D
include("models/quantum/Hubbard1D/Hubbard1D.jl")
include("models/quantum/Hubbard1D/Hubbard1D_registry.jl")  # populates REGISTRY for Hubbard1D
include("models/quantum/MajumdarGhosh/MajumdarGhosh.jl")
include("models/quantum/MajumdarGhosh/MajumdarGhosh_registry.jl")  # populates REGISTRY for MajumdarGhosh

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
