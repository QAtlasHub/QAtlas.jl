module QAtlas

export AbstractModel, Model
export BoundaryCondition, Infinite, PBC, OBC
export AbstractQuantity, Quantity
export fetch

# --- Classical Models ---
export IsingSquare, PartitionFunction, CriticalTemperature, SpontaneousMagnetization
export SixVertex
export IsingTriangular
export CurieWeissIsing                                   # complete-graph mean-field Ising
export IsingChain1D                                      # 1-D Ising chain (Ising 1925)
export SpinIce                                           # pyrochlore Pauling 1935
export TodaLattice                                       # 1-D Toda lattice (Toda 1967, integrable)
export SLEkappa                                          # Schramm-Loewner Evolution SLE_κ (Schramm 2000)
export TricriticalPotts3                                 # M(6,7) minimal model (Andrews-Baxter-Forrester 1984)
export LiouvilleCFT                                      # non-compact Liouville CFT (Polyakov 1981)
export SchwingerModel                                    # 1+1-D QED (Schwinger 1962)
export ChernSimons3D                                     # SU(N)_k Chern-Simons TQFT (Witten 1989)
export SherringtonKirkpatrick                            # mean-field Ising spin glass (Sherrington-Kirkpatrick 1975)
export KagomeHeisenbergAFM                              # Kagome S=1/2 AFM Z₂ spin liquid candidate (Yan-Huse-White 2011)
export RFIM                                              # Random-Field Ising Model (Imry-Ma 1975)
export GrossNeveu                                        # 1+1-D O(2N) 4-fermion (Gross-Neveu 1974)
export XCube                                              # fracton X-cube model (Vijay-Haah-Fu 2016)

# --- Quantum Models ---
export TFIM                                             # v0.13 concrete struct
export E8                                               # v0.13 concrete struct
export XXZ1D                                            # v0.13 new model
export KitaevHoneycomb                                  # spin-½ Kitaev honeycomb
export Kitaev1D                                         # 1D p-wave Majorana wire (Kitaev 2001)
export ToricCode                                         # Kitaev 2003 Z₂ surface code
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
export HeisenbergXYZ                                     # spin-½ XYZ chain (Baxter 1972, axis-aligned XXZ delegation)
export ShastrySutherland                                 # SrCu₂(BO₃)₂ 2D dimer GS (Shastry-Sutherland 1981)
export S1Heisenberg1D                                    # spin-1 (Haldane chain)
export AKLT1D                                            # spin-1 BLBQ at AKLT point
export KitaevHeisenberg                                  # K-J-Γ honeycomb (α-RuCl₃, Rau-Lee-Kee 2014)
export XYh1D                                          # anisotropic XY + transverse field (LSM 1961, #292)

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
export GroundStateDegeneracy, TopologicalEntanglementEntropy, AnyonStatistics  # ToricCode (#162)
export CasimirEnergyCorrection                                              # CFT 1/L correction (#150)
export ConformalWeights, PrimaryFields
export StringOrderParameter
export FermiVelocity, LuttingerVelocity, SpinWaveVelocity
export E8Spectrum
export TopologicalInvariant, EdgeModeEnergy           # Kitaev1D Pfaffian invariant + edge mode
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
export WignerSurmise, TracyWidom, MeanRatio  # RMT / Poisson level statistics (#151)
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
include("models/classical/CurieWeissIsing/CurieWeissIsing.jl")
include("models/classical/CurieWeissIsing/CurieWeissIsing_registry.jl")  # populates REGISTRY for CurieWeissIsing (#262)
include("models/classical/IsingChain1D/IsingChain1D.jl")
include("models/classical/IsingChain1D/IsingChain1D_registry.jl")        # populates REGISTRY for IsingChain1D (#262)
include("models/classical/SpinIce/SpinIce.jl")
include("models/classical/SpinIce/SpinIce_registry.jl")                  # populates REGISTRY for SpinIce (#257)
include("models/classical/TodaLattice/TodaLattice.jl")
include("models/classical/TodaLattice/TodaLattice_registry.jl")          # populates REGISTRY for TodaLattice (#254)
include("models/classical/SLEkappa/SLEkappa.jl")
include("models/classical/SLEkappa/SLEkappa_registry.jl")                # populates REGISTRY for SLEkappa (#244)
include("models/classical/TricriticalPotts3/TricriticalPotts3.jl")
include("models/classical/TricriticalPotts3/TricriticalPotts3_registry.jl")  # populates REGISTRY for TricriticalPotts3 (#245)
include("models/classical/LiouvilleCFT/LiouvilleCFT.jl")
include("models/classical/LiouvilleCFT/LiouvilleCFT_registry.jl")        # populates REGISTRY for LiouvilleCFT (#248)
include("models/quantum/SchwingerModel/SchwingerModel.jl")
include("models/quantum/SchwingerModel/SchwingerModel_registry.jl")      # populates REGISTRY for SchwingerModel (#246)
include("models/quantum/ChernSimons3D/ChernSimons3D.jl")
include("models/quantum/ChernSimons3D/ChernSimons3D_registry.jl")        # populates REGISTRY for ChernSimons3D (#250)
include("models/classical/SherringtonKirkpatrick/SherringtonKirkpatrick.jl")
include("models/classical/SherringtonKirkpatrick/SherringtonKirkpatrick_registry.jl")  # populates REGISTRY for SherringtonKirkpatrick (#260)
include("models/quantum/KagomeHeisenbergAFM/KagomeHeisenbergAFM.jl")
include("models/quantum/KagomeHeisenbergAFM/KagomeHeisenbergAFM_registry.jl")  # populates REGISTRY for KagomeHeisenbergAFM (#258)
include("models/quantum/GrossNeveu/GrossNeveu.jl")
include("models/quantum/GrossNeveu/GrossNeveu_registry.jl")              # populates REGISTRY for GrossNeveu (#247)
include("models/quantum/XCube/XCube.jl")
include("models/quantum/XCube/XCube_registry.jl")                        # populates REGISTRY for XCube (#252)
include("models/classical/RFIM/RFIM.jl")
include("models/classical/RFIM/RFIM_registry.jl")                        # populates REGISTRY for RFIM (#261)
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
include("models/quantum/Kitaev1D/Kitaev1D.jl")
include("models/quantum/Kitaev1D/Kitaev1D_registry.jl")  # populates REGISTRY for Kitaev1D
include("models/quantum/XXZ/XXZ.jl")
include("models/quantum/XXZ/XXZ_bethe.jl")     # Yang-Yang single integral, used by XXZ.jl dispatch
include("models/quantum/XXZ/XXZ_thermal.jl")
include("models/quantum/XXZ/XXZ_xx_infinite.jl")
include("models/quantum/XXZ/XXZ_xx_quench.jl")  # XX (Δ=0) Loschmidt rate at Infinite (issue #148)
include("models/quantum/XXZ/XXZ_registry.jl")  # populates REGISTRY for XXZ1D
include("models/quantum/Heisenberg/Heisenberg_registry.jl")  # populates REGISTRY for Heisenberg1D
include("models/quantum/Hubbard1D/Hubbard1D.jl")
include("models/quantum/Hubbard1D/Hubbard1D_registry.jl")  # populates REGISTRY for Hubbard1D
include("models/quantum/MajumdarGhosh/MajumdarGhosh.jl")
include("models/quantum/MajumdarGhosh/MajumdarGhosh_registry.jl")  # populates REGISTRY for MajumdarGhosh
include("models/quantum/ShastrySutherland/ShastrySutherland.jl")
include("models/quantum/ShastrySutherland/ShastrySutherland_registry.jl")  # populates REGISTRY for ShastrySutherland (#259)
include("models/quantum/HeisenbergXYZ/HeisenbergXYZ.jl")
include("models/quantum/HeisenbergXYZ/HeisenbergXYZ_registry.jl")          # populates REGISTRY for HeisenbergXYZ (#253)
include("models/quantum/KitaevHeisenberg/KitaevHeisenberg.jl")
include("models/quantum/KitaevHeisenberg/KitaevHeisenberg_registry.jl")    # populates REGISTRY for KitaevHeisenberg (#256)
include("models/quantum/XYh1D/XYh1D.jl")
include("models/quantum/XYh1D/XYh1D_registry.jl")  # populates REGISTRY for XYh1D (#292)

include("models/quantum/TFIM/TFIM_fidelity.jl")            # FidelitySusceptibility (#147)
include("models/quantum/TFIM/TFIM_quench_entanglement.jl") # VonNeumannEntropy{:quench} (#144)
include("models/quantum/ToricCode/ToricCode.jl")
include("models/quantum/ToricCode/ToricCode_registry.jl")  # populates REGISTRY for ToricCode (#162)
include("universalities/RMT.jl")                            # RMT universality class (#151)
include("universalities/Poisson.jl")                        # Poisson universality class (#151)

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
