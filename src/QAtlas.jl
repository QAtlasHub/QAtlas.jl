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
export TricriticalIsing                                  # M(5,4) unitary minimal model, c=7/10 (Belavin-Polyakov-Zamolodchikov 1984)
export LiouvilleCFT                                      # non-compact Liouville CFT (Polyakov 1981)
export SchwingerModel                                    # 1+1-D QED (Schwinger 1962)
export ChernSimons3D                                     # SU(N)_k Chern-Simons TQFT (Witten 1989)
export SherringtonKirkpatrick                            # mean-field Ising spin glass (Sherrington-Kirkpatrick 1975)
export KagomeHeisenbergAFM                              # Kagome S=1/2 AFM Z₂ spin liquid candidate (Yan-Huse-White 2011)
export RFIM                                              # Random-Field Ising Model (Imry-Ma 1975)
export GrossNeveu                                        # 1+1-D O(2N) 4-fermion (Gross-Neveu 1974)
export XCube                                              # fracton X-cube model (Vijay-Haah-Fu 2016)
export LogarithmicCFT                                      # c=0 logarithmic CFT (polymer/percolation, #235)
export BCFT                                              # Cardy boundary CFT + Affleck-Ludwig g (#237)
export ZnParafermion                                       # Z_n parafermion CFT, Fateev-Zamolodchikov 1985 (#233)
export RandomBondIsing2D                                  # 2D ±J random-bond Ising / Edwards-Anderson (#232)
export TTbar                                                # universal TT-bar deformation (#249)
export TASEP                                                # totally asymmetric simple exclusion process (#241)
export YangLee                                              # non-unitary CFT M(5,2), c=-22/5 (#234)
export ConformalBootstrap                                   # 3D Ising bootstrap exponents (#236)
export ZnClock                                              # 2D Z_n clock model (#231)

# --- Quantum Models ---
export TFIM                                             # v0.13 concrete struct
export E8                                               # v0.13 concrete struct
export XXZ1D                                            # v0.13 new model
export KitaevHoneycomb                                  # spin-½ Kitaev honeycomb
export Kitaev1D                                         # 1D p-wave Majorana wire (Kitaev 2001)
export ToricCode                                         # Kitaev 2003 Z₂ surface code
export TightBindingSpectrum
export TightBindingChecksum, TightBindingMaxEnergy  # scalar invariants for verify()
# NOTE: `Honeycomb`, `Kagome`, `Lieb`, `Triangular` are NOT exported —
# they all conflict with Lattice2D's topology types of the same name.
# Access them as `QAtlas.Honeycomb()` / `QAtlas.Kagome()` / etc. in code
# that also uses `Lattice2D`.  `Graphene` *is* exported as the
# backward-compat top-level alias for `Honeycomb` (see src/deprecate/)
# since the name does not collide with anything in Lattice2D.
export Heisenberg1D, ExactSpectrum, GroundStateEnergyDensity
export Hubbard1D                                         # Lieb-Wu Bethe ansatz half-filling
export MajumdarGhosh                                     # spin-1/2 J1-J2 chain at MG point
export HaldaneShastry, haldane_shastry_spinon_dispersion, haldane_shastry_sound_velocity         # spin-1/2 1/r^2 inverse-square Heisenberg chain
export HeisenbergXYZ                                     # spin-½ XYZ chain (Baxter 1972, axis-aligned XXZ delegation)
export ShastrySutherland                                 # SrCu₂(BO₃)₂ 2D dimer GS (Shastry-Sutherland 1981)
export S1Heisenberg1D                                    # spin-1 (Haldane chain)
export AKLT1D                                            # spin-1 BLBQ at AKLT point
export AKLT2D                                            # 2D honeycomb AKLT VBS (#239)
export KitaevHeisenberg                                  # K-J-Γ honeycomb (α-RuCl₃, Rau-Lee-Kee 2014)
export PXP1D                                            # Rydberg-blockade chain w/ scars (#300)
export S1XXZ1D                                          # S=1 XXZ chain (#303)
export Cluster1D                                         # 1D Z₂×Z₂ SPT cluster Hamiltonian (Briegel-Raussendorf 2001)
export TightBindingV1D                                # t-V spinless fermion chain (#296)
export LongRangeIsing1D                                  # 1D power-law-decaying TFIM (#293)
export TightBinding1D                                  # 1D spinless-fermion chain (#291)
export LongRangeXY1D                                   # 1D power-law XY chain (#299)
export Compass1D                                       # 1D alternating XX/YY compass chain (#295)
export S1AnisotropicD1D                               # S=1 Heisenberg + single-ion D (#302)
export DMIHeisenberg1D                                  # spin-½ Heisenberg + Dzyaloshinskii-Moriya (#298)
export J1J2Heisenberg1D                              # spin-½ J₁-J₂ chain (#297)
export MixedFieldIsing1D                            # TFIM + longitudinal field, non-integrable (#290)
export XYh1D                                          # anisotropic XY + transverse field (LSM 1961, #292)
export ExtendedHubbard1D                              # t-U-V Hubbard chain (#294)
export FibonacciAnyons                            # non-Abelian Fibonacci anyons (#240)
export PpIp2DSC                                     # 2D p+ip chiral superconductor (Read-Green 2000, #238)
export SYK                                              # Sachdev-Ye-Kitaev (#251)

# --- Core Implementation ---
include("core/alias.jl")
include("core/type.jl")
include("core/quantities.jl")
include("core/registry.jl")
include("core/realizes.jl")  # model <-> universality-class correspondence
include("core/pfaffian.jl")
include("core/dense_ed.jl")

# --- Implementation registry public API ---
export Implementation, implementation_status, implementation_status_markdown
export references_for
export definitions, validity, canonical_scheme  # multi-definition catalog / selector
export Realization, realizes!, @realizes, realizations, realized_by  # model <-> class

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
export FractalDimension                                  # SLE_κ Hausdorff dimension (Beffara 2008, #244)
export ChiralCondensate  # massless Schwinger condensate (#246)
export GroundStateDegeneracy, TopologicalEntanglementEntropy, AnyonStatistics  # ToricCode (#162)
export CasimirEnergyCorrection                                              # CFT 1/L correction (#150)
export ConformalWeights, PrimaryFields
export StringOrderParameter
export FermiVelocity,
    LuttingerVelocity,
    SpinWaveVelocity,
    LiebRobinsonVelocity,
    MutualInformation,
    EntanglementGrowthSlope,
    CardyEntropy,
    ConformalCasimirEnergy,
    LogarithmicNegativity
export SteadyStateCurrent                                # TASEP / non-equilibrium current (#241)
export E8Spectrum
export LiebRobinsonBound  # status-axis example (:bound)
export Bound              # universal-bounds namespace: Bound{:QuantumInformation}, …
export CHSHBound          # CHSH / Bell correlator bound (:bound)
export MerminGHZBound     # Mermin 3-party Bell bound (:bound)
export ChaosBound         # MSS chaos / Lyapunov bound (:bound, Dynamics)
export BekensteinBound    # Bekenstein entropy bound (:bound, Holographic)
export QuantumSpeedLimit  # Margolus-Levitin speed limit (:bound lower, Dynamics)
export OptimalCloningFidelity  # Buzek-Hillery cloning bound (:bound upper, QuantumInformation)
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
export SpectralFormFactor  # RMT spectral form factor (#243)
include("universalities/Universality.jl")            # namespace root: Universality{C} + shared types
include("universalities/E8/E8.jl")
include("universalities/MeanField/MeanField.jl")
include("universalities/Ising2D/Ising2D.jl")
include("universalities/KPZ/KPZ.jl")
include("universalities/Percolation/Percolation.jl")
include("universalities/Potts/Potts.jl")
include("universalities/ONModel/ONModel.jl")
include("universalities/MinimalModel/MinimalModel.jl")
include("universalities/WZW/WZW.jl")
include("universalities/CardyEntanglement.jl")       # shared CFT entanglement (cross-class)

# --- Universal bounds (model-independent inequalities) ---
# A bound is NOT a universality class: it is pinned by the quantity it bounds,
# its direction (:upper/:lower), and whose bound it is.  See bounds/Bounds.jl.
include("bounds/Bounds.jl")
include("bounds/QuantumInformation/QuantumInformation.jl")
include("bounds/QuantumInformation/QuantumInformation_registry.jl")
include("bounds/Dynamics/Dynamics.jl")
include("bounds/Dynamics/Dynamics_registry.jl")
include("bounds/Holographic/Holographic.jl")
include("bounds/Holographic/Holographic_registry.jl")

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
include("models/classical/TricriticalIsing/TricriticalIsing.jl")
include("models/classical/TricriticalIsing/TricriticalIsing_registry.jl")  # populates REGISTRY for TricriticalIsing
include("models/classical/LiouvilleCFT/LiouvilleCFT.jl")
include("models/classical/LiouvilleCFT/LiouvilleCFT_registry.jl")        # populates REGISTRY for LiouvilleCFT (#248)
include("models/classical/LogarithmicCFT/LogarithmicCFT.jl")
include("models/classical/LogarithmicCFT/LogarithmicCFT_registry.jl")  # populates REGISTRY for LogarithmicCFT (#235)
include("models/classical/BCFT/BCFT.jl")
include("models/classical/BCFT/BCFT_registry.jl")  # populates REGISTRY for BCFT (#237)
include("models/classical/ZnParafermion/ZnParafermion.jl")
include("models/classical/ZnParafermion/ZnParafermion_registry.jl")  # populates REGISTRY for ZnParafermion (#233)
include("models/classical/RandomBondIsing2D/RandomBondIsing2D.jl")
include("models/classical/RandomBondIsing2D/RandomBondIsing2D_registry.jl")  # populates REGISTRY for RandomBondIsing2D (#232)
include("models/classical/TTbar/TTbar.jl")
include("models/classical/TTbar/TTbar_registry.jl")  # populates REGISTRY for TTbar (#249)
include("models/classical/TASEP/TASEP.jl")
include("models/classical/TASEP/TASEP_registry.jl")  # populates REGISTRY for TASEP (#241)
include("models/classical/YangLee/YangLee.jl")
include("models/classical/YangLee/YangLee_registry.jl")        # populates REGISTRY for YangLee (#234)
include("models/classical/ConformalBootstrap/ConformalBootstrap.jl")
include("models/classical/ConformalBootstrap/ConformalBootstrap_registry.jl")  # populates REGISTRY for ConformalBootstrap (#236)
include("models/classical/ZnClock/ZnClock.jl")
include("models/classical/ZnClock/ZnClock_registry.jl")  # populates REGISTRY for ZnClock (#231)
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
include("models/quantum/tightbinding/regular/Honeycomb_registry.jl")
include("models/quantum/tightbinding/regular/Kagome_registry.jl")
include("models/quantum/tightbinding/regular/Lieb_registry.jl")
include("models/quantum/tightbinding/regular/Triangular_registry.jl")
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
include("models/quantum/Heisenberg/Heisenberg1D_thermal_cft.jl")  # c=1 CFT low-T (#521 Path B)
include("models/quantum/Heisenberg/Heisenberg_spinon.jl")
include("models/quantum/Heisenberg/HeisenbergS1.jl")
include("models/quantum/Heisenberg/HeisenbergS1_observables.jl")
include("models/quantum/Heisenberg/HeisenbergS1_registry.jl")
include("models/quantum/AKLT/AKLT1D.jl")
include("models/quantum/AKLT/AKLT1D_registry.jl")
include("models/quantum/AKLT2D/AKLT2D.jl")
include("models/quantum/AKLT2D/AKLT2D_registry.jl")
include("models/quantum/KitaevHoneycomb/KitaevHoneycomb.jl")
include("models/quantum/KitaevHoneycomb/KitaevHoneycomb_thermal.jl")
include("models/quantum/KitaevHoneycomb/KitaevHoneycomb_registry.jl")
include("models/quantum/Kitaev1D/Kitaev1D.jl")
include("models/quantum/Kitaev1D/Kitaev1D_registry.jl")  # populates REGISTRY for Kitaev1D
include("models/quantum/XXZ/XXZ.jl")
include("models/quantum/XXZ/XXZ_bethe.jl")     # Yang-Yang single integral, used by XXZ.jl dispatch
include("models/quantum/XXZ/XXZ_thermal.jl")
include("models/quantum/XXZ/XXZ_xx_infinite.jl")
include("models/quantum/XXZ/XXZ_xx_quench.jl")
include("models/quantum/XXZ/XXZ_klumper_nlie.jl")  # Klümper QTM NLIE for critical Δ ∈ (-1,1) (issue #521)
include("models/quantum/XXZ/XXZ_registry.jl")  # populates REGISTRY for XXZ1D
include("models/quantum/Heisenberg/Heisenberg_registry.jl")  # populates REGISTRY for Heisenberg1D
include("models/quantum/Hubbard1D/Hubbard1D.jl")
include("models/quantum/Hubbard1D/Hubbard1D_jks_nlie.jl")  # JKS NLIE Stage A scaffold (#523)
include("models/quantum/Hubbard1D/Hubbard1D_registry.jl")  # populates REGISTRY for Hubbard1D
include("models/quantum/MajumdarGhosh/MajumdarGhosh.jl")
include("models/quantum/MajumdarGhosh/MajumdarGhosh_registry.jl")
include("models/quantum/HaldaneShastry/HaldaneShastry.jl")
include("models/quantum/HaldaneShastry/HaldaneShastry_thermal_cft.jl")  # c=1 CFT low-T (#524 stopgap)
include("models/quantum/HaldaneShastry/HaldaneShastry_registry.jl")  # populates REGISTRY for HaldaneShastry
include("models/quantum/ShastrySutherland/ShastrySutherland.jl")
include("models/quantum/ShastrySutherland/ShastrySutherland_registry.jl")  # populates REGISTRY for ShastrySutherland (#259)
include("models/quantum/HeisenbergXYZ/HeisenbergXYZ_xy_line.jl")  # LSM XY-line GS energy
include("models/quantum/HeisenbergXYZ/HeisenbergXYZ.jl")
include("models/quantum/HeisenbergXYZ/HeisenbergXYZ_registry.jl")          # populates REGISTRY for HeisenbergXYZ (#253)
include("models/quantum/KitaevHeisenberg/KitaevHeisenberg.jl")
include("models/quantum/KitaevHeisenberg/KitaevHeisenberg_registry.jl")    # populates REGISTRY for KitaevHeisenberg (#256)
include("models/quantum/PXP1D/PXP1D.jl")
include("models/quantum/PXP1D/PXP1D_registry.jl")  # populates REGISTRY for PXP1D (#300)
include("models/quantum/S1XXZ1D/S1XXZ1D.jl")
include("models/quantum/S1XXZ1D/S1XXZ1D_registry.jl")  # populates REGISTRY for S1XXZ1D (#303)
include("models/quantum/Cluster1D/Cluster1D.jl")
include("models/quantum/Cluster1D/Cluster1D_registry.jl")  # populates REGISTRY for Cluster1D (#301)
include("models/quantum/TightBindingV1D/TightBindingV1D.jl")
include("models/quantum/TightBindingV1D/TightBindingV1D_registry.jl")  # populates REGISTRY for TightBindingV1D (#296)
include("models/quantum/LongRangeIsing1D/LongRangeIsing1D.jl")
include("models/quantum/LongRangeIsing1D/LongRangeIsing1D_registry.jl")  # populates REGISTRY for LongRangeIsing1D (#293)
include("models/quantum/TightBinding1D/TightBinding1D.jl")
include("models/quantum/TightBinding1D/TightBinding1D_registry.jl")  # populates REGISTRY for TightBinding1D (#291)
include("models/quantum/LongRangeXY1D/LongRangeXY1D.jl")
include("models/quantum/LongRangeXY1D/LongRangeXY1D_registry.jl")  # populates REGISTRY for LongRangeXY1D (#299)
include("models/quantum/Compass1D/Compass1D.jl")
include("models/quantum/Compass1D/Compass1D_registry.jl")  # populates REGISTRY for Compass1D (#295)
include("models/quantum/S1AnisotropicD1D/S1AnisotropicD1D.jl")
include("models/quantum/S1AnisotropicD1D/S1AnisotropicD1D_registry.jl")  # populates REGISTRY for S1AnisotropicD1D (#302)
include("models/quantum/DMIHeisenberg1D/DMIHeisenberg1D.jl")
include("models/quantum/DMIHeisenberg1D/DMIHeisenberg1D_registry.jl")  # populates REGISTRY for DMIHeisenberg1D (#298)
include("models/quantum/J1J2Heisenberg1D/J1J2Heisenberg1D.jl")
include("models/quantum/J1J2Heisenberg1D/J1J2Heisenberg1D_registry.jl")  # populates REGISTRY for J1J2Heisenberg1D (#297)
include("models/quantum/MixedFieldIsing1D/MixedFieldIsing1D.jl")
include("models/quantum/MixedFieldIsing1D/MixedFieldIsing1D_registry.jl")  # populates REGISTRY for MixedFieldIsing1D (#290)
include("models/quantum/XYh1D/XYh1D.jl")
include("models/quantum/XYh1D/XYh1D_registry.jl")  # populates REGISTRY for XYh1D (#292)
include("models/quantum/ExtendedHubbard1D/ExtendedHubbard1D.jl")
include("models/quantum/ExtendedHubbard1D/ExtendedHubbard1D_registry.jl")  # populates REGISTRY for ExtendedHubbard1D (#294)
include("models/quantum/FibonacciAnyons/FibonacciAnyons.jl")
include("models/quantum/FibonacciAnyons/FibonacciAnyons_registry.jl")  # populates REGISTRY for FibonacciAnyons (#240)
include("models/quantum/PpIp2DSC/PpIp2DSC.jl")
include("models/quantum/PpIp2DSC/PpIp2DSC_registry.jl")  # populates REGISTRY for PpIp2DSC (#238)
include("models/quantum/SYK/SYK.jl")
include("models/quantum/SYK/SYK_registry.jl")  # populates REGISTRY for SYK (#251)

include("models/quantum/TFIM/TFIM_fidelity.jl")            # FidelitySusceptibility (#147)
include("models/quantum/TFIM/TFIM_quench_entanglement.jl") # VonNeumannEntropy{:quench} (#144)
include("models/quantum/ToricCode/ToricCode.jl")
include("models/quantum/ToricCode/ToricCode_registry.jl")  # populates REGISTRY for ToricCode (#162)
include("universalities/RMT/RMT.jl")                        # RMT universality class (#151)
include("universalities/Poisson/Poisson.jl")                # Poisson universality class (#151)

# --- Deprecation shims (legacy API) ---
# Loaded last so they can route into any already-registered concrete
# `fetch` method.  See src/deprecate/README.md.
include("deprecate/legacy_fetch.jl")
include("deprecate/legacy_tfim.jl")
include("deprecate/legacy_e8.jl")
include("deprecate/legacy_honeycomb.jl")
export Graphene                                         # backward-compat alias
include("deprecate/legacy_xxz.jl")

# Model <-> universality-class realizations (membership) — after all models.
include("realizes_registry.jl")

end # module QAtlas
