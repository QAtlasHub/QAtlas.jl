# Model → model reductions (membership under a limit / special point).
#
# Included after all model types are defined.  Each row says a concrete model
# becomes another concrete model in a stated regime, which is exactly what
# legitimises a model→model `method=:delegation` fetch (the source delegates a
# quantity to the target there).  Query via `reductions(model)` /
# `reduced_from(model)`; consumed by the C4 coherence check.

@reduces HeisenbergXYZ XXZ1D regime = "two couplings equal (J^x = J^y); the XYZ chain becomes the XXZ chain"
@reduces LongRangeIsing1D TFIM regime = "α → ∞ (fast power-law decay); long-range Ising reduces to the short-range transverse-field Ising chain"
@reduces DMIHeisenberg1D Heisenberg1D regime = "uniform Dzyaloshinskii-Moriya vector gauged away by a local spin rotation; spectrum maps to the Heisenberg chain"
@reduces J1J2Heisenberg1D Heisenberg1D regime = "J₂ = 0; the J₁-J₂ chain reduces to the nearest-neighbour Heisenberg chain"
@reduces J1J2Heisenberg1D MajumdarGhosh regime = "J₂ = J₁/2 (Majumdar-Ghosh point); exact nearest-neighbour dimerized ground state"
@reduces MixedFieldIsing1D TFIM regime = "longitudinal field h_z = 0; the mixed-field Ising chain reduces to the transverse-field Ising chain"
@reduces ExtendedHubbard1D Hubbard1D regime = "nearest-neighbour interaction V = 0; the extended Hubbard chain reduces to the Hubbard model"
@reduces RandomBondIsing2D MinimalModel regime = "pure (zero-disorder) point; 2D Ising CFT, the M(4,3) minimal model, c = 1/2"
@reduces ZnClock MinimalModel regime = "self-dual critical point of the Z_n clock model maps onto a minimal-model CFT (e.g. n = 4 → M(4,3))"
