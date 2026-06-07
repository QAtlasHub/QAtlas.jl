# Model ↔ universality-class realizations (membership).
#
# Included after all model types are defined.  Each row says a concrete model
# flows to a universality class in a stated regime; query via
# `realizations(model)` / `realized_by(class)`.

@realizes TFIM :Ising regime = "quantum critical point h = J; (1+1)D Ising CFT, c = 1/2"
@realizes XXZ1D :XY regime = "critical line -1 < Δ ≤ 1; Luttinger liquid (free boson), c = 1"
@realizes Heisenberg1D :Heisenberg regime = "isotropic AFM point; SU(2)_1 WZW, c = 1"
@realizes IsingSquare :Ising regime = "2D classical Ising at T_c; 2D Ising universality, c = 1/2"
@realizes CurieWeissIsing :MeanField regime = "complete-graph (infinite-range) Ising; mean-field critical exponents"
@realizes TASEP :KPZ regime = "current fluctuations of the 1D exclusion process; KPZ universality"
