# Model ↔ universality-class realizations (membership).
#
# Included after all model types are defined.  Each row says a concrete model
# flows to a universality class in a stated regime; query via
# `realizations(model)` / `realized_by(class)` / `realized_class(instance)`.
#
# `at` is a predicate marking the *critical locus* (point / line / surface) in
# parameter space where the model realizes the class; `example` is a
# representative critical instance on it. For one model the loci must be
# mutually exclusive (a critical point belongs to exactly one class) — verified
# by the C8 coherence check. Classical T_c-criticality (temperature is a fetch
# kwarg, not a struct field) carries only `regime` for now.

@realizes TFIM :Ising regime = "quantum critical point h = J; (1+1)D Ising CFT, c = 1/2" at = (
    m -> m.h == m.J
) example = TFIM(; J=1.0, h=1.0)

@realizes XXZ1D :XY regime = "critical line -1 < Δ < 1; Luttinger liquid (free boson), c = 1" at = (
    m -> -1 < m.Δ < 1
) example = XXZ1D(; Δ=0.0)
@realizes XXZ1D :Heisenberg regime = "isotropic point Δ = 1; SU(2)_1 WZW, c = 1" at = (
    m -> m.Δ == 1
) example = XXZ1D(; Δ=1.0)

@realizes Heisenberg1D :Heisenberg regime = "isotropic AFM point; SU(2)_1 WZW, c = 1" at = (
    m -> true
) example = Heisenberg1D()
@realizes HaldaneShastry :Heisenberg regime = "ground state of the 1/r² inverse-square chain; SU(2)_1 WZW, c = 1" at = (
    m -> true
) example = HaldaneShastry()

# Classical T_c-criticality: regime only (temperature is a fetch kwarg).
@realizes IsingSquare :Ising regime = "2D classical Ising at T_c; 2D Ising universality, c = 1/2"
@realizes IsingTriangular :Ising regime = "ferromagnetic triangular-lattice Ising at T_c; 2D Ising universality, c = 1/2"
@realizes CurieWeissIsing :MeanField regime = "complete-graph (infinite-range) Ising; mean-field critical exponents"
@realizes TASEP :KPZ regime = "current fluctuations of the 1D exclusion process; KPZ universality"
