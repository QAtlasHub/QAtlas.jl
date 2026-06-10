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
    m -> isapprox(m.h, m.J; atol=1e-10)
) example = TFIM(; J=1.0, h=1.0)
@realizes TFIM :IsingSDRG regime = "strong-disorder limit / infinite-randomness fixed point (IRFP) under random bond/field couplings" references = ["Fisher1992", "Fisher1995", "RefaelMoore2004"]

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

@realizes Kitaev1D :Ising regime = "critical line |μ| = 2|t|; (1+1)D Ising CFT, c = 1/2" at = (
    m -> isapprox(abs(m.μ), 2*abs(m.t); atol=1e-10)
) example = Kitaev1D(; μ=2.0, t=1.0, Δ=1.0) references = ["Kitaev2001"]

@realizes ZnClock :Ising regime = "n = 2 clock model; 2D classical Ising CFT, c = 1/2" at = (
    m -> m.n == 2
) example = ZnClock(; n=2) references = ["JoseKadanoffKirkpatrickNelson1977", "ElitzurPearsonShigemitsu1979"]
@realizes ZnClock :Potts3 regime = "n = 3 clock model; 3-state Potts CFT, c = 4/5" at = (
    m -> m.n == 3
) example = ZnClock(; n=3) references = ["JoseKadanoffKirkpatrickNelson1977", "ElitzurPearsonShigemitsu1979"]

@realizes ZnParafermion :Ising regime = "n = 2 parafermions; (1+1)D Ising CFT, c = 1/2" at = (
    m -> m.n == 2
) example = ZnParafermion(; n=2) references = ["FateevZamolodchikov1985"]
@realizes ZnParafermion :Potts3 regime = "n = 3 parafermions; 3-state Potts CFT, c = 4/5" at = (
    m -> m.n == 3
) example = ZnParafermion(; n=3) references = ["FateevZamolodchikov1985"]
@realizes ZnParafermion :Potts4 regime = "n = 4 parafermions; compact free boson (c = 1)" at = (
    m -> m.n == 4
) example = ZnParafermion(; n=4) references = ["FateevZamolodchikov1985"]

@realizes SixVertex :XY regime = "disordered phase |Δ| < 1; compact free boson (Luttinger liquid / XY class), c = 1" at = (
    m -> -1 < _six_vertex_delta(m.a, m.b, m.c) < 1
) example = SixVertex(; a=1.0, b=1.0, c=1.0) references = ["Lieb1967a", "Sutherland1967"]

@realizes DimerLattice :XY regime = "close-packed dimer model; height representation is a c = 1 compact free boson (XY class)" references = ["Kasteleyn1961", "Fisher1961"]

@realizes TricriticalIsing :TricriticalIsing regime = "tricritical point of vacancy-extended Ising; M(5, 4) minimal model, c = 7/10" references = ["BelavinPolyakovZamolodchikov1984", "FriedanQiuShenker1984"]
@realizes TricriticalPotts3 :TricriticalPotts3 regime = "dilute q = 3 Potts model at criticality; M(6, 7) minimal model, c = 6/7" references = ["AndrewsBaxterForrester1984", "Huse1984"]

@realizes SSH :XY regime = "critical line |v| = |w|; (1+1)D free Dirac fermion / XY class, c = 1" at = (
    m -> isapprox(abs(m.v), abs(m.w); atol=1e-10)
) example = SSH(; v=1.0, w=1.0) references = ["SSH1979"]

@realizes YangLee :LeeYang regime = "Lee-Yang edge singularity; non-unitary minimal model M(5, 2), c = -22/5" references = ["Cardy1985"]

# ─── UniversalityClass registrations (Edges) ───────────────────────────
@register(
    TFIM,
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Emergent universality class: :Ising (clean TFIM), :IsingSDRG (IRFP / random-bond limit)"
)
@register(
    XXZ1D,
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Emergent universality class: :XY, :Heisenberg"
)
@register(
    Heisenberg1D,
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Emergent universality class: :Heisenberg"
)
@register(
    HaldaneShastry,
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Emergent universality class: :Heisenberg"
)
@register(
    IsingSquare,
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Emergent universality class: :Ising"
)
@register(
    IsingTriangular,
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Emergent universality class: :Ising"
)
@register(
    CurieWeissIsing,
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Emergent universality class: :MeanField"
)
@register(
    TASEP,
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Emergent universality class: :KPZ"
)
@register(
    Kitaev1D,
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Emergent universality class: :Ising"
)
@register(
    ZnClock,
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Emergent universality class: :Ising, :Potts3"
)
@register(
    ZnParafermion,
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Emergent universality class: :Ising, :Potts3, :Potts4"
)
@register(
    SixVertex,
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Emergent universality class: :XY"
)
@register(
    DimerLattice,
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Emergent universality class: :XY"
)
@register(
    TricriticalIsing,
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Emergent universality class: :TricriticalIsing"
)
@register(
    TricriticalPotts3,
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Emergent universality class: :TricriticalPotts3"
)
@register(
    SSH,
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Emergent universality class: :XY"
)
@register(
    YangLee,
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Emergent universality class: :LeeYang"
)

# Declarative registration for Universality{C} -> UniversalityClass
@register(
    Universality{:Ising},
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Universality class identity: :Ising"
)
@register(
    Universality{:XY},
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Universality class identity: :XY"
)
@register(
    Universality{:Heisenberg},
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Universality class identity: :Heisenberg"
)
@register(
    Universality{:MeanField},
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Universality class identity: :MeanField"
)
@register(
    Universality{:KPZ},
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Universality class identity: :KPZ"
)
@register(
    Universality{:Potts3},
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Universality class identity: :Potts3"
)
@register(
    Universality{:Potts4},
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Universality class identity: :Potts4"
)
@register(
    Universality{:TricriticalIsing},
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Universality class identity: :TricriticalIsing"
)
@register(
    Universality{:TricriticalPotts3},
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Universality class identity: :TricriticalPotts3"
)
@register(
    Universality{:LeeYang},
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Universality class identity: :LeeYang"
)
@register(
    Universality{:IsingSDRG},
    UniversalityClass,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_universality_class.jl",
    notes="Universality class identity: :IsingSDRG"
)
