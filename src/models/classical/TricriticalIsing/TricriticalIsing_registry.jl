# models/classical/TricriticalIsing/TricriticalIsing_registry.jl
#
# Declarative implementation map for the tricritical Ising CFT
# (M(5, 4) unitary minimal model).  Schema in src/core/registry.jl.

@register(
    TricriticalIsing,
    CentralCharge,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_tricritical_ising.jl",
    references=["BelavinPolyakovZamolodchikov1984", "FriedanQiuShenker1984"],
    notes="M(5,4) unitary minimal model, c=7/10; delegated to MinimalModel.",
)

@register(
    TricriticalIsing,
    ConformalWeights,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_tricritical_ising.jl",
    references=["BelavinPolyakovZamolodchikov1984"],
    notes="Kac formula h_{r,s} via MinimalModel(5,4); σ at 3/80, ε at 1/10.",
)

@register(
    TricriticalIsing,
    PrimaryFields,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_tricritical_ising.jl",
    references=["BelavinPolyakovZamolodchikov1984"],
    notes="6 primaries modulo Kac symmetry.",
)
