# models/classical/TricriticalPotts3/TricriticalPotts3_registry.jl
#
# Declarative implementation map for the tricritical 3-state Potts
# CFT (M(6, 7) minimal model).  Schema in src/core/registry.jl.

@register(
    TricriticalPotts3,
    CentralCharge,
    Infinite,
    method=:minimal_model_delegation,
    reliability=:high,
    tested_in="test/standalone/test_tricritical_potts3.jl",
    references=["AndrewsBaxterForrester1984", "Huse1984"],
    notes="Delegates to MinimalModel(6,7); c = 6/7 exact (Rational).",
)

@register(
    TricriticalPotts3,
    ConformalWeights,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_tricritical_potts3.jl",
    references=["BelavinPolyakovZamolodchikov1984", "AndrewsBaxterForrester1984"],
    notes="Delegated to MinimalModel(7, 6) Kac formula; r ∈ [1, 5], s ∈ [1, 6].",
)

@register(
    TricriticalPotts3,
    PrimaryFields,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_tricritical_potts3.jl",
    references=["BelavinPolyakovZamolodchikov1984"],
    notes="Delegated to MinimalModel(7, 6); 15 independent primaries modulo Kac symmetry.",
)
