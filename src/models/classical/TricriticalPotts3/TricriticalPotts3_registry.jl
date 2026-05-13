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
    references=["Andrews-Baxter-Forrester 1984", "Huse 1984"],
    notes="Delegates to MinimalModel(6,7); c = 6/7 exact (Rational).",
)
