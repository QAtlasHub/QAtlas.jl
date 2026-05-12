# models/quantum/SchwingerModel/SchwingerModel_registry.jl
#
# Declarative implementation map for the 1+1-D Schwinger model
# (Schwinger 1962).  Schema in src/core/registry.jl.

@register(
    SchwingerModel,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_schwinger_model.jl",
    references=["Schwinger 1962"],
    notes="Massless Schwinger m_γ = e/√π via abelian bosonisation; m ≠ 0 raises DomainError (sine-Gordon dual, Phase 2).",
)
