# models/classical/TASEP/TASEP_registry.jl
#
# Declarative implementation map for TASEP (Derrida-Evans-Hakim-Pasquier
# 1993; Derrida-Lebowitz 1998).  Schema in src/core/registry.jl.

@register(
    TASEP,
    SteadyStateCurrent,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_tasep.jl",
    references=[
        "KardarParisiZhang1986", "DerridaEvansHakimPasquier1993", "DerridaLebowitz1998"
    ],
    notes="Mean-field TASEP steady-state current j(ρ) = p ρ (1−ρ); KPZ-class non-equilibrium.",
)
