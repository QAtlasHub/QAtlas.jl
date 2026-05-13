# models/quantum/PpIp2DSC/PpIp2DSC_registry.jl
#
# Declarative implementation map for the 2-D p_x + i p_y chiral
# superconductor (Read-Green 2000).  Schema in src/core/registry.jl.

@register(
    PpIp2DSC,
    CentralCharge,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_ppip_2dsc.jl",
    references=["Read-Green 2000", "Kitaev 2006"],
    notes="Chiral Majorana edge CFT c=1/2 in weak-pairing topological phase (μ>0).",
)

@register(
    PpIp2DSC,
    TopologicalInvariant,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_ppip_2dsc.jl",
    references=["Read-Green 2000"],
    notes="First Chern number C=1 in weak-pairing topological phase.",
)
