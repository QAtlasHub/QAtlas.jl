# models/quantum/ShastrySutherland/ShastrySutherland_registry.jl
#
# Declarative implementation map for the Shastry–Sutherland model
# (2-D analog of the Majumdar–Ghosh dimer chain).  Schema in
# src/core/registry.jl.

@register(
    ShastrySutherland,
    Energy{:per_site},
    Infinite,
    method=:exact_dimer,
    reliability=:high,
    tested_in="test/standalone/test_shastry_sutherland.jl",
    references=["Shastry-Sutherland 1981", "Koga-Kawakami 2000"],
    notes="E0/N = -3 J'/8 on the exact dimer GS for J/J' ≤ α_c ≈ 0.675; DomainError outside.",
)
