# models/quantum/S1XXZ1D/S1XXZ1D_registry.jl
#
# Declarative implementation map for S1XXZ1D (spin-1 XXZ chain).
# Schema documented in src/core/registry.jl.

@register(
    S1XXZ1D,
    MassGap,
    Infinite,
    method=:s1_heisenberg_delegation,
    reliability=:medium,
    tested_in="test/models/quantum/Heisenberg/test_s1_xxz1d.jl",
    references=["WhiteHuse1993", "Schulz1986", "TzengYang2008"],
    notes="Δ=1 delegate to S1Heisenberg1D Haldane gap; Δ≠1 (XY1/large-Δ Néel) Phase 2.",
)

@register(
    S1XXZ1D,
    Energy{:per_site},
    Infinite,
    method=:s1_heisenberg_delegation,
    reliability=:medium,
    tested_in="test/models/quantum/Heisenberg/test_s1_xxz1d.jl",
    references=["WhiteHuse1993", "Schulz1986", "TzengYang2008"],
    notes="Δ=1 delegate to S1Heisenberg1D Energy{:per_site}; Δ≠1 (XY1/large-Δ Néel) Phase 2.",
)
