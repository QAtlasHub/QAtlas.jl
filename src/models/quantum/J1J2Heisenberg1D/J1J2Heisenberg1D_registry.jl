# models/quantum/J1J2Heisenberg1D/J1J2Heisenberg1D_registry.jl
#
# Spin-½ J₁-J₂ Heisenberg chain.  Phase-1 closed-form points only:
#   - j = J₂/J₁ = 0   → delegates to Heisenberg1D (Bethe-Hulthén 1938)
#   - j = 1/2          → delegates to MajumdarGhosh   (Majumdar-Ghosh 1969)
#   - generic j        → DomainError (DMRG; deferred to Phase 2)
#
# See `J1J2Heisenberg1D.jl` for the dispatch and `src/core/registry.jl`
# for the metadata schema.

@register(
    J1J2Heisenberg1D,
    Energy{:per_site},
    Infinite,
    method=:delegation,
    reliability=:high,
    tested_in="test/models/quantum/Heisenberg/test_j1j2_heisenberg1d.jl",
    references=["Hulthen1938", "MajumdarGhosh1969", "WhiteAffleck1996"],
    notes="Closed form at j=0 (Heisenberg1D delegate) and j=1/2 (MajumdarGhosh delegate); generic j deferred.",
)
