# models/quantum/LongRangeIsing1D/LongRangeIsing1D_registry.jl
#
# Declarative implementation map for the 1D long-range (power-law)
# transverse-field Ising chain.  Schema documented in
# src/core/registry.jl.

@register(
    LongRangeIsing1D,
    MassGap,
    Infinite,
    method=:delegation,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_long_range_ising1d.jl",
    references=["Pfeuty1970", "KoffelLewensteinTagliacozzo2012", "GongFossFeig2014"],
    notes="α=Inf delegated to TFIM (Δ = 2|h-J|); finite α requires DMRG, deferred to Phase 2.",
)
