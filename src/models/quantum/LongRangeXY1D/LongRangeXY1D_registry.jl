# models/quantum/LongRangeXY1D/LongRangeXY1D_registry.jl
#
# Declarative implementation map for the 1D long-range XY chain in a
# transverse field.  Schema documented in src/core/registry.jl.

@register(
    LongRangeXY1D,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_long_range_xy1d.jl",
    references=["Lieb-Schultz-Mattis 1961", "Pfeuty 1970", "Maghrebi-Gong-Gorshkov 2017"],
    notes="α=Inf NN XX limit: Δ = 2·max(0, |h|-2J); finite α deferred to Phase 2.",
)
