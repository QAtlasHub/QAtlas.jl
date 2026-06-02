# models/quantum/ExtendedHubbard1D/ExtendedHubbard1D_registry.jl
#
# Declarative implementation map for the 1D t-U-V Hubbard chain.
#
# Phase 1 only registers the V = 0 delegate to `Hubbard1D` (ChargeGap
# at half filling, Lieb-Wu 1968 integral).  The V ≠ 0 phase diagram
# (CDW / SDW / BOW / phase separation, Voit 1995; Nakamura 2000) is
# deferred to Phase 2.

@register(
    ExtendedHubbard1D,
    ChargeGap,
    Infinite,
    method=:delegation,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_extended_hubbard1d.jl",
    references=["LiebWu1968", "Voit1995", "Nakamura2000"],
    notes="V=0 delegates to Hubbard1D at half filling (Lieb-Wu integral); V≠0 (CDW/SDW/BOW phase diagram) deferred to Phase 2.",
)
