# models/quantum/SYK/SYK_registry.jl
#
# Declarative implementation map for the Sachdev-Ye-Kitaev model.
# Schema documented in src/core/registry.jl.

@register(
    SYK,
    ConformalWeights,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_syk.jl",
    references=["SachdevYe1993", "Kitaev2015", "MaldacenaStanford2016"],
    notes="Large-N IR Majorana conformal dimension Δ_ψ = 1/q; composite-operator dimensions Phase 2.",
)
