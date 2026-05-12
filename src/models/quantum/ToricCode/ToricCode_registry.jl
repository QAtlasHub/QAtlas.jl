# models/quantum/ToricCode/ToricCode_registry.jl — declarative
# implementation map for the Kitaev (2003) toric code.
#
# All ToricCode quantities are closed-form analytical results:
#   - GroundStateEnergyDensity at Infinite      = −(J_e + J_m)
#   - MassGap at Infinite                       = 2 · min(J_e, J_m)
#   - GroundStateDegeneracy at PBC              = 4^genus
#   - TopologicalEntanglementEntropy at Infinite = log 2
#
# `AnyonStatistics` is a non-BC quantity (the toric code's topological
# content is independent of any boundary tag); no `bc` is registered for
# it and `fetch(model, AnyonStatistics; type=…)` is the only call form.
# We register it against `Infinite` so the registry's `(model, quantity,
# bc)` triple is well-formed; the actual `fetch` method has no `bc`
# argument and always succeeds for any `type`.
#
# See `src/core/registry.jl` for the schema and `KitaevHoneycomb_registry.jl`
# for the structurally similar quantum-2D template.

@register(
    ToricCode,
    GroundStateEnergyDensity,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_toric_code.jl",
    references=["Kitaev 2003"],
    notes="ε₀ = −(J_e + J_m) per (vertex+plaquette) unit cell.",
)

@register(
    ToricCode,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_toric_code.jl",
    references=["Kitaev 2003"],
    notes="Δ = 2 · min(J_e, J_m) — single-anyon excitation gap.",
)

@register(
    ToricCode,
    GroundStateDegeneracy,
    PBC,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_toric_code.jl",
    references=["Kitaev 2003"],
    notes="GSD = 4^genus on a closed orientable surface; OBC has unique GS (not registered).",
)

@register(
    ToricCode,
    TopologicalEntanglementEntropy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_toric_code.jl",
    references=["Kitaev-Preskill 2006", "Levin-Wen 2006"],
    notes="γ = log 𝒟 = log 2 (Z₂ topological order, total quantum dim 𝒟 = 2).",
)

@register(
    ToricCode,
    AnyonStatistics,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_toric_code.jl",
    references=["Kitaev 2003", "Nayak-Simon-Stern-Freedman-Das Sarma 2008"],
    notes="Topological data for {1, e, m, ε} and e/m mutual braiding (phase π).",
)
