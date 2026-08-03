# models/quantum/ToricCode/ToricCode_registry.jl — declarative
# implementation map for the Kitaev (2003) toric code.
#
# All ToricCode quantities are closed-form analytical results:
#   - Energy{:per_site} at Infinite      = −(J_e + J_m)
#   - MassGap at Infinite                       = 2 · min(J_e, J_m)
#   - GroundStateDegeneracy at PBC              = 4^genus
#   - TopologicalEntanglementEntropy at Infinite = log 2
#
# The two anyon-statistics quantities are non-BC (the toric code's topological
# content is independent of any boundary tag); no `bc` is registered for them
# and `fetch(model, AnyonSelfStatistics; label=…)` /
# `fetch(model, AnyonMutualStatistics; anyons=…)` are the only call forms.
# They were one quantity whose returned SCHEMA changed with its `type` kwarg,
# which is why they are two now (#819).
# We register it against `Infinite` so the registry's `(model, quantity,
# bc)` triple is well-formed; the actual `fetch` method has no `bc`
# argument and always succeeds for any `type`.
#
# See `src/core/registry.jl` for the schema and `KitaevHoneycomb_registry.jl`
# for the structurally similar quantum-2D template.

@register(
    ToricCode,
    Energy{:per_site},
    Infinite,
    method=:analytic,
    cost=:closed_form,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_toric_code.jl",
    references=["Kitaev2003"],
    notes="ε₀ = −(J_e + J_m) per (vertex+plaquette) unit cell.",
)

@register(
    ToricCode,
    MassGap,
    Infinite,
    method=:analytic,
    cost=:closed_form,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_toric_code.jl",
    references=["Kitaev2003"],
    notes="Δ = 2 · min(J_e, J_m) — single-anyon excitation gap.",
)

@register(
    ToricCode,
    GroundStateDegeneracy,
    PBC,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_toric_code.jl",
    references=["Kitaev2003"],
    notes="GSD = 4^genus on a closed orientable surface; OBC has unique GS (not registered).",
)

@register(
    ToricCode,
    TopologicalEntanglementEntropy,
    Infinite,
    method=:analytic,
    cost=:closed_form,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_toric_code.jl",
    references=["KitaevPreskill2006", "LevinWen2006"],
    notes="γ = log 𝒟 = log 2 (Z₂ topological order, total quantum dim 𝒟 = 2).",
)

@register(
    ToricCode,
    AnyonSelfStatistics,
    Infinite,
    method=:analytic,
    cost=:closed_form,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_toric_code.jl",
    references=["Kitaev2003", "NayakSimonSternFreedmanDasSarma2008"],
    notes="Self-statistics of {1, e, m, ε}: bosons except ε, which is a fermion " *
          "(self-phase π); all quantum dimensions 1 (Abelian).",
)
@register(
    ToricCode,
    AnyonMutualStatistics,
    Infinite,
    method=:analytic,
    cost=:closed_form,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_toric_code.jl",
    references=["Kitaev2003", "NayakSimonSternFreedmanDasSarma2008"],
    notes="e/m mutual braiding phase π (Z₂ mutual semion); every other pair " *
          "braids trivially.",
)
