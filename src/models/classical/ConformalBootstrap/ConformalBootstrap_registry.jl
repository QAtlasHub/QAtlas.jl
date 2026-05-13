# models/classical/ConformalBootstrap/ConformalBootstrap_registry.jl
#
# Declarative implementation map for the 3D Ising conformal bootstrap
# (Kos-Poland-Simmons-Duffin 2014).  Schema in src/core/registry.jl.

@register(
    ConformalBootstrap,
    ConformalWeights,
    Infinite,
    method=:bootstrap_reference,
    reliability=:high,
    tested_in="test/universalities/test_conformal_bootstrap.jl",
    references=[
        "Kos-Poland-Simmons-Duffin 2014", "Simmons-Duffin 2017", "Reehorst et al 2021"
    ],
    notes="3D Ising critical exponents Δ_σ=0.51815, Δ_ε=1.41263; higher operators Phase 2.",
)
