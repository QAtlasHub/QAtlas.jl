# models/classical/ConformalBootstrap/ConformalBootstrap_registry.jl
#
# Declarative implementation map for the 3D Ising conformal bootstrap.
# Methodology: Kos-Poland-Simmons-Duffin 2014 (mixed-correlator
# bootstrap).  Precise Δ_σ, Δ_ε values: KPSD-Vichi 2016 "Precision
# Islands" (arXiv:1603.04436).  Schema in src/core/registry.jl.

@register(
    ConformalBootstrap,
    ConformalWeights,
    Infinite,
    method=:bootstrap_reference,
    reliability=:high,
    tested_in="test/universalities/test_conformal_bootstrap.jl",
    references=[
        "KosPolandSimmonsDuffin2014", "Kos2016", "SimmonsDuffin2017", "Reehorst2021"
    ],
    notes="3D Ising critical exponents Δ_σ=0.5181489(10), Δ_ε=1.412625(10) — KPSD-Vichi 2016 Precision Islands; KPSD 2014 methodology, Simmons-Duffin 2017 cross-check. Higher operators Phase 2.",
)
