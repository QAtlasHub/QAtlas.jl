# models/classical/LogarithmicCFT/LogarithmicCFT_registry.jl
#
# Declarative implementation map for the c = 0 logarithmic CFT
# (polymer / percolation universality; issue #235).
# Schema in src/core/registry.jl.

@register(
    LogarithmicCFT,
    CentralCharge,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_logarithmic_cft.jl",
    references=[
        "Saleur 1992",
        "Cardy 2001",
        "Pearce-Rasmussen-Zuber 2006",
        "Vasseur-Jacobsen-Saleur 2011",
    ],
    notes="c=0 logarithmic CFT (polymer/percolation universality); indecomposable rep structure Phase 2."
)
