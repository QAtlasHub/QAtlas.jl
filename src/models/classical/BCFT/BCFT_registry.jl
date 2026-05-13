# models/classical/BCFT/BCFT_registry.jl
#
# Declarative implementation map for the Cardy boundary CFT (Cardy 1989,
# Affleck-Ludwig 1991). Schema in src/core/registry.jl.

@register(
    BCFT,
    ResidualEntropy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_bcft.jl",
    references=["Cardy 1989", "Affleck-Ludwig 1991", "Friedan-Konechny 2004"],
    notes="Ising Cardy boundary entropy log g: fixed=-log(2)/2, free=σ=0, identity=epsilon=-log(2)/2; g-theorem.",
)
