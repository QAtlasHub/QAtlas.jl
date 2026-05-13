# models/classical/YangLee/YangLee_registry.jl
#
# Declarative implementation map for the Yang-Lee CFT — non-unitary
# Virasoro minimal model M(5, 2), c = -22/5 (Cardy 1985; Yang-Lee 1952).
# Schema in src/core/registry.jl.

@register(
    YangLee,
    CentralCharge,
    Infinite,
    method=:minimal_model_delegation,
    reliability=:high,
    tested_in="test/models/classical/test_yang_lee.jl",
    references=["Cardy 1985", "Yang-Lee 1952"],
    notes="Non-unitary minimal model M(5,2); c = -22/5 exact (Rational), delegated to MinimalModel.",
)

@register(
    YangLee,
    ConformalWeights,
    Infinite,
    method=:minimal_model_delegation,
    reliability=:high,
    tested_in="test/models/classical/test_yang_lee.jl",
    references=["Cardy 1985"],
    notes="Kac formula h_{r,s} via MinimalModel(5,2); famous Yang-Lee primary h_{1,2} = -1/5.",
)
