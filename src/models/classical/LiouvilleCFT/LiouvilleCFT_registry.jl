# models/classical/LiouvilleCFT/LiouvilleCFT_registry.jl
#
# Declarative implementation map for the Liouville CFT (Polyakov 1981).
# Schema in src/core/registry.jl.

@register(
    LiouvilleCFT,
    CentralCharge,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_liouville_cft.jl",
    references=["Polyakov1981"],
    notes="c = 1 + 6(b + 1/b)²; invariant under b↔1/b self-duality.",
)

@register(
    LiouvilleCFT,
    ConformalWeights,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_liouville_cft.jl",
    references=["Polyakov1981", "ZamolodchikovZamolodchikov1996"],
    notes="Vertex operator Δ_α = α(Q − α); reflection α↔Q−α symmetry; degenerate at α=b, 1/b.",
)
