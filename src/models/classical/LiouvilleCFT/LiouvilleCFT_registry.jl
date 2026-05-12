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
    references=["Polyakov 1981"],
    notes="c = 1 + 6(b + 1/b)²; invariant under b↔1/b self-duality.",
)
