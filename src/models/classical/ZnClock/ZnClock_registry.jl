# models/classical/ZnClock/ZnClock_registry.jl
#
# Declarative implementation map for the 2-D Z_n clock model.
# Schema in src/core/registry.jl.

@register(
    ZnClock,
    CentralCharge,
    Infinite,
    method=:delegation,
    reliability=:high,
    tested_in="test/universalities/test_zn_clock.jl",
    references=["JoseKadanoffKirkpatrickNelson1977", "ElitzurPearsonShigemitsu1979"],
    notes="n=2 (Ising c=1/2) and n=3 (Potts c=4/5) delegate to MinimalModel; n≥4 (Ashkin-Teller/BKT line) Phase 2.",
)
