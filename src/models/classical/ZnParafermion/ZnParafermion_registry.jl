# models/classical/ZnParafermion/ZnParafermion_registry.jl
#
# Declarative implementation map for the Z_n parafermion CFT
# (Fateev-Zamolodchikov 1985). Schema in src/core/registry.jl.

@register(
    ZnParafermion,
    CentralCharge,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_zn_parafermion.jl",
    references=["Fateev-Zamolodchikov 1985"],
    notes="c = 2(n-1)/(n+2); SU(2)_n / U(1) coset CFT; n=2 = Ising, n=3 = Potts.",
)
