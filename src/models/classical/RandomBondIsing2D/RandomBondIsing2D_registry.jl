# models/classical/RandomBondIsing2D/RandomBondIsing2D_registry.jl
#
# 2D ±J random-bond Ising — Phase 1: pure FM critical line (p = 1) only.

@register(
    RandomBondIsing2D,
    CentralCharge,
    Infinite,
    method=:delegation,
    reliability=:high,
    tested_in="test/models/classical/test_random_bond_ising_2d.jl",
    references=["Edwards-Anderson 1975", "Nishimori 1981", "Honecker-Picco-Pujol 2001"],
    notes="p=1 pure FM delegates to MinimalModel(4,3) (Ising c=1/2); Nishimori-line and multicritical point Phase 2.",
)
