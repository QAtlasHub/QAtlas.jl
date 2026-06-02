# models/quantum/HeisenbergXYZ/HeisenbergXYZ_registry.jl
#
# Declarative implementation map for the spin-½ XYZ chain
# (most general 1-D nearest-neighbour integrable spin model).
# Schema documented in src/core/registry.jl.

@register(
    HeisenbergXYZ,
    Energy{:per_site},
    Infinite,
    method=:xxz_delegation,
    reliability=:high,
    tested_in="test/standalone/test_heisenberg_xyz.jl",
    references=["YangYang1966", "Baxter1972"],
    notes="Delegates to XXZ1D(J=Jx, Δ=Jz/Jx) when Jx=Jy; general (Jx≠Jy) raises DomainError (Baxter elliptic deferred to Phase 2).",
)

@register(
    HeisenbergXYZ,
    LuttingerParameter,
    Infinite,
    method=:delegation,
    reliability=:high,
    tested_in="test/models/quantum/XXZ/test_heisenberg_xyz.jl",
    references=["LutherPeschel1975", "Baxter1972"],
    notes="K=1/2 at isotropic Jx=Jy=Jz; delegates to XXZ1D(Δ=1) (same target as Heisenberg1D path). Generic XYZ Phase 2 (Baxter elliptic).",
)
