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
    references=["YangYang1966", "LiebSchultzMattis1961", "Baxter1972"],
    notes="Delegates to XXZ1D(J=Jx, Δ=Jz/Jx) when Jx=Jy; XY anisotropic line (Jz=0) via Lieb-Schultz-Mattis 1961 closed form; generic XYZ (Jx≠Jy, Jz≠0) deferred to Baxter elliptic Phase 3.",
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

@register(
    HeisenbergXYZ,
    GroundStateEnergyDensity,
    Infinite,
    method=:closed_form,
    reliability=:high,
    tested_in="test/models/quantum/HeisenbergXYZ/test_heisenberg_xyz_gs.jl",
    references=["LiebSchultzMattis1961", "Baxter1972"],
    notes="XY anisotropic line (Jz=0) via Lieb-Schultz-Mattis 1961 free-fermion " *
          "closed form; axial XXZ case (Jx=Jy) delegated to XXZ1D Energy(:per_site). " *
          "Generic XYZ (Jx!=Jy, Jz!=0) deferred to Baxter elliptic Phase 3.",
)
