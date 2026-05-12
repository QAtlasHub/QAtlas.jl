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
    references=["Yang-Yang 1966", "Baxter 1972"],
    notes="Delegates to XXZ1D(J=Jx, Δ=Jz/Jx) when Jx=Jy; general (Jx≠Jy) raises DomainError (Baxter elliptic deferred to Phase 2).",
)
