# models/quantum/AKLT2D/AKLT2D_registry.jl —
# declarative implementation map for AKLT2D
# (spin-3/2 AKLT VBS on the hexagonal / honeycomb lattice, refs #239).
#
# Phase 1 exposes only the frustration-free closed-form GS energy
# density. The Hamiltonian is a sum of non-negative spin-3 projectors
# annihilating the VBS state (AKLT 1988; Verstraete-Cirac 2004 PEPS),
# so the result is exact and J-independent — `reliability=:high`.

@register(AKLT2D, Energy{:per_site}, Infinite,
    method=:analytic, reliability=:high,
    tested_in="test/models/quantum/misc/test_aklt2d.jl",
    references=["Affleck-Kennedy-Lieb-Tasaki 1988", "Verstraete-Cirac 2004"],
    notes="Honeycomb-lattice spin-3/2 AKLT; frustration-free → exact zero GS energy density.")
