# models/quantum/DMIHeisenberg1D/DMIHeisenberg1D_registry.jl
#
# Spin-½ Heisenberg chain with Dzyaloshinskii-Moriya interaction.
# Phase-1 closed-form point only:
#   - D = 0   → delegates to Heisenberg1D (Bethe-Hulthén 1938)
#   - D ≠ 0   → DomainError (twisted-XXZ / spiral; deferred to Phase 2,
#               Affleck-Oshikawa 1999)
#
# See `DMIHeisenberg1D.jl` for the dispatch and `src/core/registry.jl`
# for the metadata schema.

@register(
    DMIHeisenberg1D,
    Energy{:per_site},
    Infinite,
    method=:delegation,
    reliability=:high,
    tested_in="test/models/quantum/Heisenberg/test_dmi_heisenberg1d.jl",
    references=[
        "Bethe-Hulthén 1938", "Dzyaloshinskii 1958", "Moriya 1960", "Affleck-Oshikawa 1999"
    ],
    notes="D = 0 delegated to Heisenberg1D Bethe-Hulthén; D ≠ 0 spiral/twisted-XXZ deferred to Phase 2.",
)
