# models/quantum/Compass1D/Compass1D_registry.jl — declarative implementation map.
#
# Compass1D Phase 1 exposes the closed-form bulk gap of the 1D
# alternating-bond compass chain (dual to a dimerised Kitaev wire):
#
#   `MassGap` at `Infinite()`  →  Δ = 2 |J_x − J_y|.

@register(
    Compass1D,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_compass1d.jl",
    references=["BrzezickiDziarmagaOles2007", "KugelKhomskii1982", "Kitaev2006"],
    notes="Δ = 2|J_x − J_y| from JW-dual dimerised Kitaev chain; Δ=0 at J_x=J_y is a first-order QPT.",
)
