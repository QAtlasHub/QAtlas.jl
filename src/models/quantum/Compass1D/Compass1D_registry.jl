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
    references=[
        "Brzezicki-Dziarmaga-Oles PRB 75, 134415 (2007)",
        "Kugel-Khomskii Sov. Phys. Usp. 25, 231 (1982)",
        "Kitaev Annals of Physics 321, 2 (2006)",
    ],
    notes="Δ = 2|J_x − J_y| from JW-dual dimerised Kitaev chain; Δ=0 at J_x=J_y is a first-order QPT.",
)
