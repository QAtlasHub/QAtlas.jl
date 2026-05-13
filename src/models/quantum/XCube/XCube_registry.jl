# models/quantum/XCube/XCube_registry.jl
#
# Fracton X-cube ground-state degeneracy at PBC.

@register(
    XCube,
    GroundStateDegeneracy,
    PBC,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_x_cube.jl",
    references=["Vijay-Haah-Fu 2016", "Slagle-Kim 2017"],
    notes="Subextensive log_2 GSD = 2(Lx+Ly+Lz) - 3 on closed cubic torus.",
)
