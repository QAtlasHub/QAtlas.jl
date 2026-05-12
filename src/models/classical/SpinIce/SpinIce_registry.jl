# models/classical/SpinIce/SpinIce_registry.jl
#
# Declarative implementation map for the classical spin-ice model on
# the pyrochlore lattice (Pauling 1935).  Schema in src/core/registry.jl.

@register(
    SpinIce,
    ResidualEntropy,
    Infinite,
    method=:analytic,
    reliability=:medium,
    tested_in="test/standalone/test_spin_ice.jl",
    references=["Pauling 1935"],
    notes="Pauling 1935 mean-tetrahedron closed form S/N = (1/2) log(3/2); a few percent below Nagle 1966.",
)
