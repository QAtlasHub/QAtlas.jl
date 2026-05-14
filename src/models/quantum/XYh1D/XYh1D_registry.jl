# models/quantum/XYh1D/XYh1D_registry.jl
#
# Declarative implementation map for the anisotropic XY chain in a
# transverse field (Lieb-Schultz-Mattis 1961).  Schema documented in
# src/core/registry.jl.

@register(
    XYh1D,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_xyh1d.jl",
    references=["Lieb-Schultz-Mattis 1961", "Pfeuty 1970"],
    notes="Phase 1: isotropic XX limit (Jx = Jy) only; MassGap = 2·max(0, |h| − 2J). Anisotropic Jx ≠ Jy deferred to Phase 2.",
)

@register(
    XYh1D,
    Energy{:per_site},
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_xyh1d.jl",
    references=["Lieb-Schultz-Mattis 1961", "Pfeuty 1970"],
    notes="Phase 1: isotropic XX limit (Jx = Jy) only; closed-form free-fermion ground-state energy density E/N = -h + (2h/π)·arccos(h/2J) - (4J/π)·√(1-(h/2J)²) for |h|<2J, E/N = -|h| for |h|≥2J. h=0 gives E/N = -4J/π. Anisotropic Jx ≠ Jy deferred to Phase 2.",
)
