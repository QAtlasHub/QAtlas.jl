# models/quantum/XYh1D/XYh1D_registry.jl
#
# Declarative implementation map for the anisotropic XY chain in a
# transverse field (Lieb-Schultz-Mattis 1961).  Schema documented in
# src/core/registry.jl.

# ── Mass Gap ──────────────────────────────────────────────────────────
@register(
    XYh1D,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_xyh1d.jl",
    references=["Lieb-Schultz-Mattis 1961", "Pfeuty 1970"],
    notes="Exact closed-form dispersion minimization for arbitrary Jx, Jy, h.",
)
@register(
    XYh1D,
    MassGap,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_xyh1d.jl",
    references=["Lieb-Schultz-Mattis 1961", "Pfeuty 1970"],
    notes="Lowest BdG quasiparticle energy from 2N×2N diagonalization.",
)

