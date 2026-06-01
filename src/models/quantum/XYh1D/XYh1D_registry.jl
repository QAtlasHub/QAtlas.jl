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

# ── Energy ────────────────────────────────────────────────────────────
@register(
    XYh1D,
    Energy{:per_site},
    Infinite,
    method=:quadgk,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_xyh1d.jl",
    references=["Lieb-Schultz-Mattis 1961", "Pfeuty 1970"],
    notes="Ground-state and finite-T per-site energy via Gauss-Kronrod integration over dispersion.",
)
@register(
    XYh1D,
    Energy{:total},
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_xyh1d.jl",
    references=["Lieb-Schultz-Mattis 1961", "Pfeuty 1970"],
    notes="Total ground-state and finite-T energy via BdG spectrum sum.",
)

# ── Scalar Thermodynamic Potentials ───────────────────────────────────
for QTy in (FreeEnergy, ThermalEntropy, SpecificHeat)
    register!(
        XYh1D,
        QTy,
        Infinite;
        method=:quadgk,
        reliability=:high,
        tested_in="test/models/quantum/misc/test_xyh1d.jl",
        references=["Lieb-Schultz-Mattis 1961", "Pfeuty 1970"],
        notes="Per-site thermal quantity via QuadGK over dispersion.",
    )
    register!(
        XYh1D,
        QTy,
        OBC;
        method=:bdg,
        reliability=:high,
        tested_in="test/models/quantum/misc/test_xyh1d.jl",
        references=["Lieb-Schultz-Mattis 1961", "Pfeuty 1970"],
        notes="Per-site thermal quantity from OBC BdG spectrum sum.",
    )
end

# ── Magnetization & Susceptibility (Phase 2, #292) ─────────────────────
for QTy in (MagnetizationZ, SusceptibilityZZ)
    register!(
        XYh1D,
        QTy,
        Infinite;
        method=:quadgk,
        reliability=:high,
        tested_in="test/models/quantum/misc/test_xyh1d.jl",
        references=["Lieb-Schultz-Mattis 1961", "Pfeuty 1970"],
        notes="Per-site value via QuadGK over dispersion.",
    )
    register!(
        XYh1D,
        QTy,
        OBC;
        method=:bdg,
        reliability=:high,
        tested_in="test/models/quantum/misc/test_xyh1d.jl",
        references=["Lieb-Schultz-Mattis 1961", "Pfeuty 1970"],
        notes="Per-site value from OBC BdG spectrum / Majorana covariance.",
    )
end
register!(
    XYh1D,
    MagnetizationZLocal,
    OBC;
    method=:bdg,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_xyh1d.jl",
    references=["Lieb-Schultz-Mattis 1961", "Pfeuty 1970"],
    notes="Site-resolved ⟨σᶻ_i⟩ from Majorana thermal covariance on OBC chain.",
)
