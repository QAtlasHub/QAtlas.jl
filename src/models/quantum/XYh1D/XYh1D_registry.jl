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

# ── Transverse Magnetization & Susceptibility ─────────────────────────
@register(
    XYh1D,
    MagnetizationZ,
    Infinite,
    method=:quadgk,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_xyh1d.jl",
    references=["Lieb-Schultz-Mattis 1961", "Pfeuty 1970"],
    notes="transverse magnetization per site (along the field direction z) via QuadGK.",
)
@register(
    XYh1D,
    MagnetizationZ,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_xyh1d.jl",
    references=["Lieb-Schultz-Mattis 1961", "Pfeuty 1970"],
    notes="transverse magnetization per site from the thermal Majorana covariance matrix.",
)
@register(
    XYh1D,
    SusceptibilityZZ,
    Infinite,
    method=:quadgk,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_xyh1d.jl",
    references=["Lieb-Schultz-Mattis 1961", "Pfeuty 1970"],
    notes="transverse susceptibility per site via QuadGK.",
)
@register(
    XYh1D,
    SusceptibilityZZ,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_xyh1d.jl",
    references=["Lieb-Schultz-Mattis 1961", "Pfeuty 1970"],
    notes="transverse susceptibility per site from exact Wick contraction over covariance.",
)

# ── Site-local Equilibrium Observables ────────────────────────────────
@register(
    XYh1D,
    MagnetizationZLocal,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_xyh1d.jl",
    notes="Site-local expectation value vector of σᶻ.",
)
@register(
    XYh1D,
    MagnetizationXLocal{:equilibrium},
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_xyh1d.jl",
    notes="Identically zero by Z₂ symmetry.",
)
@register(
    XYh1D,
    MagnetizationYLocal,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_xyh1d.jl",
    notes="Identically zero by Z₂ symmetry.",
)
@register(
    XYh1D,
    EnergyLocal,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_xyh1d.jl",
    notes="Symmetric split of bond energies ε_i such that sum(ε) == ⟨H⟩.",
)

# ── PBC ──────────────────────────────────────────────────────────────
@register(
    XYh1D,
    MassGap,
    PBC,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_xyh1d_pbc.jl",
    references=["Lieb-Schultz-Mattis 1961", "Pfeuty 1970"],
    notes="Minimum quasiparticle energy over both PBC sectors (AP and P).",
)
