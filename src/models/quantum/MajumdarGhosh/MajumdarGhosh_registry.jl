# models/quantum/MajumdarGhosh/MajumdarGhosh_registry.jl — declarative
# implementation map for the Majumdar–Ghosh chain (S = 1/2 J₁–J₂ chain at
# J₂/J₁ = 1/2).  Mirrors the layout of `Heisenberg_registry.jl` and
# `KitaevHoneycomb_registry.jl`.  See `src/core/registry.jl` for the
# metadata schema.

# ── Closed-form ground state energy density ───────────────────────────
@register(
    MajumdarGhosh,
    GroundStateEnergyDensity,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_majumdar_ghosh.jl",
    references=["Majumdar-Ghosh 1969"],
    notes="E₀/N = -3J/8 from the exact dimer-product ground state.",
)
@register(
    MajumdarGhosh,
    GroundStateEnergyDensity,
    PBC,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_majumdar_ghosh.jl",
    references=["Majumdar-Ghosh 1969"],
    notes="Size-independent: dimer-product state is an exact eigenstate for any even N.",
)

# ── Spectrum / criticality ────────────────────────────────────────────
@register(
    MajumdarGhosh,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_majumdar_ghosh.jl",
    references=["Shastry-Sutherland 1981", "White-Affleck 1996"],
    notes="method=:lower_bound → J/4 (Shastry-Sutherland); :numerical → 0.234 J (White-Affleck DMRG).",
)

@register(
    MajumdarGhosh,
    SpinGap,
    Infinite,
    method=:dmrg_reference,
    reliability=:medium,
    tested_in="test/models/quantum/Heisenberg/test_majumdar_ghosh.jl",
    references=["White-Affleck 1996", "Eggert 1996"],
    notes="Δ ≈ 0.234 J DMRG; Shastry-Sutherland 1981 bound Δ ≥ J/4 is trimer-sector only.",
)
