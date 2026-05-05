# models/quantum/Heisenberg/Heisenberg_registry.jl — declarative implementation map.
#
# Heisenberg1D (spin-1/2) reuses the XXZ1D(Δ=1) finite-N OBC implementations
# as thin delegators (see `Heisenberg.jl`).  Every row below is therefore a
# reflection of the corresponding XXZ1D row, with `:dense_ed` reliability
# inherited.  The thermodynamic-limit ground-state energy density at the
# isotropic point is the original Hulthén result.

# ── Closed-form ground state in the thermodynamic limit ───────────────
@register(
    Heisenberg1D,
    GroundStateEnergyDensity,
    Infinite,
    method=:bethe_ansatz,
    reliability=:high,
    tested_in="test/standalone/test_bethe_ansatz.jl",
    references=["Hulthén 1938", "Bethe 1931"],
    notes="e₀ = J(1/4 - ln 2) at the isotropic AF point.",
)

# ── Energy (delegates to XXZ1D(Δ=1)) ──────────────────────────────────
@register(
    Heisenberg1D,
    Energy{:total},
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl",
    notes="Delegates to XXZ1D(Δ=1.0); J passed via kwargs.",
)

# ── Spectrum / criticality ────────────────────────────────────────────
@register(
    Heisenberg1D,
    MassGap,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl",
    notes="Delegates to XXZ1D(Δ=1.0).",
)
@register(
    Heisenberg1D,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl",
    notes="Gapless (0.0) at the isotropic critical point.",
)

# ── Finite-T thermodynamic scalars (per-site at OBC) ──────────────────
@register(
    Heisenberg1D,
    FreeEnergy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl"
)
@register(
    Heisenberg1D,
    ThermalEntropy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl"
)
@register(
    Heisenberg1D,
    SpecificHeat,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl"
)

# ── Magnetisations ────────────────────────────────────────────────────
@register(
    Heisenberg1D,
    MagnetizationX,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl"
)
@register(
    Heisenberg1D,
    MagnetizationY,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl"
)
@register(
    Heisenberg1D,
    MagnetizationZ,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl"
)

# ── Local site-resolved observables ───────────────────────────────────
@register(
    Heisenberg1D,
    MagnetizationXLocal,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl"
)
@register(
    Heisenberg1D,
    MagnetizationYLocal,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl"
)
@register(
    Heisenberg1D,
    MagnetizationZLocal,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl"
)
@register(
    Heisenberg1D,
    EnergyLocal,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl"
)

# ── Susceptibilities ──────────────────────────────────────────────────
@register(
    Heisenberg1D,
    SusceptibilityXX,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl"
)
@register(
    Heisenberg1D,
    SusceptibilityYY,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl"
)
@register(
    Heisenberg1D,
    SusceptibilityZZ,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl"
)

# ── Two-point correlators (static + connected) ────────────────────────
for CorrTy in (XXCorrelation, YYCorrelation, ZZCorrelation), mode in (:static, :connected)
    register!(
        Heisenberg1D,
        CorrTy{mode},
        OBC;
        method=:dense_ed,
        reliability=:high,
        tested_in="test/models/test_Heisenberg1D_thermal.jl",
    )
end

# ── Entanglement ──────────────────────────────────────────────────────
@register(
    Heisenberg1D,
    VonNeumannEntropy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl"
)
@register(
    Heisenberg1D,
    RenyiEntropy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_Heisenberg1D_thermal.jl"
)
