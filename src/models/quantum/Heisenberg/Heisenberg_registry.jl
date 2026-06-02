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
    MagnetizationXLocal{:equilibrium},
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
    VonNeumannEntropy{:equilibrium},
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

# ── Spinon kinematics (issue #154 phase 1, Infinite) ──────────────────
# Single-spinon dispersion + des Cloizeaux–Pearson 2-spinon continuum
# edges are exposed as top-level helpers
# (heisenberg_spinon_dispersion, heisenberg_two_spinon_lower_edge,
#  heisenberg_two_spinon_upper_edge); only the dynamic structure factor
# is registered as a Quantity here.
@register(
    Heisenberg1D,
    ZZStructureFactor,
    Infinite,
    method=:muller_ansatz,
    reliability=:medium,
    tested_in="test/standalone/test_heisenberg_spinon.jl",
    references=["des Cloizeaux–Pearson 1962", "Müller-Thomas-Beck-Bonner 1981"],
    notes="Phase 1 closed-form Müller ansatz for S^{zz}(q,ω); exact Caux–Hagemans 2006 result reserved for Phase 2.",
)

# ── Luttinger-liquid parameter at SU(2)-symmetric point (Phase 2) ─────
@register(
    Heisenberg1D,
    LuttingerParameter,
    Infinite,
    method=:delegation,
    reliability=:high,
    tested_in="test/models/quantum/Heisenberg/test_heisenberg1d_luttinger.jl",
    references=["Luther-Peschel 1975", "Affleck 1989", "Haldane 1980"],
    notes="K = 1/2 (SU(2)-symmetric Heisenberg AFM); delegate to XXZ1D(Δ=1).",
)

# ── Finite-T at Infinite() via c=1 CFT low-T expansion (#521 Path B) ──
# Stopgap: valid for β > 5/J only. The Klümper Δ → 1 limit (Path A) will
# replace these once implemented.

@register(
    Heisenberg1D,
    FreeEnergy,
    Infinite,
    method=:cft_low_T,
    reliability=:medium,
    tested_in="test/models/quantum/Heisenberg/test_heisenberg1d_thermal_cft.jl",
    references=[
        "Affleck PRL 56, 746 (1986)",
        "Blöte-Cardy-Nightingale PRL 56, 742 (1986)",
        "Eggert-Affleck-Takahashi PRL 73, 332 (1994)",
    ],
    notes="f = e₀ - π T² / (6 v_s), v_s = π J / 2. Valid β > 5/J; NaN+warn otherwise.",
)

@register(
    Heisenberg1D,
    ThermalEntropy,
    Infinite,
    method=:cft_low_T,
    reliability=:medium,
    tested_in="test/models/quantum/Heisenberg/test_heisenberg1d_thermal_cft.jl",
    references=["Affleck PRL 56, 746 (1986)"],
    notes="s = π T / (3 v_s) = 2T / (3J). Valid β > 5/J.",
)

@register(
    Heisenberg1D,
    SpecificHeat,
    Infinite,
    method=:cft_low_T,
    reliability=:medium,
    tested_in="test/models/quantum/Heisenberg/test_heisenberg1d_thermal_cft.jl",
    references=["Affleck PRL 56, 746 (1986)"],
    notes="c_v = π T / (3 v_s) = 2T / (3J). Equals s(T) at LO CFT. Valid β > 5/J.",
)
