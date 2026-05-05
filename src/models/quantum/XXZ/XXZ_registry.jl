# models/quantum/XXZ/XXZ_registry.jl — declarative implementation map.
#
# One `@register` line per natively-implemented (model, quantity, bc)
# triple.  See `src/core/registry.jl` for the metadata schema and
# `src/models/quantum/TFIM/TFIM_registry.jl` for the canonical example.
#
# `:dense_ed` rows are tagged `:high` reliability — finite-N ED is exact
# for any N ≤ `_MAX_ED_SITES`, so the only failure mode is a cap miss
# (caught by an explicit `ArgumentError`).

# ── Energy (granularity-aware) ─────────────────────────────────────────
@register(
    XXZ1D,
    Energy{:total},
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_thermal.jl",
    references=["Yang Yang 1966", "Takahashi 1999"],
    notes="Total ⟨H⟩(β) by dense ED of the 2^N × 2^N XXZ Hamiltonian.",
)
@register(
    XXZ1D,
    Energy{:per_site},
    Infinite,
    method=:bethe_ansatz,
    reliability=:high,
    tested_in="test/models/test_XXZ1D.jl",
    references=["Hulthén 1938", "Yang Yang 1966"],
    notes="Closed form at Δ ∈ {-1, 0, 1}; general-Δ Bethe-ansatz integral deferred.",
)

# ── Spectrum / criticality ────────────────────────────────────────────
@register(
    XXZ1D,
    CentralCharge,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/test_XXZ1D.jl",
    references=["Giamarchi 2004"],
    notes="c = 1 in the critical regime -1 < Δ < 1.",
)
@register(
    XXZ1D,
    LuttingerParameter,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/test_XXZ1D.jl",
    references=["Giamarchi 2004"],
    notes="K = π / (2(π − arccos Δ)) for -1 < Δ ≤ 1.",
)
@register(
    XXZ1D,
    LuttingerVelocity,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/test_XXZ1D.jl",
    references=["Giamarchi 2004", "des Cloizeaux Pearson 1962"],
    notes="u = (πJ/2) sin γ / γ, γ = arccos Δ; -1 < Δ ≤ 1.",
)
@register(
    XXZ1D,
    MassGap,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl",
    notes="E₁ - E₀ from full-spectrum dense ED.",
)
@register(
    XXZ1D,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl",
    references=["Giamarchi 2004"],
    notes="0 in the critical regime -1 < Δ ≤ 1; gapped regime returns NaN with a warning.",
)

# ── Finite-T thermodynamic scalars (per-site at OBC) ──────────────────
@register(
    XXZ1D,
    FreeEnergy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_thermal.jl"
)
@register(
    XXZ1D,
    ThermalEntropy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_thermal.jl"
)
@register(
    XXZ1D,
    SpecificHeat,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_thermal.jl"
)

# ── Magnetisations (Pauli convention) ─────────────────────────────────
@register(
    XXZ1D,
    MagnetizationX,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)
@register(
    XXZ1D,
    MagnetizationY,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)
@register(
    XXZ1D,
    MagnetizationZ,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)

# ── Site-resolved local observables ───────────────────────────────────
@register(
    XXZ1D,
    MagnetizationXLocal,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)
@register(
    XXZ1D,
    MagnetizationYLocal,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)
@register(
    XXZ1D,
    MagnetizationZLocal,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)
@register(
    XXZ1D,
    EnergyLocal,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl",
    notes="Bonds split symmetrically: Σᵢ ε_i = ⟨H⟩.",
)

# ── Susceptibilities (β · Var(M_α) / N) ───────────────────────────────
@register(
    XXZ1D,
    SusceptibilityXX,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)
@register(
    XXZ1D,
    SusceptibilityYY,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)
@register(
    XXZ1D,
    SusceptibilityZZ,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)

# ── Two-point correlators (static + connected) ────────────────────────
for CorrTy in (XXCorrelation, YYCorrelation, ZZCorrelation), mode in (:static, :connected)
    register!(
        XXZ1D,
        CorrTy{mode},
        OBC;
        method=:dense_ed,
        reliability=:high,
        tested_in="test/models/test_XXZ1D_observables.jl",
        notes="(i,j) ⟨σᵅᵢ σᵅⱼ⟩_β; mode=:connected subtracts the disconnected piece.",
    )
end

# ── Entanglement (β = Inf default → ground-state pure-state entropy) ──
@register(
    XXZ1D,
    VonNeumannEntropy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl",
    notes="Pass subsystem length ℓ; β=Inf gives ground-state EE.",
)
@register(
    XXZ1D,
    RenyiEntropy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl",
    notes="S_α = log Tr ρ_A^α / (1 - α); pass subsystem ℓ and order α.",
)
