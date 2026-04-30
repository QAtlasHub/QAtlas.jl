# models/quantum/Heisenberg/HeisenbergS1_registry.jl —
# declarative implementation map for S1Heisenberg1D (Haldane chain).
#
# Every native fetch method declared in `HeisenbergS1.jl` and
# `HeisenbergS1_observables.jl` gets one `@register` line here.  Method
# tag conventions match `TFIM_registry.jl` (`:dense_ed`,
# `:literature_value`, `:high`/`:medium` reliability).

# ── Energy / FreeEnergy / Entropy / SpecificHeat (OBC, dense ED) ──────
@register(
    S1Heisenberg1D,
    Energy{:total},
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_thermal.jl",
    notes="Total ⟨H⟩(β) by dense ED on the 3^N Hilbert space; N ≤ 8.",
)
@register(
    S1Heisenberg1D,
    FreeEnergy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_thermal.jl",
)
@register(
    S1Heisenberg1D,
    ThermalEntropy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_thermal.jl",
)
@register(
    S1Heisenberg1D,
    SpecificHeat,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_thermal.jl",
)

# ── Per-site magnetisations (OBC) ─────────────────────────────────────
@register(
    S1Heisenberg1D,
    MagnetizationX,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
)
@register(
    S1Heisenberg1D,
    MagnetizationY,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
)
@register(
    S1Heisenberg1D,
    MagnetizationZ,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
)

# ── Susceptibilities (OBC) ────────────────────────────────────────────
@register(
    S1Heisenberg1D,
    SusceptibilityXX,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
)
@register(
    S1Heisenberg1D,
    SusceptibilityYY,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
)
@register(
    S1Heisenberg1D,
    SusceptibilityZZ,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
)

# ── Two-point correlators (static / connected, OBC) ───────────────────
@register(
    S1Heisenberg1D,
    XXCorrelation{:static},
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
)
@register(
    S1Heisenberg1D,
    YYCorrelation{:static},
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
)
@register(
    S1Heisenberg1D,
    ZZCorrelation{:static},
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
)
@register(
    S1Heisenberg1D,
    XXCorrelation{:connected},
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
)
@register(
    S1Heisenberg1D,
    YYCorrelation{:connected},
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
)
@register(
    S1Heisenberg1D,
    ZZCorrelation{:connected},
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
)

# ── Local one-site / one-bond observables (OBC) ───────────────────────
@register(
    S1Heisenberg1D,
    MagnetizationXLocal,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
)
@register(
    S1Heisenberg1D,
    MagnetizationZLocal,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
)
@register(
    S1Heisenberg1D,
    EnergyLocal,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
)

# ── Spectrum / criticality ────────────────────────────────────────────
@register(
    S1Heisenberg1D,
    MassGap,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
    notes="E₁ - E₀ from dense ED on 3^N space; finite-N + edge-state corrections to Δ_∞.",
)
@register(
    S1Heisenberg1D,
    MassGap,
    Infinite,
    method=:literature_value,
    reliability=:medium,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
    references=["White-Huse 1993"],
    notes="Haldane gap Δ ≈ 0.41048 J (DMRG; no closed form).",
)
@register(
    S1Heisenberg1D,
    Energy{:per_site},
    Infinite,
    method=:literature_value,
    reliability=:medium,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
    references=["White-Huse 1993"],
    notes="GS energy density e₀ ≈ -1.401484 J (DMRG; no closed form).",
)

# ── Entanglement (T = 0; β kwarg defaults to Inf) ─────────────────────
@register(
    S1Heisenberg1D,
    VonNeumannEntropy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
    notes="Partial trace of dense thermal ρ; subsystem length ℓ ∈ [1, N-1].",
)
@register(
    S1Heisenberg1D,
    RenyiEntropy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_S1Heisenberg1D_observables.jl",
    notes="S_α = log Tr ρ_A^α / (1-α); same partial-trace path as VonNeumann.",
)
