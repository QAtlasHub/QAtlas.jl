# models/quantum/Kitaev1D/Kitaev1D_registry.jl — declarative implementation
# map for the 1D Majorana wire (Kitaev 2001).  See `src/core/registry.jl`
# for the metadata schema.

# ── Energy (granularity-aware) ─────────────────────────────────────────
@register(
    Kitaev1D,
    Energy{:per_site},
    Infinite,
    method=:bdg,
    reliability=:high,
    tested_in="test/standalone/test_kitaev1d.jl",
    references=["Kitaev2001"],
    notes="Per-site ε₀ by Gauss-Kronrod over the PBC dispersion E(k).",
)

# ── Free-fermion (BdG) finite-T thermodynamics — Kitaev1D_thermal.jl ───
@register(
    Kitaev1D,
    FreeEnergy,
    Infinite,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_kitaev1d_thermo.jl",
    references=["Kitaev2001"],
    notes="BdG f(β) = -(1/2π)∫[E/2 + β⁻¹ln(1+e^{-βE})]dk over the PBC dispersion; β→∞ → ε₀.",
)
@register(
    Kitaev1D,
    SpecificHeat,
    Infinite,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_kitaev1d_thermo.jl",
    references=["Kitaev2001"],
    notes="c_v = β²Var(H)/L = (β²/8π)∫E²sech²(βE/2)dk (free-fermion energy FDT).",
)
@register(
    Kitaev1D,
    ThermalEntropy,
    Infinite,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_kitaev1d_thermo.jl",
    references=["Kitaev2001"],
    notes="s(β) = β(ε-f); bounded 0 (T→0) … ln2 (T→∞, two states per spinless mode).",
)

# ── Spectrum / criticality ────────────────────────────────────────────
@register(
    Kitaev1D,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_kitaev1d.jl",
    references=["Kitaev2001", "Alicea2012"],
    notes="Closed-form min over k of √((2t cos k + μ)² + 4Δ² sin² k).",
)
@register(
    Kitaev1D,
    MassGap,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/standalone/test_kitaev1d.jl",
    references=["Kitaev2001"],
    notes="Smallest non-negative BdG eigenvalue (Majorana edge mode in topological phase).",
)
@register(
    Kitaev1D,
    EdgeModeEnergy,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/standalone/test_kitaev1d.jl",
    references=["Kitaev2001", "Alicea2012"],
    notes="Same value as MassGap@OBC; named for the Majorana boundary-mode interpretation.",
)
@register(
    Kitaev1D,
    CorrelationLength,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_kitaev1d.jl",
    references=["Kitaev2001"],
    notes="ξ = 1/Δ_gap; Inf on the critical line |μ| = 2|t|.",
)
@register(
    Kitaev1D,
    TopologicalInvariant,
    Infinite,
    method=:pfaffian,
    reliability=:high,
    tested_in="test/standalone/test_kitaev1d.jl",
    references=["Kitaev2001", "AsbothOroszlanyPalyi2016"],
    notes="ν = sgn[Pf A(k=0)·Pf A(k=π)] = sgn(μ² - 4t²); ν=-1 topological, +1 trivial.",
)
