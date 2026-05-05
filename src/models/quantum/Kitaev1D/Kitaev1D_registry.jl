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
    references=["Kitaev 2001"],
    notes="Per-site ε₀ by Gauss-Kronrod over the PBC dispersion E(k).",
)

# ── Spectrum / criticality ────────────────────────────────────────────
@register(
    Kitaev1D,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_kitaev1d.jl",
    references=["Kitaev 2001", "Alicea 2012"],
    notes="Closed-form min over k of √((2t cos k + μ)² + 4Δ² sin² k).",
)
@register(
    Kitaev1D,
    MassGap,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/standalone/test_kitaev1d.jl",
    references=["Kitaev 2001"],
    notes="Smallest non-negative BdG eigenvalue (Majorana edge mode in topological phase).",
)
@register(
    Kitaev1D,
    EdgeModeEnergy,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/standalone/test_kitaev1d.jl",
    references=["Kitaev 2001", "Alicea 2012"],
    notes="Same value as MassGap@OBC; named for the Majorana boundary-mode interpretation.",
)
@register(
    Kitaev1D,
    CorrelationLength,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_kitaev1d.jl",
    references=["Kitaev 2001"],
    notes="ξ = 1/Δ_gap; Inf on the critical line |μ| = 2|t|.",
)
@register(
    Kitaev1D,
    TopologicalInvariant,
    Infinite,
    method=:pfaffian,
    reliability=:high,
    tested_in="test/standalone/test_kitaev1d.jl",
    references=["Kitaev 2001", "Asboth-Oroszlany-Palyi 2016"],
    notes="ν = sgn[Pf A(k=0)·Pf A(k=π)] = sgn(μ² - 4t²); ν=-1 topological, +1 trivial.",
)
