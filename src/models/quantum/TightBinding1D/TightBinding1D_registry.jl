# models/quantum/TightBinding1D/TightBinding1D_registry.jl
#
# Declarative implementation map for the 1D non-interacting spinless-fermion
# chain TightBinding1D (Phase 1 — three closed-form quantities at Infinite).
# Schema documented in src/core/registry.jl.

@register(
    TightBinding1D,
    Energy{:per_site},
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_tight_binding1d.jl",
    references=["Ashcroft-Mermin 1976"],
    notes="E/N = -(2t/π) sin(k_F) - (μ/π) k_F at partial filling; E/N = 0 (empty) / -μ (full).",
)

@register(
    TightBinding1D,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_tight_binding1d.jl",
    references=["Ashcroft-Mermin 1976"],
    notes="Δ = max(0, |μ| - 2t); gapless for |μ| ≤ 2t (metallic).",
)

@register(
    TightBinding1D,
    FermiVelocity,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_tight_binding1d.jl",
    references=["Ashcroft-Mermin 1976"],
    notes="v_F = 2t sin(k_F) = 2t√(1 - μ²/(4t²)); returns 0 in gapped phase |μ| ≥ 2t (no Fermi surface; convention).",
)
