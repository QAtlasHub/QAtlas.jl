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
    references=["AshcroftMermin1976"],
    notes="E/N = -(2t/π) sin(k_F) - (μ/π) k_F at partial filling; E/N = 0 (empty) / -μ (full).",
)

@register(
    TightBinding1D,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_tight_binding1d.jl",
    references=["AshcroftMermin1976"],
    notes="Δ = max(0, |μ| - 2t); gapless for |μ| ≤ 2t (metallic).",
)

@register(
    TightBinding1D,
    FermiVelocity,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_tight_binding1d.jl",
    references=["AshcroftMermin1976"],
    notes="v_F = 2t sin(k_F) = 2t√(1 - μ²/(4t²)); returns 0 in gapped phase |μ| ≥ 2t (no Fermi surface; convention).",
)
@register(
    TightBinding1D,
    FreeEnergy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_tight_binding1d.jl",
    references=["Mahan2000", "Coleman2015"],
    notes="ω(β;t,μ) = -(πβ)⁻¹ ∫₀^π log(1+e^{-βε}) dk, ε(k)=-2t cos k - μ; QuadGK rtol=1e-10.",
)

@register(
    TightBinding1D,
    ThermalEntropy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_tight_binding1d.jl",
    references=["Mahan2000"],
    notes="s(β;t,μ) = β(u-ω); high-T limit log 2 per site, Sommerfeld linear in T at low T.",
)

@register(
    TightBinding1D,
    SpecificHeat,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_tight_binding1d.jl",
    references=["Mahan2000"],
    notes="c_μ(β;t,μ) = (β²/π) ∫₀^π ε² n_F(1-n_F) dk; QuadGK rtol=1e-10.",
)

@register(
    TightBinding1D,
    NMRSpinRelaxationRate,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_tight_binding1d.jl",
    notes="1/T_1(β, η) = 1/π³ ∫₀^π dk₁ ∫₀^π dk₂ f(ε(k₁)) (1-f(ε(k₂))) η / ((ε(k₁)-ε(k₂))² + η²); regularized 1D Korringa rate; QuadGK nested rtol=1e-6.",
)
