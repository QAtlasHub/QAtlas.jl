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
    status=:approx,
    reliability=:medium,
    valid_domain="η > 0 Lorentzian broadening regulator; the physical η→0⁺ limit is cut off by the 1D band-edge van Hove singularities, so the rate is η-dependent and is reported at finite η.",
    tested_in="test/models/quantum/misc/test_tight_binding1d.jl",
    references=["Korringa1950"],
    notes="η-regularized free-fermion (Korringa-type) golden-rule rate 1/T₁(β,η) = (1/π³)∫₀^π∫₀^π f(ε₁)(1-f(ε₂)) η/((ε₁-ε₂)²+η²) dk₁dk₂ — the q-summed dynamic structure factor S(q,ω→0) of the 1D tight-binding band, Lorentzian-broadened by η. Evaluated by nested QuadGK (rtol=1e-6); verified against an independent discrete N-mode k-sum and the β→0 ¼-factorisation.",
)
