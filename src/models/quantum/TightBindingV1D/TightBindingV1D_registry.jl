# models/quantum/TightBindingV1D/TightBindingV1D_registry.jl
#
# Declarative implementation map for the 1D spinless-fermion t-V chain.
# Schema documented in src/core/registry.jl.

@register(
    TightBindingV1D,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_tight_binding_v1d.jl",
    references=["YangYang1966", "Ashcroft-Mermin 1976"],
    notes="V=0 free-fermion gap = max(0, |μ|-2t); V≠0 (JW XXZ) deferred to Phase 2.",
)

@register(
    TightBindingV1D,
    FermiVelocity,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_tight_binding_v1d.jl",
    references=["Ashcroft-Mermin 1976"],
    notes="V=0 free-fermion v_F = 2t sin(k_F); |μ|≥2t (no Fermi surface) DomainError.",
)

@register(
    TightBindingV1D,
    Energy{:per_site},
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_tight_binding_v1d.jl",
    references=["Mahan2000", "Ashcroft-Mermin 1976"],
    notes="V=0 free-fermion e₀ = -(2t/π)sin(k_F) - (μ/π)k_F, k_F = arccos(-μ/2t); band-edge piecewise.",
)
@register(
    TightBindingV1D,
    FreeEnergy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_tight_binding_v1d.jl",
    references=["Mahan2000", "Coleman2015"],
    notes="V=0 free-fermion ω(β;t,μ) = -(πβ)⁻¹ ∫₀^π log(1+e^{-βε}) dk; V≠0 (JW XXZ, Yang-Yang 1966) deferred to Phase 2.",
)

@register(
    TightBindingV1D,
    ThermalEntropy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_tight_binding_v1d.jl",
    references=["Mahan2000"],
    notes="V=0 s(β;t,μ) = β(u-ω); high-T limit log 2 per site; V≠0 deferred to Phase 2.",
)

@register(
    TightBindingV1D,
    SpecificHeat,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_tight_binding_v1d.jl",
    references=["Mahan2000"],
    notes="V=0 c_μ(β;t,μ) = (β²/π) ∫₀^π ε² n_F(1-n_F) dk; V≠0 deferred to Phase 2.",
)
