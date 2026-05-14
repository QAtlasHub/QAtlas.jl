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
    references=["Yang-Yang 1966", "Ashcroft-Mermin 1976"],
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
    references=["Mahan 2000", "Ashcroft-Mermin 1976"],
    notes="V=0 free-fermion e₀ = -(2t/π)sin(k_F) - (μ/π)k_F, k_F = arccos(-μ/2t); band-edge piecewise.",
)
