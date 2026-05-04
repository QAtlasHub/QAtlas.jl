# models/quantum/TFIM/S1TFIM_registry.jl —
# declarative implementation map for S1TFIM (spin-1 transverse-field Ising chain).
#
# Method tag conventions match `TFIM_registry.jl` and
# `HeisenbergS1_registry.jl` (`:dense_ed`, `:high` reliability for the
# finite-N exact path).

@register(
    S1TFIM,
    Energy{:total},
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/standalone/test_s1_tfim.jl",
    notes="Total ⟨H⟩(β) by dense ED on the 3^N Hilbert space; N ≤ 8.",
)
@register(
    S1TFIM,
    FreeEnergy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/standalone/test_s1_tfim.jl",
)
@register(
    S1TFIM,
    ThermalEntropy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/standalone/test_s1_tfim.jl",
)
@register(
    S1TFIM,
    SpecificHeat,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/standalone/test_s1_tfim.jl",
)

# Per-site magnetisations
@register(
    S1TFIM,
    MagnetizationX,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/standalone/test_s1_tfim.jl",
)
@register(
    S1TFIM,
    MagnetizationY,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/standalone/test_s1_tfim.jl",
)
@register(
    S1TFIM,
    MagnetizationZ,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/standalone/test_s1_tfim.jl",
)

# Susceptibilities
@register(
    S1TFIM,
    SusceptibilityXX,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/standalone/test_s1_tfim.jl",
)
@register(
    S1TFIM,
    SusceptibilityYY,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/standalone/test_s1_tfim.jl",
)
@register(
    S1TFIM,
    SusceptibilityZZ,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/standalone/test_s1_tfim.jl",
)

# ZZ correlator
@register(
    S1TFIM,
    ZZCorrelation{:static},
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/standalone/test_s1_tfim.jl",
)
@register(
    S1TFIM,
    ZZCorrelation{:connected},
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/standalone/test_s1_tfim.jl",
)

# Spectrum
@register(
    S1TFIM,
    MassGap,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/standalone/test_s1_tfim.jl",
    notes="E₁ - E₀ from dense ED on the 3^N space.",
)
@register(
    S1TFIM,
    ExactSpectrum,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/standalone/test_s1_tfim.jl",
    notes="Sorted eigenvalues of the dense 3^N × 3^N Hamiltonian.",
)