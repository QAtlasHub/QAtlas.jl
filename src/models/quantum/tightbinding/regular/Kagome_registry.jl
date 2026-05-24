# models/quantum/tightbinding/regular/Kagome_registry.jl
#
# Declarative implementation map for the regular kagome tight-binding
# tile.  Schema documented in src/core/registry.jl.

@register(
    Kagome,
    TightBindingChecksum,
    Infinite,
    method=:bloch_diagonalization,
    reliability=:high,
    tested_in="test/models/quantum/tightbinding/test_kagome_tb_scalar_invariants.jl",
    notes="Σ λᵢ² = tr(H²); forwarded through TightBindingSpectrum so it agrees with the spectrum closed-form by construction; verify cards pin it against the chiral-symmetry identity 2 t² · z · Lx · Ly and against real-space ED.",
)

@register(
    Kagome,
    TightBindingMaxEnergy,
    Infinite,
    method=:bloch_diagonalization,
    reliability=:high,
    tested_in="test/models/quantum/tightbinding/test_kagome_tb_scalar_invariants.jl",
    notes="max(λᵢ); forwarded through TightBindingSpectrum.  Sister scalar to TightBindingChecksum: together they discriminate single-eigenvalue perturbations.",
)
