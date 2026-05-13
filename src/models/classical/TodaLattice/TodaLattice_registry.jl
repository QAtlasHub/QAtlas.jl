# models/classical/TodaLattice/TodaLattice_registry.jl
#
# Declarative implementation map for the 1-D Toda lattice
# (Toda 1967, classical integrable).  Schema in src/core/registry.jl.

@register(
    TodaLattice,
    MassGap,
    Infinite,
    method=:linear_phonon,
    reliability=:high,
    tested_in="test/standalone/test_toda_lattice.jl",
    references=["Toda 1967", "Flaschka 1974"],
    notes="Linearised acoustic phonon ω(k) = 2√(ab)|sin(k/2)| ⇒ MassGap = 0.  Soliton / quantum-Toda spectrum tracked as Phase 2.",
)
