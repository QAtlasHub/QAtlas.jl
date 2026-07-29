# models/classical/TodaLattice/TodaLattice_registry.jl
#
# Declarative implementation map for the 1-D Toda lattice
# (Toda 1967, classical integrable).  Schema in src/core/registry.jl.

@register(
    TodaLattice,
    MassGap,
    Infinite,
    method=:linear_phonon,
    cost=:closed_form,
    status=:approx,
    valid_domain="harmonic (linearised) branch only -- small oscillations about the " *
                 "equilibrium chain",
    error_order="the LINEARISATION is the approximation, not a small parameter: the " *
                "Toda lattice is integrable and its spectrum also carries a SOLITON " *
                "sector that this branch does not see, so MassGap = 0 is a statement " *
                "about the acoustic phonons alone (quantum Toda tracked as Phase 2).",
    reliability=:high,
    tested_in="test/standalone/test_toda_lattice.jl",
    references=["Toda1967", "Flaschka1974"],
    notes="Linearised acoustic phonon ω(k) = 2√(ab)|sin(k/2)| ⇒ MassGap = 0.  Soliton / quantum-Toda spectrum tracked as Phase 2.",
)
