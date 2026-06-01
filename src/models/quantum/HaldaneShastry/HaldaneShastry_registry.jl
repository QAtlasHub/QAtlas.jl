# models/quantum/HaldaneShastry/HaldaneShastry_registry.jl — declarative
# implementation map for the spin-1/2 Haldane-Shastry chain.
#
# Phase 1 (this file) covers thermodynamic-limit closed forms at
# Infinite(): ground-state energy density. Finite-T thermodynamics
# via free-spinon semion gas is tracked in #524.

@register(
    HaldaneShastry,
    GroundStateEnergyDensity,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/HaldaneShastry/test_haldane_shastry.jl",
    references=["Haldane PRL 60, 635 (1988)", "Shastry PRL 60, 639 (1988)"],
    notes="E_0/N = -π² J / 24, finite-N exact (Gutzwiller-projected free-fermion eigenstate of the chord-distance Hamiltonian).",
)
