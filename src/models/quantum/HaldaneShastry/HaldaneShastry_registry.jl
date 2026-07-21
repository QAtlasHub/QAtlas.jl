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
    references=["Haldane1988", "Shastry1988"],
    notes="E_0/N = -π² J / 24, finite-N exact (Gutzwiller-projected free-fermion eigenstate of the chord-distance Hamiltonian).",
)

# ── Finite-T at Infinite() via c=1 CFT low-T (#524 stopgap) ────────────
# Mirrors the Heisenberg1D Path B (#521 / PR #526): same v_s = π J / 2,
# different e_0. Valid β > 5/J. The full Haldane semion gas will replace
# these once #524 is closed.

@register(
    HaldaneShastry,
    FreeEnergy,
    Infinite,
    method=:cft_low_T,
    reliability=:medium,
    tested_in="test/models/quantum/HaldaneShastry/test_haldane_shastry_thermal_cft.jl",
    references=["Haldane1988", "Affleck1986"],
    notes="f = -π² J / 24 - π T² / (6 v_s), v_s = π J / 2. Valid β > 5/J.",
)

@register(
    HaldaneShastry,
    ThermalEntropy,
    Infinite,
    method=:cft_low_T,
    reliability=:medium,
    tested_in="test/models/quantum/HaldaneShastry/test_haldane_shastry_thermal_cft.jl",
    references=["Haldane1988", "Affleck1986"],
    notes="s = π T / (3 v_s) = 2T / (3J). Valid β > 5/J.",
)

@register(
    HaldaneShastry,
    SpecificHeat,
    Infinite,
    method=:cft_low_T,
    reliability=:medium,
    tested_in="test/models/quantum/HaldaneShastry/test_haldane_shastry_thermal_cft.jl",
    references=["Haldane1988", "Affleck1986"],
    notes="c_v = π T / (3 v_s) = 2T / (3J). Equals s at LO CFT. Valid β > 5/J.",
)

# ── CC entanglement at Infinite via Universality(:Heisenberg) (#580 Phase 2)
@register(
    HaldaneShastry,
    VonNeumannEntropy,
    Infinite,
    method=:delegation,
    reliability=:high,
    tested_in="test/models/quantum/HaldaneShastry/test_haldaneshastry_cft_entanglement.jl",
    references=["CalabreseCardy2004"],
    notes="Delegates to Universality(:Heisenberg) c=1 Calabrese-Cardy form. HS is gapless free-spinon SU(2)_1 WZW, same c=1 class as Heisenberg1D.",
)

@register(
    HaldaneShastry,
    RenyiEntropy,
    Infinite,
    method=:delegation,
    reliability=:high,
    tested_in="test/models/quantum/HaldaneShastry/test_haldaneshastry_cft_entanglement.jl",
    references=["CalabreseCardy2004"],
    notes="Delegates to Universality(:Heisenberg) with c -> c*(1+1/alpha)/2 substitution. Reduces to VN at alpha=1.",
)
