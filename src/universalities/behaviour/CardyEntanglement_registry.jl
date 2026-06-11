# universalities/CardyEntanglement_registry.jl — the `predicts` edges of the
# knowledge graph for the cross-class 1+1D CFT quantities.
#
# `Universality{class} --predicts--> quantity` rows (status=:universal), one per
# (CFT class × quantity).  These make the universality layer first-class in the
# registry (queryable via `implementation_status` / `predicts`), symmetric with
# models/ and bounds/, and feed the graph-derived coherence checks.  The fetch
# bodies live in CardyEntanglement.jl; here we register what each class predicts.

# 1+1D CFT classes carrying a central charge (the classes the cross-class
# Calabrese–Cardy machinery applies to).
const _CFT_CLASSES = (:Ising, :XY, :Heisenberg, :Potts3, :Potts4)

# (quantity, references, note) — all at the Infinite (thermodynamic) BC.
const _CFT_PREDICTIONS = (
    (
        MutualInformation,
        ["CalabreseCardy2009"],
        "Adjacent-interval mutual information (Calabrese-Cardy 2009).",
    ),
    (
        EntanglementGrowthSlope,
        ["CalabreseCardy2005"],
        "Linear entanglement-growth slope pi c v/(3 beta_eff) (Calabrese-Cardy 2005).",
    ),
    (
        EntanglementSaturationDensity,
        ["CalabreseCardy2005"],
        "Post-quench entanglement saturation density pi c/(6 beta_eff) (Calabrese-Cardy 2005).",
    ),
    (
        ConformalCasimirEnergy,
        ["Affleck1986"],
        "Universal Casimir ground-state energy -pi c/(6 L) (Affleck/BCN/Cardy 1986).",
    ),
    (
        ThermalEnergyDensity,
        ["Affleck1986"],
        "Leading thermal energy density pi c/(6 beta^2) (Affleck 1986).",
    ),
    (
        CFTThermalEntropyDensity,
        ["BloteCardyNightingale1986"],
        "Thermal entropy density pi c/(3 beta) (Bloete-Cardy-Nightingale 1986).",
    ),
    (
        CardyEntropy,
        String[],
        "Asymptotic high-energy entropy 2 pi sqrt(c E/6) (Cardy 1986).",
    ),
    (
        LogarithmicNegativity,
        String[],
        "Adjacent-interval logarithmic negativity (Calabrese-Cardy-Tonni 2012).",
    ),
)

for C in _CFT_CLASSES, (q, refs, note) in _CFT_PREDICTIONS
    register!(
        Universality{C},
        q,
        Infinite;
        method=:analytic,  # status=:universal derived by construction (register!)
        reliability=:high,
        references=refs,
        notes=note,
    )
end

# Register post-quench entanglement dynamics across all boundary conditions
for C in _CFT_CLASSES, BC in (Infinite, PBC, OBC)
    register!(
        Universality{C},
        VonNeumannEntropy{:quench},
        BC;
        method=:analytic,
        reliability=:high,
        references=["CalabreseCardy2005"],
        notes="Universal post-quench von Neumann entanglement entropy time evolution (Calabrese-Cardy 2005).",
    )
end

# Register Affleck-Ludwig boundary entropy log g at Infinite BC
for C in (:Ising, :XY, :Heisenberg)
    register!(
        Universality{C},
        BoundaryEntropy,
        Infinite;
        method=:analytic,
        reliability=:high,
        references=["AffleckLudwig1991", "Cardy1989"],
        notes="Affleck-Ludwig universal boundary entropy log g from the Cardy-state " *
              "modular S-matrix g_a = S_0a/sqrt(S_00) (Affleck-Ludwig 1991; Cardy 1989).",
    )
end
