# response_registry.jl — declared RESPONSE edges (core/response.jl).
#
# These are the AbstractQAtlas relations that were reachable on QAtlas hubs but
# unusable, because one slot is a derivative rather than a fetchable value.  With
# core/derivative.jl supplying it they become ordinary generated checks.
#
# Both are exact thermodynamic identities that hold at every N and β, so a small
# finite_N loses no coverage — the same argument :gibbs makes.

# ── S = -∂F/∂T ────────────────────────────────────────────────────────
@response(
    :entropy_response,
    relation = EntropyResponse,
    derived = (dF_dT=∂(FreeEnergy, :T),),
    sweep = (beta=[0.5, 1.0, 2.0],),
    finite_N = 6,
    exclusions = _THERMO_DERIVATIVE_EXCLUSIONS,
    notes = "S = -∂F/∂T — the Maxwell relation between entropy and the free energy.",
)

# ── C = T ∂S/∂T ───────────────────────────────────────────────────────
@response(
    :specific_heat_from_entropy,
    relation = SpecificHeatFromEntropy,
    derived = (dS_dT=∂(ThermalEntropy, :T),),
    sweep = (beta=[0.5, 1.0, 2.0],),
    finite_N = 6,
    exclusions = _THERMO_DERIVATIVE_EXCLUSIONS,
    # Onsager's C and S are quadrature evaluations: measured, they reproduce
    # C = T ∂S/∂T to 2.5e-5 on IsingSquare with AD and with finite differences
    # alike.  That is the accuracy of the FETCHES, so AD's 1e-6 would report a
    # quadrature residue as a physics failure.
    rtol_floor = 1e-4,
    notes = "C = T ∂S/∂T — the caloric definition, independent of the C = β² Var(E) route.",
)
