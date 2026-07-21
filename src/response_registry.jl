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

# ── U = ∂(βF)/∂β ──────────────────────────────────────────────────────
# The Gibbs–Helmholtz relation.  What is differentiated is the PRODUCT βF, not
# F alone — that is what `of` is for.  Independent of :entropy_response despite
# relating the same two potentials: this is the β-derivative of βF, that one is
# the T-derivative of F, and an implementation can satisfy one while breaking
# the other.
@response(
    :gibbs_helmholtz,
    relation = GibbsHelmholtz,
    derived = (dβF_dβ=∂(FreeEnergy, :β; of=(F, β) -> β * F),),
    sweep = (beta=[0.5, 1.0, 2.0],),
    finite_N = 6,
    exclusions = _THERMO_DERIVATIVE_EXCLUSIONS,
    rtol_floor = 1e-4,
    notes = "U = ∂(βF)/∂β — Gibbs–Helmholtz; equivalently U = -∂ln Z/∂β.",
)

# ── C = β² Var(E) ─────────────────────────────────────────────────────
# The energy-fluctuation route to the specific heat.  `var_E` is not a fetchable
# quantity, but it does not need to be: Var(E) = -∂⟨E⟩/∂β exactly, so `then`
# supplies it by negating the derivative of the energy.  With `Energy{:per_site}`
# the relation's `N` stays 1 — per-site C against per-site variance.
#
# This is a genuinely INDEPENDENT route to C: :specific_heat_from_entropy gets it
# from the entropy, this one from the energy.  A model that computes C by one
# formula and S or U by another is exactly what these two disagree on.
@response(
    :specific_heat_fdt,
    relation = SpecificHeatFDT,
    derived = (var_E=∂(Energy{:per_site}, :β; then=d -> -d),),
    sweep = (beta=[0.5, 1.0, 2.0],),
    finite_N = 6,
    exclusions = _THERMO_DERIVATIVE_EXCLUSIONS,
    rtol_floor = 1e-4,
    notes = "C = β² Var(E), Var(E) = -∂⟨E⟩/∂β — the fluctuation route, independent of C = T ∂S/∂T.",
)
