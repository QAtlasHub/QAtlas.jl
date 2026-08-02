# bound_registry.jl — declared BOUND edges (core/bound.jl), the inequality
# sibling of identity_registry.jl.
#
# Each names an AbstractQAtlas `@bound` and lets the generator materialize
# it on every hub that implements the participants.  Nothing here restates the
# criterion: the slack, its sign convention, and the physics all live upstream.
#
# WHY THESE TWO.  Most of AbstractQAtlas's reachable bounds need something
# `fetch` does not return.  `FreeEnergyLegendre` is already covered by the
# `:gibbs` identity; the rest wait on a supplier for a derived input (`var_E`,
# `dS_dT`, the region entropies `S_A`/`S_AB`/…), which is the same gap
# identity_registry.jl's deletion criterion names for the derivative-form
# identities.  The eight universal bounds adopted in AbstractQAtlas 0.6.1 are
# unwirable for a different and permanent reason — their BOUNDED side is an
# untyped slot by design — and are allow-listed in test_abq_conformance.jl.
#
# The exact counts that used to stand here ("129 relations, 93 type-keyed, 22
# reachable") went stale within two releases and were never checked by anything.
# The reasoning above does not depend on them; test_abq_conformance.jl computes
# the current numbers and fails when one is unaccounted for.
#
# These are STABILITY bounds, and they fail differently from an equality: a
# sign error, a bad analytic continuation, or a mis-normalized thermal state
# shows up here even when two equally-wrong quantities would still satisfy an
# identity between themselves.

# Shared by the bound and response registries: the hubs whose finite-T
# thermodynamics the 1D generator cannot materialize, or that self-declare a
# validity window.  Same list :gibbs carries — kept in one place so the three
# edge kinds cannot drift apart on which hubs are honestly out of scope.
const _THERMO_DERIVATIVE_EXCLUSIONS = [
    (IsingSquare, PBC) => "2D PBC fetches take Lx/Ly kwargs, not bc.N — generator finite-size materialization is 1D-only",
    (KitaevHoneycomb, OBC) => "2D OBC fetches take Lx/Ly kwargs, not bc.N — generator finite-size materialization is 1D-only",
    IsingTriangular => "default J > 0 is the frustrated AFM branch (no Houtappel closed form); finite-T requires J < 0",
    AKLT1D => "finite-β canonical thermodynamics supports β = ∞ only (HTSE is a separate :approx scheme, #506)",
    Heisenberg1D => "Infinite-BC thermodynamics is a c=1 CFT low-T expansion valid only for β > 5/J; it warns and returns NaN on this sweep (#521 Path A will replace it)",
    HaldaneShastry => "Infinite-BC thermodynamics returns NaN outside its low-T validity window, same CFT-expansion guard as Heisenberg1D",
]

# ── C_v ≥ 0 ───────────────────────────────────────────────────────────
# Thermodynamic stability: the specific heat is a variance (β²·Var(E)) and so
# cannot be negative for any equilibrium state at any β.  Exact at every N, so
# a small finite_N loses no coverage — same argument as :gibbs.
@bound_edge(
    :specific_heat_positivity,
    inequality = SpecificHeatPositivity,
    sweep = (beta=[0.5, 1.0, 2.0],),
    finite_N = 6,
    exclusions = _THERMO_DERIVATIVE_EXCLUSIONS,
    notes = "C_v = β² Var(E) ≥ 0 — thermodynamic stability; holds at every N and β.",
)

# ── χ_T ≥ 0 ───────────────────────────────────────────────────────────
# The isothermal susceptibility is likewise a variance (β·Var(M)), so it is
# non-negative for every axis pair.  `SusceptibilityPositivity` types its slot
# as the parametric FAMILY `Susceptibility`, so the generator expands it to one
# check per concrete axis pair the hub implements.
@bound_edge(
    :susceptibility_positivity,
    inequality = SusceptibilityPositivity,
    sweep = (beta=[0.5, 1.0],),
    finite_N = 6,
    exclusions = _THERMO_DERIVATIVE_EXCLUSIONS,
    notes = "χ_αα = β Var(M_α) ≥ 0 for every axis pair — thermodynamic stability.",
)
