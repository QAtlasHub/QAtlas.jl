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

# ── v ≥ 0 ─────────────────────────────────────────────────────────────
# The first edge over a QUANTIFIED GROUP slot (`EachOf{AbstractVelocity}`,
# AbstractQAtlas design §8c).  Every member of that group is a propagation SPEED,
# and every relation consuming one requires it positive — `ξ = v/Δ` would give a
# negative correlation length, the CFT finite-size forms a negative gap.  A measured
# `v < 0` is a sign error in a dispersion derivative.
#
# The group is expanded by `_group_members` rather than `_family_members`: an
# abstract group's members are DIFFERENT quantities, most with no `component` to
# declare, and the component filter returns empty for them (#823).
#
# MEASURED: 3 group members are registered somewhere (Velocity{:luttinger},
# FermiVelocity, LiebRobinsonVelocity) across 4 hubs, and all 4 pass — no
# exclusions.  No sweep: these are T = 0 characteristic velocities.
@bound_edge(
    :velocity_positivity,
    inequality = VelocityPositivity,
    finite_N = 6,
    notes = "Every characteristic velocity is a propagation speed, so v ≥ 0; a negative value is a sign error in a dispersion derivative.",
)

# ── Δ ≥ 0 ─────────────────────────────────────────────────────────────
# The spectral gap is non-negative because `Δ = E₁ − E₀` and `E₀` is by definition
# the lowest level, so this is airtight rather than a modelling assumption — and it
# fails in a way no equality identity can catch: a negative gap means the reported
# "ground state" was not one (a variational state above a level the solver missed,
# or two levels ordered wrongly), which two equally-wrong quantities would still
# satisfy an identity between themselves.
#
# The widest edge in this file: 37 hubs register a MassGap.  MEASURED before
# declaring (see the exclusion below) — 36 of the 37 pass at default parameters,
# none returns a non-finite or negative value, and the one failure is a missing
# required kwarg rather than a physics disagreement.
#
# NO SWEEP, deliberately.  MassGap is a T = 0 spectral quantity; measured with
# `beta = 0.5` and `beta = 2.0` the results are identical to the no-kwarg call on
# every hub, so a beta sweep would triple the check count for the same numbers.
@bound_edge(
    :mass_gap_positivity,
    inequality = MassGapPositivity,
    finite_N = 6,
    exclusions = [
        GrossNeveu => "MassGap fetch takes a required cutoff kwarg `Λ` (UndefKeywordError without it); the generator's sweep grid is hub-generic and has no value to supply",
    ],
    notes = "Δ = E₁ − E₀ ≥ 0 by definition of the ground state; a negative value means the reported ground state was not one.",
)

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
