# derivation_registry.jl — the scaling-plane catalog (core/derivation.jl).
#
# One declaration per hub with a `CriticalExponents` method, swept or refused.
# `check_derivation_coverage()` enforces that: a hub missing from here is a
# table the AbstractQAtlas `:scaling` algebra never sees, and silence about it
# would read exactly like coverage.
#
# `derived_from` is transcribed from each table's own provenance header, and it
# is the difference between a check and a tautology: where the header says a
# value follows FROM a scaling relation, that route is skipped with the reason
# rather than counted as agreement.

# ── Ising ────────────────────────────────────────────────────────────
# d=2 provenance (src/universalities/Ising2D/Ising2D.jl header): α Onsager 1944,
# β Yang 1952, γ Fisher 1964, η Kadanoff 1966, ν den Nijs 1979 — each
# individually sourced.  δ = 15 is the one the header marks as obtained from the
# scaling relation, so Widom cannot judge it.
@exponent_sweep(
    Universality{:Ising},
    Infinite,
    sweep = (d=[2],),
    derived_from = [:δ => [:Widom]],
    references = ["Onsager1944", "Yang1952"],
)

# d=3 is one bootstrap table (Kos-Poland-Simmons-Duffin-Vichi 2016, Table 2) and
# d=4 delegates to the mean-field table; neither header derives a value from a
# scaling relation, so every route judges.
@exponent_sweep(
    Universality{:Ising}, Infinite, sweep = (d=[3, 4],), references = ["Kos2016"],
)

# ── Percolation ──────────────────────────────────────────────────────
# d=2 provenance (src/universalities/Percolation/Percolation.jl header): β, ν, η
# from the Coulomb gas (Nienhuis 1982, den Nijs 1979); γ, δ and α are each
# marked as following from Fisher, Widom and Rushbrooke respectively.
@exponent_sweep(
    Universality{:Percolation},
    Infinite,
    sweep = (d=[2],),
    derived_from = [:γ => [:Fisher], :δ => [:Widom], :α => [:Rushbrooke]],
)

# d=3 is the Monte-Carlo table of Wang-Zhou-Zhang-Garoni-Deng 2013; d≥6 is the
# mean-field row above the upper critical dimension.
@exponent_sweep(Universality{:Percolation}, Infinite, sweep = (d=[3, 6],),)

# ── Potts ────────────────────────────────────────────────────────────
# Both tables come from the Coulomb-gas formulas in the coupling g, not from the
# scaling relations, so every route judges.
@exponent_sweep(Universality{:Potts3}, Infinite, sweep = (d=[2],))
@exponent_sweep(Universality{:Potts4}, Infinite, sweep = (d=[2],))

# ── O(n) ─────────────────────────────────────────────────────────────
# d=2 is the BKT point and returns η alone, which the generator refuses as too
# few to close — a visible skip, which is the honest report for a transition with
# no power-law exponents.  d=3 is one bootstrap table, d≥4 mean-field.
@exponent_sweep(
    Universality{:XY}, Infinite, sweep = (d=[2, 3, 4],), references = ["Chester2020"]
)
@exponent_sweep(
    Universality{:Heisenberg}, Infinite, sweep = (d=[3, 4],), references = ["Chester2021"],
)

# ── Mean field ───────────────────────────────────────────────────────
# `d` is not a fetch kwarg here and the table is d-independent, but hyperscaling
# is not: `2 − α = dν` holds at the upper critical dimension and nowhere else, so
# the dimension is pinned to d_c = 4 rather than left to a sweep.
@exponent_sweep(
    Universality{:MeanField},
    Infinite,
    sweep = (;),
    dimension = 4,
    references = ["LandauLifshitz1980"],
)

# ── Refusals ─────────────────────────────────────────────────────────
# A `consistency_report` on a NamedTuple keys on formula LETTERS, so a table
# whose letters mean something else must be kept out rather than fed in and
# reported as a contradiction.

@refuse_exponents(
    KPZ1D,
    Infinite,
    reason = "KPZ1D/CriticalExponents returns GrowthExponents: its α is the \
               roughness exponent and its β the growth exponent, not the \
               specific-heat and order-parameter exponents the :scaling \
               relations are written on. Same letters, different quantities.",
)

@refuse_exponents(
    Universality{:IsingSDRG},
    Infinite,
    reason = "the infinite-randomness table carries β and ν beside ψ, φ and \
               x_m, and its d=2 is QAtlas's Euclidean convention while the \
               relations' d is the spatial dimension (d=1 for the chain). Two \
               different d and a different fixed point; ActivatedExponent and \
               the :scaling algebra are not the same statement.",
)

@refuse_exponents(
    Ising2D,
    Infinite,
    reason = "pinned-d alias of Universality(:Ising) at d=2, which is swept \
               above; the same table under a second name.",
)

@refuse_exponents(
    MeanField,
    Infinite,
    reason = "legacy alias of Universality(:MeanField), which is swept above; \
               the same table under a second name.",
)
