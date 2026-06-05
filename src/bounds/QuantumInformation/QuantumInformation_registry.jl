# Registry rows for the QuantumInformation bounds domain.
# One row per theory regime (scheme); `canonical=true` marks the bare-fetch
# default (the quantum bound).  Per-scheme `references` fix whose bound each is.

# CHSH correlator S — upper bounds per regime (classical 2 / quantum 2√2 / PR 4).
@register(
    Bound{:QuantumInformation},
    CHSHBound,
    Infinite,
    scheme=:classical,
    method=:analytic,
    status=:bound,
    direction=:upper,
    canonical=false,
    reliability=:high,
    references=["CHSH1969"],
    tested_in="test/bounds/test_chsh_bound.jl",
    notes="Local-hidden-variable (classical) CHSH bound: S ≤ 2.",
)
@register(
    Bound{:QuantumInformation},
    CHSHBound,
    Infinite,
    scheme=:quantum,
    method=:analytic,
    status=:bound,
    direction=:upper,
    canonical=true,
    reliability=:high,
    references=["Tsirelson1980"],
    tested_in="test/bounds/test_chsh_bound.jl",
    notes="Tsirelson (quantum) CHSH bound: S ≤ 2√2, saturated by the optimal Bell state.",
)
@register(
    Bound{:QuantumInformation},
    CHSHBound,
    Infinite,
    scheme=:no_signalling,
    method=:analytic,
    status=:bound,
    direction=:upper,
    canonical=false,
    reliability=:high,
    references=["PopescuRohrlich1994"],
    tested_in="test/bounds/test_chsh_bound.jl",
    notes="No-signalling (Popescu-Rohrlich) CHSH bound: S ≤ 4.",
)

# Mermin 3-party operator |<M3>| — upper bounds per regime (classical 2 / quantum 4).
@register(
    Bound{:QuantumInformation},
    MerminGHZBound,
    Infinite,
    scheme=:classical,
    method=:analytic,
    status=:bound,
    direction=:upper,
    canonical=false,
    reliability=:high,
    references=["Mermin1990"],
    tested_in="test/bounds/test_mermin_ghz_bound.jl",
    notes="Local-realistic Mermin bound: |<M3>| ≤ 2.",
)
@register(
    Bound{:QuantumInformation},
    MerminGHZBound,
    Infinite,
    scheme=:quantum,
    method=:analytic,
    status=:bound,
    direction=:upper,
    canonical=true,
    reliability=:high,
    references=["Mermin1990"],
    tested_in="test/bounds/test_mermin_ghz_bound.jl",
    notes="Quantum Mermin bound: |<M3>| ≤ 4, saturated by the GHZ state.",
)

# Universal 1->2 qubit cloning fidelity — upper bound (Buzek-Hillery 1996).
@register(
    Bound{:QuantumInformation},
    OptimalCloningFidelity,
    Infinite,
    method=:analytic,
    status=:bound,
    direction=:upper,
    reliability=:high,
    references=["BuzekHillery1996"],
    tested_in="test/bounds/test_optimal_cloning_fidelity.jl",
    notes="Buzek-Hillery 1996: universal 1->2 qubit cloning fidelity F ≤ 5/6, saturated by the optimal cloner.",
)

# BB84 asymptotic secret-key rate R(e) = 1 - 2 H2(e) — achievable lower bound.
@register(
    Bound{:QuantumInformation},
    BB84KeyRate,
    Infinite,
    method=:analytic,
    status=:bound,
    direction=:lower,
    reliability=:high,
    references=["ShorPreskill2000"],
    tested_in="test/bounds/test_bb84_key_rate.jl",
    notes="Shor-Preskill 2000 BB84 secret-key rate R = 1 - 2 H2(qber): achievable rate, lower bound on key fraction; positive for qber < ~11%.",
)
