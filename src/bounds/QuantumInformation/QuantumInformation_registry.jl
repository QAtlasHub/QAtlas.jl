# Registry rows for the QuantumInformation bounds domain.
# (populates REGISTRY for Bound{:QuantumInformation}.)

@register(
    Bound{:QuantumInformation},
    CHSHBound,
    Infinite,
    method=:analytic,
    status=:bound,
    direction=:upper,
    reliability=:high,
    references=["CHSH1969", "Tsirelson1980", "PopescuRohrlich1994"],
    tested_in="test/bounds/test_chsh_bound.jl",
    notes="Upper bound on the CHSH correlator S; select whose bound with source=:bell/:tsirelson/:popescu_rohrlich (2 / 2√2 / 4). Quantum value 2√2 saturated by the optimal Bell state.",
)

@register(
    Bound{:QuantumInformation},
    MerminGHZBound,
    Infinite,
    method=:analytic,
    status=:bound,
    direction=:upper,
    reliability=:high,
    references=["Mermin1990"],
    tested_in="test/bounds/test_mermin_ghz_bound.jl",
    notes="Upper bound on the Mermin 3-party operator |<M3>|; source=:classical (2, local-realistic) / :mermin (4, GHZ-saturated quantum).",
)
