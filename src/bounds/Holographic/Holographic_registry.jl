# Registry rows for the Holographic bounds domain.

@register(
    Bound{:Holographic},
    BekensteinBound,
    Infinite,
    method=:analytic,
    status=:bound,
    direction=:upper,
    reliability=:high,
    references=["Bekenstein1981"],
    tested_in="test/bounds/test_bekenstein_bound.jl",
    notes="Bekenstein 1981 entropy bound: S ≤ 2π R E, saturated (up to O(1)) by black holes.",
)
