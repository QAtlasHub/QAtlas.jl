# Registry rows for the Dynamics bounds domain.

@register(
    Bound{:Dynamics},
    ChaosBound,
    Infinite,
    method=:analytic,
    status=:bound,
    direction=:upper,
    reliability=:high,
    references=["MaldacenaShenkerStanford2016"],
    tested_in="test/bounds/test_chaos_bound.jl",
    notes="MSS 2016 chaos bound: Lyapunov λ_L ≤ 2π/β, saturated by holographic / large-N SYK.",
)
