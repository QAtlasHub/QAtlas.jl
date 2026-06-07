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

@register(
    Bound{:Dynamics},
    QuantumSpeedLimit,
    Infinite,
    scheme=:margolus_levitin,
    method=:analytic,
    status=:bound,
    direction=:lower,
    canonical=true,
    reliability=:high,
    references=["MargolusLevitin1998"],
    tested_in="test/bounds/test_quantum_speed_limit.jl",
    notes="Margolus-Levitin 1998 quantum speed limit: orthogonalization time τ ≥ π/(2E).",
)
