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

@register(
    Bound{:Dynamics},
    QuantumSpeedLimit,
    Infinite,
    scheme=:mandelstam_tamm,
    method=:analytic,
    status=:bound,
    direction=:lower,
    canonical=false,
    reliability=:high,
    references=["MandelstamTamm1945"],
    tested_in="test/bounds/test_quantum_speed_limit.jl",
    notes="Mandelstam-Tamm 1945 quantum speed limit: orthogonalization time τ ≥ π/(2ΔE), with ΔE the energy uncertainty.",
)

@register(
    Bound{:Dynamics},
    ScramblingTime,
    Infinite,
    method=:analytic,
    status=:bound,
    direction=:lower,
    reliability=:high,
    references=["SekinoSusskind2008"],
    tested_in="test/bounds/test_scrambling_time.jl",
    notes="Sekino-Susskind 2008 fast-scrambling conjecture: t_* = (β/2π) log N, a lower bound saturated by black holes.",
)
