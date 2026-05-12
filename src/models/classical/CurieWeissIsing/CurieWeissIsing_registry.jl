# models/classical/CurieWeissIsing/CurieWeissIsing_registry.jl
#
# Declarative implementation map for the classical mean-field (Curie-Weiss)
# Ising model on the complete graph.  Schema documented in
# `src/core/registry.jl`.

@register(
    CurieWeissIsing,
    CriticalTemperature,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_curie_weiss_ising.jl",
    references=["Landau-Lifshitz §149"],
    notes="Mean-field T_c = J for J > 0 (k_B = 1); 0 otherwise.",
)

@register(
    CurieWeissIsing,
    SpontaneousMagnetization,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_curie_weiss_ising.jl",
    references=["Landau-Lifshitz §149"],
    notes="m = tanh(βJ m) self-consistency; bisection on g(m)=m-tanh(betaJ m) over [m_min, 1); 0 in paramagnetic phase.",
)
