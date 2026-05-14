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

@register(
    CurieWeissIsing,
    CriticalExponents,
    Infinite,
    method=:delegation,
    reliability=:high,
    tested_in="test/models/classical/test_curie_weiss_ising.jl",
    references=["Landau 1937", "Stanley 1971"],
    notes="Delegated to MeanField — α=0, β=1/2, γ=1, δ=3, ν=1/2, η=0 (mean-field universality).",
)
