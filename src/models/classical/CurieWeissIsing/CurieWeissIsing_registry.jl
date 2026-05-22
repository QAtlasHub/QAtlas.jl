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
    tested_in="test/models/classical/test_curie_weiss_ising.jl",
    references=["Landau-Lifshitz §149"],
    notes="Mean-field T_c = J for J > 0 (k_B = 1); 0 otherwise. Zero-field reference value (no sharp transition at h ≠ 0).",
)

@register(
    CurieWeissIsing,
    SpontaneousMagnetization,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_curie_weiss_ising.jl",
    references=["Landau-Lifshitz §149"],
    notes="m*(β) = lim_{h→0⁺} m(β,J,h); positive root of m = tanh(βJm) by bisection; 0 in paramagnetic phase. Independent of model.h.",
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

@register(
    CurieWeissIsing,
    FreeEnergy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_curie_weiss_ising.jl",
    references=["Landau-Lifshitz §149"],
    notes="f(β;J,h) = Jm*²/2 − β⁻¹ log[2cosh(β(Jm*+h))] at the SCE solution; J≤0 reduces to single-spin -β⁻¹ log[2cosh(βh)].",
)

@register(
    CurieWeissIsing,
    Energy{:per_site},
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_curie_weiss_ising.jl",
    references=["Landau-Lifshitz §149"],
    notes="u(β;J,h) = -Jm*²/2 - h m*; J≤0 reduces to single-spin -h tanh(βh); T→0 saturation -J/2 - |h|.",
)

@register(
    CurieWeissIsing,
    ThermalEntropy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_curie_weiss_ising.jl",
    references=["Landau-Lifshitz §149"],
    notes="s(β;J,h) = log[2cosh(β(Jm*+h))] − β(Jm*+h)m* (Gibbs); bounded [0, log 2].",
)

@register(
    CurieWeissIsing,
    SpecificHeat,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_curie_weiss_ising.jl",
    references=["Landau-Lifshitz §149"],
    notes="c_v(β;J,h) = β²(Jm*+h)²(1-m*²)/[1 - βJ(1-m*²)]; J>0,h=0: 0 above T_c, jump 3/2 at T_c⁻; J≤0: (βh sech(βh))².",
)

@register(
    CurieWeissIsing,
    SusceptibilityZZ,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_curie_weiss_ising.jl",
    references=["Landau-Lifshitz §149"],
    notes="χ(β;J,h) = β(1-m*²)/[1 - βJ(1-m*²)] (∂m/∂h); Curie-Weiss law β/(1-βJ) at h=0, T>T_c; J≤0: β sech²(βh).",
)
