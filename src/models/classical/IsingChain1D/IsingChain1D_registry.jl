# models/classical/IsingChain1D/IsingChain1D_registry.jl
#
# Declarative implementation map for the classical 1-D Ising chain
# (Ising 1925).  Schema documented in `src/core/registry.jl`.

@register(
    IsingChain1D,
    CriticalTemperature,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_ising_chain_1d.jl",
    references=["Ising1925"],
    notes="No finite-temperature phase transition in 1-D; T_c = 0.",
)

@register(
    IsingChain1D,
    FreeEnergy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_ising_chain_1d.jl",
    references=["Ising1925"],
    notes="f(β,h) = -β^{-1} log λ_+; at h=0 reduces to -β^{-1} log(2 cosh βJ).",
)

@register(
    IsingChain1D,
    CorrelationLength,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_ising_chain_1d.jl",
    references=["Ising1925"],
    notes="ξ(β,h) = 1/log(λ_+/λ_-); at h=0 reduces to 1/log(coth βJ).",
)

@register(
    IsingChain1D,
    Energy{:per_site},
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_ising_chain_1d.jl",
    references=["Ising1925"],
    notes="u(b,h=0) = -J tanh(bJ); h=0 only (textbook scope).",
)

@register(
    IsingChain1D,
    SpecificHeat,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_ising_chain_1d.jl",
    references=["Ising1925"],
    notes="c_v(b,h=0) = (bJ)^2 sech^2(bJ); h=0 only.",
)

@register(
    IsingChain1D,
    ThermalEntropy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_ising_chain_1d.jl",
    references=["Ising1925"],
    notes="s(b,h=0) = log(2 cosh bJ) - bJ tanh(bJ); h=0 only.",
)

@register(
    IsingChain1D,
    SusceptibilityZZ,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_ising_chain_1d.jl",
    references=["Ising1925", "Brush1967"],
    notes="chi(b,h=0) = b exp(2bJ) per site; h=0 only.",
)

@register(
    IsingChain1D,
    SpontaneousMagnetization,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_ising_chain_1d.jl",
    references=["Ising1925"],
    notes="m_spont = 0 for all T > 0 (no LRO in 1D Ising); independent of h.",
)
