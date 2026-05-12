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
    references=["Ising 1925"],
    notes="No finite-temperature phase transition in 1-D; T_c = 0.",
)

@register(
    IsingChain1D,
    FreeEnergy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_ising_chain_1d.jl",
    references=["Ising 1925"],
    notes="f(β,h) = -β^{-1} log λ_+; at h=0 reduces to -β^{-1} log(2 cosh βJ).",
)

@register(
    IsingChain1D,
    CorrelationLength,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_ising_chain_1d.jl",
    references=["Ising 1925"],
    notes="ξ(β,h) = 1/log(λ_+/λ_-); at h=0 reduces to 1/log(coth βJ).",
)
