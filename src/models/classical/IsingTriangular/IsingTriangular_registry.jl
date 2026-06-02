# models/classical/IsingTriangular/IsingTriangular_registry.jl
#
# Declarative implementation map for the classical 2D Ising model on the
# triangular lattice (Wannier 1950 / Houtappel 1950).  Schema documented in
# `src/core/registry.jl`.

@register(
    IsingTriangular,
    CriticalTemperature,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_ising_triangular.jl",
    references=["Wannier1950", "Houtappel1950"],
    notes="AFM (J>0): T_c = 0 (frustrated, no order). FM (J<0): T_c = 4|J|/ln 3.",
)

@register(
    IsingTriangular,
    ResidualEntropy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_ising_triangular.jl",
    references=["Wannier1950"],
    notes="AFM (J>0): S/N = (2/π) ∫₀^{π/3} log(2 cos θ) dθ ≈ 0.3230659669 (QuadGK). FM (J<0): 0.",
)

@register(
    IsingTriangular,
    CriticalExponents,
    Infinite,
    method=:delegation,
    reliability=:high,
    tested_in="test/models/classical/test_ising_triangular.jl",
    references=["Onsager1944", "Houtappel1950"],
    notes="2D Ising universality (Onsager exponents) shared with IsingSquare; delegated to Universality(:Ising) d=2.",
)
