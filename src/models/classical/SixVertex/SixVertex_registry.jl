# models/classical/SixVertex/SixVertex_registry.jl
#
# Declarative implementation map for the classical six-vertex model on
# the square lattice (Lieb 1967a/b/c, Sutherland 1967, Baxter 1982).
# Schema documented in `src/core/registry.jl`.

@register(
    SixVertex,
    ResidualEntropy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_six_vertex.jl",
    references=["Lieb1967a", "Lieb1967b", "Sutherland1967", "Baxter1982"],
    notes="Configurational entropy density S = E - f. Square ice (a=b=c): S/N = (3/2) log(4/3) ≈ 0.4315231. FE phase: 0. Calculated off-diagonal using exact FreeEnergy and Energy.",
)

@register(
    SixVertex,
    FreeEnergy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_six_vertex.jl",
    references=["Lieb1967a", "Lieb1967b", "Lieb1967c", "Sutherland1967", "Baxter1982"],
    notes="FE phase (Δ > 1): f = -log max(a, b). Disordered phase (|Δ| <= 1): Lieb-Sutherland trigonometric integral. AFE phase (Δ < -1): Lieb elliptic sum.",
)

@register(
    SixVertex,
    Energy{:per_site},
    Infinite,
    method=:numerical,
    reliability=:high,
    tested_in="test/models/classical/test_six_vertex.jl",
    references=["Lieb1967a", "Lieb1967b", "Sutherland1967", "Baxter1982"],
    notes="Thermodynamic energy per site computed via central finite difference (h=1e-6) of the exact free energy.",
)

@register(
    SixVertex,
    Polarization,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_six_vertex.jl",
    references=["Baxter1982"],
    notes="Spontaneous bulk polarization density. FE phase: 1.0. Disordered phase: 0.0. AFE phase: staggered polarization given by the infinite product prod_{n=1}^inf tanh^2(n*lambda).",
)
