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
    tested_in="test/standalone/test_six_vertex.jl",
    references=["Lieb1967a"],
    notes="Square ice (a=b=c): S/N = (3/2) log(4/3) ≈ 0.4315231. FE phase (Δ > 1): 0. Generic disordered & AFE residual entropy deferred (issue #163 phase 2/3).",
)

@register(
    SixVertex,
    FreeEnergy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_six_vertex.jl",
    references=["Lieb1967a", "Lieb1967c", "Baxter 1982"],
    notes="FE phase (Δ > 1): f = -log max(a, b) (Lieb 1967c). Square-ice point (a = b = c): f = -(3/2) log(4/3) (Lieb 1967a). Generic disordered (off-diagonal) deferred to issue #163 phase 2; AFE phase (Δ < -1, Lieb 1967b elliptic) deferred to phase 3.",
)
