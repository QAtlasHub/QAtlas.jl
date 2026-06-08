# models/classical/DimerLattice/DimerLattice_registry.jl — declarative
# implementation map for the close-packed square-lattice dimer model
# (Kasteleyn-Temperley-Fisher).  See `src/core/registry.jl` for the schema.

@register(
    DimerLattice,
    PartitionFunction,
    OBC,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_dimer_lattice.jl",
    references=["Kasteleyn1961", "TemperleyFisher1961"],
    notes="Number of perfect matchings of the open Lx×Ly grid via the KTF product ∏(4cos²(jπ/(m+1))+4cos²(kπ/(n+1)))^{1/4}; 0 for odd Lx·Ly.",
)
@register(
    DimerLattice,
    ResidualEntropy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_dimer_lattice.jl",
    references=["Fisher1961"],
    notes="Entropy per site s = lim (ln Z)/N = G/π ≈ 0.29156 (Catalan/π); per dimer 2G/π.",
)
@register(
    DimerLattice,
    FreeEnergy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_dimer_lattice.jl",
    references=["Fisher1961"],
    notes="Free-energy density per site f = −s = −G/π (unit weights ⇒ zero internal energy).",
)
