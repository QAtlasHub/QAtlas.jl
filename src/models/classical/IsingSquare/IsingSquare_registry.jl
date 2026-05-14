# models/classical/IsingSquare/IsingSquare_registry.jl
#
# Declarative implementation map for the classical 2D Ising model.
# Schema documented in `src/core/registry.jl`.

# ── Onsager / Yang closed forms ──────────────────────────────────────
@register(
    IsingSquare,
    PartitionFunction,
    PBC,
    method=:transfer_matrix,
    reliability=:high,
    tested_in="test/util/classical_partition.jl",
    references=["Onsager 1944"],
    notes="Z = tr(T^Lx) on the symmetric transfer matrix; bond-counting differs at Lx,Ly=2.",
)
@register(
    IsingSquare,
    CriticalTemperature,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_ising2d_observables.jl",
    references=["Onsager 1944"],
    notes="T_c = 2J/log(1+√2).",
)
@register(
    IsingSquare,
    SpontaneousMagnetization,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_ising2d_observables.jl",
    references=["Yang 1952"],
    notes="M(T) = (1 - sinh⁻⁴(2βJ))^(1/8) for T < T_c.",
)

# ── Tier 4: finite-temperature thermodynamic potentials ──────────────
@register(
    IsingSquare,
    FreeEnergy,
    PBC,
    method=:transfer_matrix,
    reliability=:high,
    tested_in="test/models/test_IsingSquare_thermal.jl",
    notes="f = -log tr(T^Lx) / (β Lx Ly).",
)
@register(
    IsingSquare,
    Energy{:per_site},
    PBC,
    method=:central_diff,
    reliability=:medium,
    tested_in="test/models/test_IsingSquare_thermal.jl",
    notes="ε = -∂(log Z)/∂β / N via central diff (O(δ²) truncation).",
)
@register(
    IsingSquare,
    ThermalEntropy,
    PBC,
    method=:central_diff,
    reliability=:medium,
    tested_in="test/models/test_IsingSquare_thermal.jl",
    notes="s = β(ε - f).",
)
@register(
    IsingSquare,
    SpecificHeat,
    PBC,
    method=:central_diff,
    reliability=:medium,
    tested_in="test/models/test_IsingSquare_thermal.jl",
    notes="c_v = β² ∂²(log Z)/∂β² / N via 3-point stencil.",
)
@register(
    IsingSquare,
    FreeEnergy,
    Infinite,
    method=:onsager,
    reliability=:high,
    tested_in="test/models/test_IsingSquare_thermal.jl",
    references=["Onsager 1944"],
    notes="-βf = log(2) + (1/2π)∫log[(A+√(A²-B²))/2]dφ; bond-counting matches PBC limit.",
)
@register(
    IsingSquare,
    Energy{:per_site},
    Infinite,
    method=:central_diff,
    reliability=:high,
    tested_in="test/models/test_IsingSquare_thermal.jl",
    notes="ε = -∂(log Z/N)/∂β via central diff on Onsager log Z.",
)
@register(
    IsingSquare,
    ThermalEntropy,
    Infinite,
    method=:central_diff,
    reliability=:high,
    tested_in="test/models/test_IsingSquare_thermal.jl",
    notes="s = β(ε - f).",
)
@register(
    IsingSquare,
    SpecificHeat,
    Infinite,
    method=:central_diff,
    reliability=:medium,
    tested_in="test/models/test_IsingSquare_thermal.jl",
    notes="Diverges at T_c (ln(1+√2)/2 ≈ 0.4407); finite elsewhere.",
)

# ── Phase 2: critical exponents at T_c (universality delegation) ─────
@register(
    IsingSquare,
    CriticalExponents,
    Infinite,
    method=:delegation,
    reliability=:high,
    tested_in="test/models/classical/test_IsingSquare_critical.jl",
    references=["Onsager 1944"],
    notes="2D Ising Onsager exponents (β=1/8, γ=7/4, ν=1) delegated to Universality(:Ising) at d=2.",
)
