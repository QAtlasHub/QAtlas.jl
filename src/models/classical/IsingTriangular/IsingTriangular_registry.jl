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

# ── Ferromagnetic (J<0) finite-T thermodynamics — Houtappel 1950 ──────
@register(
    IsingTriangular,
    FreeEnergy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_ising_triangular_thermo.jl",
    references=["Houtappel1950"],
    notes="FM (J<0): -βf = ln2 + (1/8π²)∫∫ ln[cosh³2K+sinh³2K-sinh2K(cosθ+cosφ+cos(θ+φ))], K=β|J|.",
)

@register(
    IsingTriangular,
    Energy{:per_site},
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_ising_triangular_thermo.jl",
    references=["Houtappel1950"],
    notes="FM (J<0): ε = -∂(logZ/N)/∂β (central diff of Houtappel f); ε → -3|J| at T→0.",
)

@register(
    IsingTriangular,
    SpecificHeat,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_ising_triangular_thermo.jl",
    references=["Houtappel1950"],
    notes="FM (J<0): c_v = β²∂²(logZ/N)/∂β²; singular at T_c=4|J|/ln3 (2D-Ising universality).",
)

@register(
    IsingTriangular,
    ThermalEntropy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_ising_triangular_thermo.jl",
    references=["Houtappel1950"],
    notes="FM (J<0): s = β(ε-f); bounded 0 (T→0) … ln2 (T→∞).",
)

@register(
    IsingTriangular,
    SpontaneousMagnetization,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_ising_triangular_thermo.jl",
    references=["Houtappel1950", "Baxter1982"],
    notes="FM (J<0): Potts–Domb M=[1-16x³/((1-x)³(1+3x))]^{1/8}, x=e^{-4β|J|}, 0 above T_c. AFM (J>0): 0 (frustrated).",
)

@register(
    IsingTriangular,
    ZZCorrelation{:static},
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_ising_triangular_thermo.jl",
    references=["Wannier1950"],
    notes="AFM (J>0) T=0 nearest-neighbour ⟨σσ⟩ = -1/3 (every triangle: Σσσ=-1 ⇒ -1/3 per bond). r>1 (Stephenson 1964) deferred.",
)
