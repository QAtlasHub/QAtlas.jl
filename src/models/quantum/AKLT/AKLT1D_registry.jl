# models/quantum/AKLT/AKLT1D_registry.jl —
# declarative implementation map for AKLT1D (S=1 BLBQ at the AKLT point).
#
# Closed-form analytical rows have `reliability=:high` because they
# follow from the AKLT 1988 projector construction and the
# Kennedy–Tasaki 1992 hidden-symmetry argument.  The Haldane gap row is
# `:medium` because it is an García-Saez–Murg–Verstraete 2013 numerical-exact value
# with no closed form.

# ── Infinite analytical rows (closed form) ───────────────────────────
@register(
    AKLT1D,
    GroundStateEnergyDensity,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_aklt.jl",
    references=["AKLT1988"],
    notes="Closed form e₀ = -2J/3 from the bond-projector decomposition of H.",
)
@register(
    AKLT1D,
    Energy{:per_site},
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_aklt.jl",
    references=["AKLT1988"],
    notes="Same analytic e₀ = -2J/3 routed through Energy(:per_site).",
)
@register(
    AKLT1D,
    CorrelationLength,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_aklt.jl",
    references=["AKLT1988"],
    notes="Closed form ξ = 1/log 3 ≈ 0.910 from VBS transfer matrix.",
)
@register(
    AKLT1D,
    StringOrderParameter,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_aklt.jl",
    references=["AKLT1988", "KennedyTasaki1992"],
    notes="Closed form O_str = 4/9; detects hidden Z₂×Z₂ in the Haldane phase.",
)

# ── Infinite numerical-exact (no closed form) ────────────────────────
@register(
    AKLT1D,
    MassGap,
    Infinite,
    method=:literature_value,
    reliability=:medium,
    tested_in="test/standalone/test_aklt.jl",
    references=["GarciaSaez2013"],
    notes="Haldane gap Δ ≈ 0.350 J; DMRG numerical-exact, no closed form.",
)

# ── OBC dense ED (cap N ≤ 8 from _MAX_ED_SITES_S1) ───────────────────
@register(
    AKLT1D,
    ExactSpectrum,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/standalone/test_aklt.jl",
    notes="Full sorted spectrum from 3^N dense ED; N ≤ 8 (3^8 = 6561).",
)

# ── VBS ground-state correlations (Infinite, closed form) ────────────
@register(
    AKLT1D,
    ZZCorrelation{:static},
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_aklt.jl",
    references=["AKLT1988"],
    notes="Exact VBS ⟨Sᶻ₀Sᶻ_r⟩ = (-1)^r (4/3) 3^{-|r|} (r≠0), 2/3 (r=0); J-independent. ED in tests verifies finite-N convergence.",
)
@register(
    AKLT1D,
    ZZStructureFactor,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_aklt.jl",
    references=["Arovas1988"],
    notes="Exact static structure factor S_zz(q) = 2(1-cos q)/(5+3 cos q); S(0)=0, S(π)=2; J-independent.",
)

# ── β = ∞ (T = 0) thermodynamic limits (closed form from bond-projector
# decomposition, AKLT 1988).  Finite β throws DomainError — no analytic
# reduction known for the AKLT chain (not Bethe-ansatz integrable).
# OBC × SusceptibilityZZ is intentionally not registered: it diverges at
# β = ∞ (edge-mode Curie tail); the fetch method exists only to throw a
# descriptive DomainError.
@register(
    AKLT1D,
    FreeEnergy,
    OBC,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_aklt_thermal_limits.jl",
    references=["AKLT1988"],
    notes="β=∞ only: f_OBC(∞) = -(2J/3)(N-1)/N from bond-projector E_GS_OBC = -(2J/3)(N-1). Finite β throws DomainError.",
)
@register(
    AKLT1D,
    ThermalEntropy,
    OBC,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_aklt_thermal_limits.jl",
    references=["AKLT1988"],
    notes="β=∞ only: s_OBC(∞) = log(4)/N from 4-fold edge-mode GS degeneracy (two free spin-½). Finite β throws DomainError.",
)
@register(
    AKLT1D,
    SpecificHeat,
    OBC,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_aklt_thermal_limits.jl",
    references=["AKLT1988"],
    notes="β=∞ only: c_OBC(∞) = 0 (pure GS manifold, zero energy variance). Finite β throws DomainError.",
)
@register(
    AKLT1D,
    FreeEnergy,
    PBC,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_aklt_thermal_limits.jl",
    references=["AKLT1988"],
    notes="β=∞ only: f_PBC(∞) = -2J/3 (unique VBS GS, all N bond projectors annihilated). Finite β throws DomainError.",
)
@register(
    AKLT1D,
    ThermalEntropy,
    PBC,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_aklt_thermal_limits.jl",
    references=["AKLT1988"],
    notes="β=∞ only: s_PBC(∞) = 0 (unique gapped GS). Finite β throws DomainError.",
)
@register(
    AKLT1D,
    SpecificHeat,
    PBC,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_aklt_thermal_limits.jl",
    references=["AKLT1988"],
    notes="β=∞ only: c_PBC(∞) = 0. Finite β throws DomainError.",
)
@register(
    AKLT1D,
    SusceptibilityZZ,
    PBC,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_aklt_thermal_limits.jl",
    references=["AKLT1988", "GarciaSaez2013"],
    notes="β=∞ only: χ_PBC(∞) = 0 (Haldane gap exponential suppression). Finite β throws DomainError.",
)
@register(
    AKLT1D,
    FreeEnergy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_aklt_thermal_limits.jl",
    references=["AKLT1988"],
    notes="β=∞ only: f(∞) = -2J/3, matches GroundStateEnergyDensity at Infinite. Finite β throws DomainError.",
)
@register(
    AKLT1D,
    ThermalEntropy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_aklt_thermal_limits.jl",
    references=["AKLT1988"],
    notes="β=∞ only: s(∞) = 0 (unique bulk GS in the Haldane phase). Finite β throws DomainError.",
)
@register(
    AKLT1D,
    SpecificHeat,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_aklt_thermal_limits.jl",
    references=["AKLT1988"],
    notes="β=∞ only: c(∞) = 0. Finite β throws DomainError.",
)
@register(
    AKLT1D,
    SusceptibilityZZ,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_aklt_thermal_limits.jl",
    references=["AKLT1988", "GarciaSaez2013"],
    notes="β=∞ only: χ(∞) = 0 (Haldane gap suppression). Finite β throws DomainError.",
)
