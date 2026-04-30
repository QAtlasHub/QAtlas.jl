# models/quantum/TFIM/TFIM_registry.jl — declarative implementation map.
#
# One `@register` line per natively-implemented (model, quantity, bc)
# triple.  Conversion fallbacks (e.g. `Energy(:per_site)` at OBC routed
# through the `Energy(:total)` native + `÷ N`) are *not* listed here:
# the registry tracks native implementations and the routing is
# automatic.  See `src/core/registry.jl` for the metadata schema.
#
# When you add a new fetch method to TFIM, add a sibling `@register`
# line here.  `test/core/test_registry.jl` will fail loudly if the
# registry and the dispatch table drift apart.

# ── Energy (granularity-aware) ─────────────────────────────────────────
@register(
    TFIM,
    Energy{:total},
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_thermal.jl",
    references=["Pfeuty 1970"],
    notes="Total ⟨H⟩(β) via the BdG spectrum; ground state when no β kwarg.",
)
@register(
    TFIM,
    Energy{:per_site},
    Infinite,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_thermal.jl",
    references=["Pfeuty 1970"],
    notes="Per-site ε(β) by QuadGK over the PBC dispersion Λ(k).",
)

# ── Spectrum / criticality ────────────────────────────────────────────
@register(
    TFIM,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/test_TFIM_massgap.jl",
    references=["Pfeuty 1970"],
    notes="Δ_∞(J,h) = 2|h - J| — closed-form Ising gap.",
)
@register(
    TFIM,
    MassGap,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_massgap.jl",
    references=["Pfeuty 1970"],
    notes="Smallest positive BdG eigenvalue of the OBC chain.",
)
@register(
    TFIM,
    CentralCharge,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/test_TFIM_central_charge.jl",
    references=["Belavin-Polyakov-Zamolodchikov 1984"],
    notes="c = 1/2 at the critical point (h = J), 0 otherwise.",
)

# ── Free-fermion thermal (per-site) — meta-defined in TFIM_thermal.jl ──
@register(
    TFIM,
    FreeEnergy,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_thermal.jl"
)
@register(
    TFIM,
    FreeEnergy,
    Infinite,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_thermal.jl"
)
@register(
    TFIM,
    ThermalEntropy,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_thermal.jl"
)
@register(
    TFIM,
    ThermalEntropy,
    Infinite,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_thermal.jl"
)
@register(
    TFIM,
    SpecificHeat,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_thermal.jl"
)
@register(
    TFIM,
    SpecificHeat,
    Infinite,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_thermal.jl"
)
@register(
    TFIM,
    MagnetizationX,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_thermal.jl"
)
@register(
    TFIM,
    MagnetizationX,
    Infinite,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_thermal.jl"
)
@register(
    TFIM,
    SusceptibilityXX,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_thermal.jl"
)
@register(
    TFIM,
    SusceptibilityXX,
    Infinite,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_thermal.jl"
)

# ── Local one-site observables (per-site index, not bulk-averaged) ────
@register(
    TFIM,
    EnergyLocal,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_local.jl"
)
@register(
    TFIM,
    MagnetizationXLocal,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_local.jl"
)
@register(
    TFIM,
    MagnetizationZLocal,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_local.jl"
)

# ── Z-axis dynamics + correlations (BdG time evolution) ───────────────
@register(
    TFIM,
    SusceptibilityZZ,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_dynamics.jl"
)
@register(
    TFIM,
    ZZStructureFactor,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_dynamics.jl"
)
@register(
    TFIM,
    ZZCorrelation{:static},
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_dynamics.jl"
)
@register(
    TFIM,
    ZZCorrelation{:dynamic},
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_dynamics.jl",
    notes="Single (i, j, t) point; for sweeps loop the kwargs.",
)
@register(
    TFIM,
    ZZCorrelation{:lightcone},
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_dynamics.jl",
    notes="C(r,t) lightcone slice for fixed center; takes a `times` vector.",
)
@register(
    TFIM,
    XXCorrelation{:dynamic},
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_dynamics.jl"
)

# ── Entanglement (T = 0; β kwarg defaults to Inf) ─────────────────────
@register(
    TFIM,
    VonNeumannEntropy,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_entanglement.jl",
    references=["Peschel 2003"],
    notes="Free-fermion correlation-matrix method; pass subsystem length ℓ.",
)

# ── PBC free-fermion thermal (per-site) — TFIM_pbc_thermal.jl ─────────
@register(
    TFIM,
    Energy{:per_site},
    PBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_pbc_thermal.jl",
    references=["Lieb-Schultz-Mattis 1961", "Sachdev 2011"],
    notes="Per-site ε(β) with parity-projected fermion sectors (NS + R).",
)
@register(
    TFIM,
    FreeEnergy,
    PBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_pbc_thermal.jl",
    references=["Lieb-Schultz-Mattis 1961"],
)
@register(
    TFIM,
    ThermalEntropy,
    PBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_pbc_thermal.jl",
)
@register(
    TFIM,
    SpecificHeat,
    PBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_pbc_thermal.jl",
)
@register(
    TFIM,
    MagnetizationX,
    PBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_pbc_thermal.jl",
)
@register(
    TFIM,
    SusceptibilityXX,
    PBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_pbc_thermal.jl",
)
@register(
    TFIM,
    MassGap,
    PBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_pbc_thermal.jl",
    notes="Smallest excitation across NS (two-mode flip) and R (one-mode flip) sectors.",
)

# ── Z-axis Infinite (TFIM_zaxis.jl) ──────────────────────────────────
@register(
    TFIM,
    MagnetizationZ,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/test_TFIM_zaxis.jl",
    references=["Pfeuty 1970"],
    notes="m_z = (1 - (h/J)²)^(1/8) for h < J, else 0 (T = 0 spontaneous).",
)
@register(
    TFIM,
    SpontaneousMagnetization,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/test_TFIM_zaxis.jl",
    references=["Pfeuty 1970"],
    notes="Same value as MagnetizationZ; order-parameter alias.",
)
@register(
    TFIM,
    CorrelationLength,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/test_TFIM_zaxis.jl",
    notes="ξ = 1/(2|h-J|) (gapped phase); Inf at criticality.",
)
@register(
    TFIM,
    SusceptibilityZZ,
    Infinite,
    method=:bdg,
    reliability=:medium,
    tested_in="test/models/test_TFIM_zaxis.jl",
    notes="OBC large-N proxy (N_proxy kwarg controls precision).",
)
@register(
    TFIM,
    ZZStructureFactor,
    Infinite,
    method=:bdg,
    reliability=:medium,
    tested_in="test/models/test_TFIM_zaxis.jl",
    notes="Static S_zz(q) via OBC large-N proxy of correlator Fourier sum.",
)

# ── ZZ correlator connected mode (OBC) ───────────────────────────────
@register(
    TFIM,
    ZZCorrelation{:connected},
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_zaxis.jl",
    notes="Connected = static for TFIM by Z₂ symmetry; explicit method for clarity.",
)

# ── Tier 2: XX static + connected via Pfaffian Wick contraction ─────
@register(
    TFIM,
    XXCorrelation{:static},
    OBC,
    method=:pfaffian,
    reliability=:high,
    tested_in="test/models/test_TFIM_xx_static.jl",
    references=["Lieb-Schultz-Mattis 1961", "Sachdev 2011"],
    notes="t=0 limit of dynamic XX correlator; real Pfaffian.",
)
@register(
    TFIM,
    XXCorrelation{:connected},
    OBC,
    method=:pfaffian,
    reliability=:high,
    tested_in="test/models/test_TFIM_xx_static.jl",
)
@register(
    TFIM,
    XXCorrelation{:static},
    Infinite,
    method=:pfaffian,
    reliability=:medium,
    tested_in="test/models/test_TFIM_xx_static.jl",
    notes="OBC large-N proxy (N_proxy kwarg).",
)
@register(
    TFIM,
    XXCorrelation{:connected},
    Infinite,
    method=:pfaffian,
    reliability=:medium,
    tested_in="test/models/test_TFIM_xx_static.jl",
)

# ── Tier 2: Renyi entropy at OBC via Peschel correlation matrix ──────
@register(
    TFIM,
    RenyiEntropy,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_renyi.jl",
    references=["Peschel 2003", "Calabrese-Cardy 2009"],
    notes="Free-fermion correlation-matrix Renyi α ≠ 1.",
)

# ── Tier 2: CC entanglement at Infinite ──────────────────────────────
@register(
    TFIM,
    VonNeumannEntropy,
    Infinite,
    method=:cft,
    reliability=:high,
    tested_in="test/models/test_TFIM_cft_entanglement.jl",
    references=["Calabrese-Cardy 2004", "Calabrese-Cardy 2009"],
    notes="Closed-form CC; T=0 critical/gapped + T>0 critical (gapped + T>0 errors).",
)
@register(
    TFIM,
    RenyiEntropy,
    Infinite,
    method=:cft,
    reliability=:high,
    tested_in="test/models/test_TFIM_cft_entanglement.jl",
    references=["Calabrese-Cardy 2009"],
    notes="CC Renyi: prefactor (c/6)(1 + 1/α). Same case coverage as VN.",
)

# ── Tier 2: dynamic structure factor at Infinite (proxy) ────────────
# Note: the static `ZZStructureFactor, Infinite` is registered via the
# router method in TFIM_infinite_dynamics.jl (when ω === nothing).
# The dynamic branch (ω::Real) is the new content.
