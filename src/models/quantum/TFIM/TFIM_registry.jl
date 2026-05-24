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
    MagnetizationXLocal{:equilibrium},
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/models/test_TFIM_local.jl"
)
@register(
    TFIM,
    MagnetizationXLocal{:quench},
    OBC,
    method=:majorana_evolution,
    reliability=:high,
    tested_in="test/standalone/test_tfim_sigma_x_quench.jl",
    references=["Barouch-McCoy-Dresden 1970", "Calabrese-Essler-Fagotti 2012"],
    notes="Sudden h_0 -> h_f quench; Sigma(t) = R(t) Sigma_0 R(t)^T; sigma^x_i(t) = Sigma(t)[2i-1, 2i].",
)
@register(
    TFIM,
    MagnetizationXLocal{:quench},
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_tfim_sigma_x_quench.jl",
    references=["Barouch-McCoy-Dresden 1970", "Calabrese-Essler-Fagotti 2012"],
    notes="Closed-form k-integral over the Bogoliubov angles theta_k(h_0,f); QuadGK rtol=1e-12.",
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
    VonNeumannEntropy{:equilibrium},
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
    tested_in="test/verification/test_tfim_fdt.jl",
    notes="OBC large-N proxy.  Static (ω = nothing) → uniform χ_zz(β); dynamic (ω::Real, q required) → χ''_zz(q,ω;β) via Kubo commutator.  N_proxy kwarg controls precision.",
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
    VonNeumannEntropy{:equilibrium},
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

# ── Tier 3: σʸ correlators + MagnetizationY/SusceptibilityYY OBC ────
# (defined in TFIM_yy.jl; closes the YY gap left by PR #130 Tier 2)
@register(
    TFIM,
    YYCorrelation{:static},
    OBC,
    method=:pfaffian,
    reliability=:high,
    tested_in="test/models/test_TFIM_yy.jl",
    notes="σʸ_i = -(-i)^{i-1} γ_1 … γ_{2i-2} γ_{2i}; same Pfaffian machinery as σᶻ.",
)
@register(
    TFIM,
    YYCorrelation{:connected},
    OBC,
    method=:pfaffian,
    reliability=:high,
    tested_in="test/models/test_TFIM_yy.jl",
    notes="Connected = static for TFIM since ⟨σʸ⟩ = 0 by parity (odd Majorana product).",
)
@register(
    TFIM,
    YYCorrelation{:dynamic},
    OBC,
    method=:pfaffian,
    reliability=:high,
    tested_in="test/models/test_TFIM_yy.jl",
)
@register(
    TFIM,
    MagnetizationY,
    OBC,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/test_TFIM_yy.jl",
    notes="Identically 0 in any Gaussian state (odd-Majorana product).",
)
@register(
    TFIM,
    SusceptibilityYY,
    OBC,
    method=:pfaffian,
    reliability=:high,
    tested_in="test/models/test_TFIM_yy.jl",
    notes="Per-site β·Var(M_y)/N via Wick contraction over O(N²) pairs.",
)

# ── Tier 3: XX / YY static structure factors at OBC + Infinite proxy ─
@register(
    TFIM,
    XXStructureFactor,
    OBC,
    method=:pfaffian,
    reliability=:high,
    tested_in="test/models/test_TFIM_xx_yy_structure_factor.jl",
    notes="(1/N) Σ_{ij} e^{-iq(i-j)} ⟨σˣᵢ σˣⱼ⟩ from t=0 Pfaffian correlator.",
)
@register(
    TFIM,
    YYStructureFactor,
    OBC,
    method=:pfaffian,
    reliability=:high,
    tested_in="test/models/test_TFIM_xx_yy_structure_factor.jl",
)
@register(
    TFIM,
    XXStructureFactor,
    Infinite,
    method=:pfaffian,
    reliability=:medium,
    tested_in="test/models/test_TFIM_xx_yy_structure_factor.jl",
    notes="OBC large-N proxy (N_proxy kwarg).",
)
@register(
    TFIM,
    YYStructureFactor,
    Infinite,
    method=:pfaffian,
    reliability=:medium,
    tested_in="test/models/test_TFIM_xx_yy_structure_factor.jl",
)

# ── Quench dynamics: Loschmidt echo + DQPT rate function ─────────────
@register(
    TFIM,
    LoschmidtEcho{:amplitude},
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/standalone/test_tfim_loschmidt.jl",
    references=[
        "Heyl-Polkovnikov-Kehrein PRL 110, 135704 (2013)",
        "Heyl Rep. Prog. Phys. 81, 054001 (2018)",
    ],
    notes="L(t) = ∏_n |cos²θ_n + sin²θ_n e^{-2iΛ_n t}|² via OBC BdG diagonalisation of H_0, H_f.",
)
@register(
    TFIM,
    LoschmidtEcho{:rate},
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/standalone/test_tfim_loschmidt.jl",
    references=["Heyl-Polkovnikov-Kehrein PRL 110, 135704 (2013)"],
    notes="λ(t) = -log L(t) / N from the OBC BdG product.",
)
@register(
    TFIM,
    LoschmidtEcho{:rate},
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_tfim_loschmidt.jl",
    references=[
        "Heyl-Polkovnikov-Kehrein PRL 110, 135704 (2013)",
        "Heyl Rep. Prog. Phys. 81, 054001 (2018)",
    ],
    notes="λ(t) = -(1/2π) ∫₀^π log|cos²Δθ_k + sin²Δθ_k e^{-2iΛ_k(h_f) t}|² dk via QuadGK.",
)
# ── GGE stationary values for h-quench (TFIM_gge.jl) ─────────────────
@register(
    TFIM,
    GGEValue{Energy{:per_site}},
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_tfim_gge.jl",
    references=[
        "Rigol et al. PRL 98 (2007)", "Calabrese-Essler-Fagotti J. Stat. Mech. (2012)"
    ],
    notes="Per-site ε_GGE via QuadGK over the post-quench dispersion with Bogoliubov-mismatch occupations n_k = sin²(θ_k(h₀)−θ_k(h_f)).",
)
@register(
    TFIM,
    GGEValue{MagnetizationX},
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_tfim_gge.jl",
    references=[
        "Rigol et al. PRL 98 (2007)", "Calabrese-Essler-Fagotti J. Stat. Mech. (2012)"
    ],
    notes="⟨σˣ⟩_GGE via QuadGK over (h_f − J cos k)/Λ_k(h_f) · (1 − 2 n_k).",
)

# ── Quench entanglement entropy (issue #144) ─────────────────────────
@register(
    TFIM,
    VonNeumannEntropy{:quench},
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/standalone/test_tfim_quench_entanglement.jl",
    references=["Calabrese-Cardy 2005", "Peschel 2003"],
    notes="S(ℓ, t) after a global quench from initial::TFIM ground state; Peschel on time-evolved Σ(t) = R(t) Σ_0 R(t)ᵀ.",
)

# ── Fidelity susceptibility (BdG analytical, issue #147) ─────────────
@register(
    TFIM,
    FidelitySusceptibility,
    OBC,
    method=:bdg,
    reliability=:high,
    tested_in="test/standalone/test_tfim_fidelity_susceptibility.jl",
    references=["Gu IJMPB 24 4371 (2010)", "Damski PRB 87 165101 (2013)"],
    notes="χ_F = Σ_{p<q} 4 X_{pq}² / (Λ_p+Λ_q)² from Bogoliubov amplitudes.",
)
@register(
    TFIM,
    FidelitySusceptibility,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_tfim_fidelity_susceptibility.jl",
    references=["Gu IJMPB 24 4371 (2010)", "Damski PRB 87 165101 (2013)"],
    notes="χ_F/L = 1/(16(J²-h²)) (h<J), J²/(16h²(h²-J²)) (h>J); QuadGK, divergent at h=J.",
)

# ── Critical exponents (Phase 2 delegation to 2D Ising universality) ──
@register(
    TFIM,
    CriticalExponents,
    Infinite,
    method=:delegation,
    reliability=:high,
    tested_in="test/models/quantum/TFIM/test_tfim_critical_exponents.jl",
    references=["Onsager 1944", "Pfeuty 1970"],
    notes="2D-Ising Onsager exponents (β=1/8, γ=7/4, ν=1) via TFIM↔2D-Ising mapping.",
)

# ── Loschmidt rate function (quench dynamics, infinite chain) ─────────
@register(
    TFIM,
    LoschmidtRateFunction,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/TFIM/test_tfim_loschmidt.jl",
    references=["Heyl Polkovnikov Kehrein PRL 110 135704 (2013)"],
    notes="λ(t) = -lim_{L→∞} (1/L) log |⟨ψ₀|e^{-iH_f t}|ψ₀⟩|² for TFIM h-quench; closed-form via Bogoliubov mode amplitudes.",
)
