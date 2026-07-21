# models/quantum/XXZ/XXZ_registry.jl — declarative implementation map.
#
# One `@register` line per natively-implemented (model, quantity, bc)
# triple.  See `src/core/registry.jl` for the metadata schema and
# `src/models/quantum/TFIM/TFIM_registry.jl` for the canonical example.
#
# `:dense_ed` rows are tagged `:high` reliability — finite-N ED is exact
# for any N ≤ `_MAX_ED_SITES`, so the only failure mode is a cap miss
# (caught by an explicit `ArgumentError`).

# ── Energy (granularity-aware) ─────────────────────────────────────────
@register(
    XXZ1D,
    Energy{:total},
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_thermal.jl",
    references=["YangYang1966", "Takahashi1999"],
    notes="Total ⟨H⟩(β) by dense ED of the 2^N × 2^N XXZ Hamiltonian.",
)
@register(
    XXZ1D,
    Energy{:per_site},
    Infinite,
    method=:bethe_ansatz,
    reliability=:high,
    tested_in="test/models/test_XXZ1D.jl",
    references=["Hulthen1938", "YangYang1966"],
    notes="Closed form at Δ ∈ {-1, 0, 1}; Yang-Yang single integral via QuadGK for general -1 < Δ < 1; |Δ| > 1 (gapped) deferred.",
)

# ── Spectrum / criticality ────────────────────────────────────────────
@register(
    XXZ1D,
    CentralCharge,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/test_XXZ1D.jl",
    references=["Giamarchi2003"],
    notes="c = 1 in the critical regime -1 < Δ < 1.",
)
@register(
    XXZ1D,
    LuttingerParameter,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/test_XXZ1D.jl",
    references=["Giamarchi2003"],
    notes="K = π / (2(π − arccos Δ)) for -1 < Δ ≤ 1.",
)
@register(
    XXZ1D,
    LuttingerVelocity,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/test_XXZ1D.jl",
    references=["Giamarchi2003", "desCloizeauxPearson1962"],
    notes="u = (πJ/2) sin γ / γ, γ = arccos Δ; -1 < Δ ≤ 1.",
)

@register(
    XXZ1D,
    NMRRelaxationExponent,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/XXZ/test_XXZ1D.jl",
    references=["ChitraGiamarchi1997", "Giamarchi2003"],
    notes="Leading θ_NMR = 1/(2K) - 1 (dominant transverse staggered channel, Δ_op = 1/(4K)) in the critical Luttinger liquid regime -1 < Δ ≤ 1; ChitraGiamarchi1997 Eq.27. Subdominant longitudinal channel gives T^{2K-1}.",
)
@register(
    XXZ1D,
    MassGap,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl",
    notes="E₁ - E₀ from full-spectrum dense ED.",
)
@register(
    XXZ1D,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl",
    references=["Giamarchi2003"],
    notes="0 in the critical regime -1 < Δ ≤ 1; gapped regime returns NaN with a warning.",
)

# ── Finite-T thermodynamic scalars (per-site at OBC) ──────────────────
@register(
    XXZ1D,
    FreeEnergy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_thermal.jl"
)
@register(
    XXZ1D,
    ThermalEntropy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_thermal.jl"
)
@register(
    XXZ1D,
    SpecificHeat,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_thermal.jl"
)

# ── Magnetisations (Pauli convention) ─────────────────────────────────
@register(
    XXZ1D,
    MagnetizationX,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)
@register(
    XXZ1D,
    MagnetizationY,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)
@register(
    XXZ1D,
    MagnetizationZ,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)

# ── Site-resolved local observables ───────────────────────────────────
@register(
    XXZ1D,
    LocalMagnetization{:x},
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)
@register(
    XXZ1D,
    LocalMagnetization{:y},
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)
@register(
    XXZ1D,
    LocalMagnetization{:z},
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)
@register(
    XXZ1D,
    EnergyLocal,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl",
    notes="Bonds split symmetrically: Σᵢ ε_i = ⟨H⟩.",
)

# ── Susceptibilities (β · Var(M_α) / N) ───────────────────────────────
@register(
    XXZ1D,
    SusceptibilityXX,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)
@register(
    XXZ1D,
    SusceptibilityYY,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)
@register(
    XXZ1D,
    SusceptibilityZZ,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl"
)

# ── Two-point correlators (static + connected) ────────────────────────
for CorrT in (
    SpinCorrelation{:x,:x},
    SpinCorrelation{:y,:y},
    SpinCorrelation{:z,:z},
    ConnectedSpinCorrelation{:x,:x},
    ConnectedSpinCorrelation{:y,:y},
    ConnectedSpinCorrelation{:z,:z},
)
    register!(
        XXZ1D,
        CorrT,
        OBC;
        method=:dense_ed,
        reliability=:high,
        tested_in="test/models/test_XXZ1D_observables.jl",
        notes="(i,j) ⟨σᵅᵢ σᵅⱼ⟩_β; the connected variant subtracts the disconnected piece.",
    )
end

# ── Entanglement (β = Inf default → ground-state pure-state entropy) ──
@register(
    XXZ1D,
    VonNeumannEntropy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl",
    notes="Pass subsystem length ℓ; β=Inf gives ground-state EE.",
)
@register(
    XXZ1D,
    RenyiEntropy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/test_XXZ1D_observables.jl",
    notes="S_α = log Tr ρ_A^α / (1 - α); pass subsystem ℓ and order α.",
)

# ── Free-fermion thermal at Infinite (Δ = 0 only; general Δ pending #108) ──
@register(
    XXZ1D,
    FreeEnergy,
    Infinite,
    method=:klumper_nlie,
    reliability=:medium,
    tested_in="test/models/quantum/XXZ/test_xxz_klumper_nlie.jl",
    references=["Mahan2000", "Coleman2015", "Takahashi1999", "Klumper1993"],
    notes="Δ = 0: XX free-fermion f(β) by QuadGK (exact); -1 < Δ < 1, |Δ| < 0.99: Klümper QTM NLIE (issue #521); |Δ| ≥ 0.99 or gapped: NaN + warn.",
)
@register(
    XXZ1D,
    ThermalEntropy,
    Infinite,
    method=:free_fermion_quadgk_or_klumper_nlie,
    reliability=:high,
    tested_in="test/standalone/test_xxz_xx_infinite.jl",
    references=["Mahan2000", "Coleman2015", "Klumper1993"],
    notes="XX free-fermion s(β); -1<Δ<1 (Δ≠0) routes through Klümper NLIE finite-diff (issue #521); NaN+warn at |Δ|≥0.99.",
)
@register(
    XXZ1D,
    SpecificHeat,
    Infinite,
    method=:free_fermion_quadgk_or_klumper_nlie,
    reliability=:high,
    tested_in="test/standalone/test_xxz_xx_infinite.jl",
    references=["Mahan2000", "Klumper1993"],
    notes="XX free-fermion closed form; -1<Δ<1 (Δ≠0) via Klümper NLIE finite-diff (issue #521); NaN+warn at |Δ|≥0.99.",
)
# ── Quench observables (Δ = 0 / XX free fermion only; issue #148 phase 1) ──
@register(
    XXZ1D,
    LoschmidtRateFunction,
    Infinite,
    method=:free_fermion_analytic,
    reliability=:high,
    tested_in="test/standalone/test_xxz_xx_quench.jl",
    references=["CalabreseEsslerFagotti2012", "Heyl2013", "EsslerFagotti2016"],
    notes="XX → XX quench Loschmidt rate λ(t) at Δ = 0 only; same-sign J ⇒ λ ≡ 0 (Fermi sea preserved), sign-flip ⇒ 0 (|GS(J₀)⟩ is a number eigenstate of H_f; Anderson orthogonality does not apply to the Loschmidt amplitude).  Δ ≠ 0 throws DomainError.",
)

# ── Ground-state energy density (Bethe-ansatz / ferromagnetic limit) ──
@register(
    XXZ1D,
    GroundStateEnergyDensity,
    Infinite,
    method=:bethe_ansatz,
    reliability=:high,
    tested_in="test/identities/test_identities_XXZ1D.jl",
    references=["YangYang1969", "desCloizeauxPearson1962"],
    notes="e₀(Δ) at the isotropic AF point (Δ=1) reduces to Heisenberg1D Hulthén value J(1/4-ln 2); at the FM point (Δ=-1) the aligned state is exact, e₀ = -J/4.",
)

# ── CC entanglement at Infinite via Universality(:XY) / (:Heisenberg) (#580 Phase 2)
@register(
    XXZ1D,
    VonNeumannEntropy,
    Infinite,
    method=:delegation,
    reliability=:high,
    tested_in="test/models/quantum/XXZ/test_xxz1d_cft_entanglement.jl",
    references=["CalabreseCardy2004"],
    notes="Critical regime -1 < Δ < 1: delegate to Universality(:XY) (c=1). Δ=1: route via Universality(:Heisenberg). |Δ|>1 gapped: DomainError.",
)

@register(
    XXZ1D,
    RenyiEntropy,
    Infinite,
    method=:delegation,
    reliability=:high,
    tested_in="test/models/quantum/XXZ/test_xxz1d_cft_entanglement.jl",
    references=["CalabreseCardy2004"],
    notes="Same critical-regime guard as VN. Standard c -> c*(1+1/alpha)/2 substitution. Reduces to VN at alpha=1.",
)

@register(
    XXZ1D,
    DynamicalSpinStructureFactor{:z,:z},
    Infinite,
    method=:exact_2spinon,
    reliability=:high,
    tested_in="test/models/quantum/XXZ/test_xxz_spinon.jl",
    references=["PerezCastillo2020"],
    notes="Exact longitudinal two-spinon dynamical structure factor S^{zz}(q, ω) in the massive regime Δ > 1.",
)
