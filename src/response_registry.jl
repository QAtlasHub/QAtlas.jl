# response_registry.jl — declared RESPONSE edges (core/response.jl).
#
# These are the AbstractQAtlas relations that were reachable on QAtlas hubs but
# unusable, because one slot is a derivative rather than a fetchable value.  With
# core/derivative.jl supplying it they become ordinary generated checks.
#
# Both are exact thermodynamic identities that hold at every N and β, so a small
# finite_N loses no coverage — the same argument :gibbs makes.

# ── S = -∂F/∂T ────────────────────────────────────────────────────────
@response(
    :entropy_response,
    relation = EntropyResponse,
    derived = (dF_dT=∂(FreeEnergy, :T),),
    sweep = (beta=[0.5, 1.0, 2.0],),
    finite_N = 6,
    exclusions = _THERMO_DERIVATIVE_EXCLUSIONS,
    notes = "S = -∂F/∂T — the Maxwell relation between entropy and the free energy.",
)

# ── C = T ∂S/∂T ───────────────────────────────────────────────────────
@response(
    :specific_heat_from_entropy,
    relation = SpecificHeatFromEntropy,
    derived = (dS_dT=∂(ThermalEntropy, :T),),
    sweep = (beta=[0.5, 1.0, 2.0],),
    finite_N = 6,
    exclusions = _THERMO_DERIVATIVE_EXCLUSIONS,
    # Onsager's C and S are quadrature evaluations: measured, they reproduce
    # C = T ∂S/∂T to 2.5e-5 on IsingSquare with AD and with finite differences
    # alike.  That is the accuracy of the FETCHES, so AD's 1e-6 would report a
    # quadrature residue as a physics failure.
    rtol_floor = 1e-4,
    notes = "C = T ∂S/∂T — the caloric definition, independent of the C = β² Var(E) route.",
)

# Hubs whose `Energy{:per_site}` does not respond to the swept β.  Split out from
# the list above because the reason is different and the affected edges are the
# ones with an ENERGY subject or an energy derivative — a different hub set from
# the entropy-subject edges, which is why they surfaced only when
# :gibbs_helmholtz arrived.  The first three carry the wording
# identity_registry.jl already uses for :gibbs.
const _BETA_PINNED_ENERGY_MODELS = [
    SSH => "Energy fetch returns T=0 ground-state energy (beta swallowed); thermal ε not implemented — the relation does not apply as stated (#508 kwargs-swallow audit)",
    TightBinding1D => "Energy fetch returns T=0 ground-state energy (beta swallowed); thermal ε not implemented (#508 kwargs-swallow audit)",
    TightBindingV1D => "Energy fetch returns T=0 ground-state energy (beta swallowed); thermal ε not implemented (#508 kwargs-swallow audit)",
    # CIRCULAR, not merely β-pinned: SixVertex's Energy{:per_site} IS ∂(βF)/∂β,
    # computed internally by a hard-coded central difference of f(a^β, b^β, c^β)
    # at β = 1.  Checking Gibbs–Helmholtz against it would verify the model with
    # the very relation it already assumes — the circularity the atlas's
    # independence axis exists to prevent — and it ignores the swept β besides.
    SixVertex => "Energy{:per_site} is itself computed as ∂(βF)/∂β by an internal finite difference pinned at β = 1; checking this relation against it would be circular",
]

# ── U = ∂(βF)/∂β ──────────────────────────────────────────────────────
# The Gibbs–Helmholtz relation.  What is differentiated is the PRODUCT βF, not
# F alone — that is what `of` is for.  Independent of :entropy_response despite
# relating the same two potentials: this is the β-derivative of βF, that one is
# the T-derivative of F, and an implementation can satisfy one while breaking
# the other.
@response(
    :gibbs_helmholtz,
    relation = GibbsHelmholtz,
    derived = (dβF_dβ=∂(FreeEnergy, :β; of=(F, β) -> β * F),),
    sweep = (beta=[0.5, 1.0, 2.0],),
    finite_N = 6,
    exclusions = vcat(_THERMO_DERIVATIVE_EXCLUSIONS, _BETA_PINNED_ENERGY_MODELS),
    rtol_floor = 1e-4,
    notes = "U = ∂(βF)/∂β — Gibbs–Helmholtz; equivalently U = -∂ln Z/∂β.",
)

# ── C = β² Var(E) ─────────────────────────────────────────────────────
# The energy-fluctuation route to the specific heat.  `var_E` is not a fetchable
# quantity, but it does not need to be: Var(E) = -∂⟨E⟩/∂β exactly, so `then`
# supplies it by negating the derivative of the energy.  With `Energy{:per_site}`
# the relation's `N` stays 1 — per-site C against per-site variance.
#
# This is a genuinely INDEPENDENT route to C: :specific_heat_from_entropy gets it
# from the entropy, this one from the energy.  A model that computes C by one
# formula and S or U by another is exactly what these two disagree on.
@response(
    :specific_heat_fdt,
    relation = SpecificHeatFDT,
    derived = (var_E=∂(Energy{:per_site}, :β; then=d -> -d),),
    sweep = (beta=[0.5, 1.0, 2.0],),
    finite_N = 6,
    exclusions = vcat(_THERMO_DERIVATIVE_EXCLUSIONS, _BETA_PINNED_ENERGY_MODELS),
    rtol_floor = 1e-4,
    notes = "C = β² Var(E), Var(E) = -∂⟨E⟩/∂β — the fluctuation route, independent of C = T ∂S/∂T.",
)

# ── M = −∂F/∂h  and  χ = ∂M/∂h ────────────────────────────────────────
# MODEL-axis derivatives: `h` is a field of the model, not a fetch kwarg, so
# these rebuild the model at each step and are finite-difference only (see
# `_diff_target`).  Reported at the finite-difference tolerance accordingly.
#
# WHY AN ALLOW-LIST AND NOT AN EXCLUSION LIST.  AbstractQAtlas types the subject
# `Magnetization{:z}` / `Susceptibility{(:z,:z)}`, so these relations hold only
# where `h` is the LONGITUDINAL field.  It is not, in this atlas:
#
#     -h Σ σᶻ  (longitudinal, valid)   IsingChain1D, CurieWeissIsing, LongRangeXY1D
#     -h Σ σˣ  (TRANSVERSE, invalid)   TFIM, LongRangeIsing1D
#
# For a transverse-field model −∂F/∂h is ⟨σˣ⟩, not M_z. Both differentiation
# backends would agree on that wrong number, so `derivative_agreement` cannot
# save us — the cross-check catches a bad METHOD, never a relation applied to a
# Hamiltonian it does not describe. Hence opt-in: a model added later is skipped
# rather than checked against physics that does not apply to it.
const _LONGITUDINAL_FIELD_MODELS = [IsingChain1D, CurieWeissIsing]

@response(
    :magnetization_response,
    relation = MagnetizationResponse,
    derived = (dF_dh=∂(FreeEnergy, :h; then=d -> -d),),
    models = _LONGITUDINAL_FIELD_MODELS,
    sweep = (beta=[0.5, 1.0, 2.0],),
    finite_N = 6,
    exclusions = _THERMO_DERIVATIVE_EXCLUSIONS,
    rtol_floor = 1e-4,
    notes = "M_z = -∂F/∂h — valid only where h is the longitudinal field.",
)

@response(
    :susceptibility_response,
    relation = SusceptibilityResponse,
    derived = (dM_dh=∂(Magnetization{:z}, :h),),
    models = _LONGITUDINAL_FIELD_MODELS,
    sweep = (beta=[0.5, 1.0, 2.0],),
    finite_N = 6,
    exclusions = _THERMO_DERIVATIVE_EXCLUSIONS,
    rtol_floor = 1e-4,
    notes = "χ_zz = ∂M_z/∂h — the isothermal susceptibility as a field response.",
)
