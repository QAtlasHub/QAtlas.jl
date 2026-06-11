# identity_registry.jl — quantity ↔ quantity identities (@identity, #698).
#
# Migration discipline (per #698): the Gibbs relation is ported END-TO-END
# first; the hand-written harness (test/util/thermodynamic_identities.jl,
# exercised by test/identities/) stays in place as the reference
# implementation until full parity is confirmed.  The derivative-form
# identities (c_v = -β² ∂ε/∂β, m_x = -∂f/∂h) need a derivative engine
# (ForwardDiff / central differences) that src/ deliberately does not depend
# on — they remain harness-only until the generator grows a pluggable
# derivative backend.
#
# The isotropy identities are the #690 × #700 integration: declared ONCE
# against the family supertype, gated on the @symmetry profiles — they
# replace the harness's hand-coded `is_su2_symmetric` SU(2) axis-equality
# rules with registry edges, and extend automatically to every model that
# later declares `internal=:SU2`.

# ── Gibbs relation: ε = f + T·s ───────────────────────────────────────
# Pure algebra over three independent fetches — no derivative engine needed.
# Catches per-site/total drift, entropy sign errors, missing T factors on
# every (model, bc) hub that implements the thermal triple, current and
# future, with zero hand-written tests.
# Exclusions (visible :skip, never silent) — what the FIRST generated run of
# this edge surfaced across the 24 implementing hubs:
#   * SSH / TightBinding1D / TightBindingV1D: the Energy fetch SWALLOWS the
#     beta kwarg and returns the T=0 ground-state energy (−2/π for the
#     uniform chain) while f and s are genuinely finite-T (f + T·s converges
#     to that value only as β → ∞) — exactly the convention drift this layer
#     exists to catch.  Skip until those models implement thermal ε; tracked
#     with the kwargs-swallow audit (#508).
#   * AKLT1D: finite-β canonical path supports β = ∞ only (HTSE is a separate
#     :approx scheme, #506).
#   * IsingTriangular: the default instance (J > 0) is the frustrated AFM
#     branch with no Houtappel closed form (finite-T needs J < 0).
#   * IsingSquare@PBC / KitaevHoneycomb@OBC: 2D fetches take Lx/Ly kwargs,
#     not bc.N — the generator's 1D finite_N materialization does not apply
#     (their Infinite hubs DO run and pass).
@identity(
    :gibbs,
    quantities = (e=Energy{:per_site}, f=FreeEnergy, s=ThermalEntropy),
    check = (v, p) -> (v.e, v.f + v.s / p.beta),
    sweep = (beta=[0.5, 1.0, 2.0],),
    exclusions = [
        SSH => "Energy fetch returns T=0 ground-state energy (beta swallowed); thermal ε not implemented — Gibbs does not apply as stated (#508 kwargs-swallow audit)",
        TightBinding1D => "Energy fetch returns T=0 ground-state energy (beta swallowed); thermal ε not implemented — Gibbs does not apply as stated (#508 kwargs-swallow audit)",
        TightBindingV1D => "Energy fetch returns T=0 ground-state energy (beta swallowed); thermal ε not implemented — Gibbs does not apply as stated (#508 kwargs-swallow audit)",
        AKLT1D => "finite-β canonical thermodynamics supports β = ∞ only (HTSE is a separate :approx scheme, #506)",
        IsingTriangular => "default J > 0 is the frustrated AFM branch (no Houtappel closed form); finite-T requires J < 0",
        (IsingSquare, PBC) => "2D PBC fetches take Lx/Ly kwargs, not bc.N — generator finite-size materialization is 1D-only",
        (KitaevHoneycomb, OBC) => "2D OBC fetches take Lx/Ly kwargs, not bc.N — generator finite-size materialization is 1D-only",
    ],
    notes = "ε(β) = f(β) + T·s(β); per-site convention throughout (Energy granularity routing applies).",
)

# ── SU(2) isotropy of the susceptibility family ───────────────────────
# χ_xx = χ_yy = χ_zz for every model whose @symmetry profile declares
# internal=:SU2 — pairwise component equality over AbstractSusceptibility.
@identity(
    :su2_susceptibility_isotropy,
    family = AbstractSusceptibility,
    requires_internal = :SU2,
    sweep = (beta=[0.5, 1.0],),
    notes = "SU(2) invariance forces equal diagonal susceptibilities along all spin axes.",
)

# ── SU(2) isotropy of the magnetization family ────────────────────────
# m_x = m_y = m_z (all zero in a finite-N canonical ensemble, but the edge
# asserts the symmetry statement — equality — not the value).
@identity(
    :su2_magnetization_isotropy,
    family = AbstractMagnetization,
    requires_internal = :SU2,
    sweep = (beta=[0.5, 1.0],),
    notes = "SU(2) invariance forces equal (vanishing) magnetization along all spin axes.",
)
