# identity_registry.jl — quantity ↔ quantity identities (@identity_edge, #698).
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
# finite_N=6 matches the hand-written identity files (Heisenberg1D at OBC(6),
# S1Heisenberg1D at OBC(4)): the Gibbs relation is an internal-consistency
# statement of the SAME finite system, exact at every N, so small N loses no
# coverage — while N=8 put the spin-1 hubs at 3^8 = 6561-dimensional dense ED
# per fetch and made the generated suite CI's long pole (52 min → seconds).
# The ARITHMETIC is not restated here: it is AbstractQAtlas's `FreeEnergyLegendre`
# (`F - (U - S/β)`, relations/fundamental.jl), and this edge asks it to derive `f`
# from the hub's own `e` and `s`.  That makes the generated check the conformance
# statement the split was for (#734): QAtlas's fetched FreeEnergy must equal what
# the universal relation implies from QAtlas's fetched Energy and ThermalEntropy.
# `solve` is generic over any variable the relation is affine in, so there is no
# QAtlas-side rearrangement to drift out of sync.  Both sides stay physical
# numbers (fetched f vs derived f), keeping failures diagnosable — a residual-vs-0
# pair would not.  Granularity: the relation is type-keyed per-site, matching the
# `Energy{:per_site}` participant below.
@identity_edge(
    :gibbs,
    quantities = (e=Energy{:per_site}, f=FreeEnergy, s=ThermalEntropy),
    check = (v, p) -> (v.f, solve(FreeEnergyLegendre(), Val(:F); U=v.e, S=v.s, β=p.beta)),
    sweep = (beta=[0.5, 1.0, 2.0],),
    finite_N = 6,
    exclusions = [
        SSH => "Energy fetch returns T=0 ground-state energy (beta swallowed); thermal ε not implemented — Gibbs does not apply as stated (#508 kwargs-swallow audit)",
        TightBinding1D => "Energy fetch returns T=0 ground-state energy (beta swallowed); thermal ε not implemented — Gibbs does not apply as stated (#508 kwargs-swallow audit)",
        TightBindingV1D => "Energy fetch returns T=0 ground-state energy (beta swallowed); thermal ε not implemented — Gibbs does not apply as stated (#508 kwargs-swallow audit)",
        AKLT1D => "finite-β canonical thermodynamics supports β = ∞ only (HTSE is a separate :approx scheme, #506)",
        # Same shape as the three above: this hub's `Energy{:per_site}` is a T = 0 closed
        # form (E_0/N = -pi^2 J / 24) and swallows `beta`. It became visible to this
        # identity only when `GroundStateEnergyDensity` folded into `Energy{:per_site}` --
        # before that the hub had no Energy row at all, so no gibbs edge was generated.
        # The finite-T rows it does have (FreeEnergy, ThermalEntropy, SpecificHeat) are
        # the c = 1 CFT low-T stopgap; there is no thermal <H> to close the identity with.
        HaldaneShastry => "Energy fetch returns the T=0 ground-state energy (beta swallowed); the finite-T rows are the c=1 CFT stopgap and carry no thermal ε — Gibbs does not apply as stated",
        # Heisenberg1D, same shape, but keyed per-bc: its `Energy{:per_site}` at Infinite
        # is the Hulthén T = 0 closed form `e₀ = J(1/4 - ln 2)` and takes no `beta` kwarg
        # at all (so the identity errors rather than mismatching). At OBC the hub's Energy
        # row is `Energy{:total}`, which IS the beta-dependent thermal energy delegated to
        # XXZ1D(Δ=1) — a bare `Heisenberg1D` key would cover every bc and hide that edge.
        (Heisenberg1D, Infinite) => "Infinite-bc Energy fetch is the Hulthén T=0 ground-state energy density and takes no beta; the finite-T rows at this bc carry no thermal ε — Gibbs does not apply as stated (#508 kwargs-swallow audit)",
        IsingTriangular => "default J > 0 is the frustrated AFM branch (no Houtappel closed form); finite-T requires J < 0",
        (IsingSquare, PBC) => "2D PBC fetches take Lx/Ly kwargs, not bc.N — generator finite-size materialization is 1D-only",
        (KitaevHoneycomb, OBC) => "2D OBC fetches take Lx/Ly kwargs, not bc.N — generator finite-size materialization is 1D-only",
    ],
    notes = "ε(β) = f(β) + T·s(β); per-site convention throughout (Energy granularity routing applies).",
)

# ── SU(2) isotropy of the susceptibility family ───────────────────────
# χ_xx = χ_yy = χ_zz for every model whose @symmetry profile declares
# internal=:SU2 — pairwise component equality over AbstractSusceptibility.
@identity_edge(
    :su2_susceptibility_isotropy,
    family = AbstractSusceptibility,
    requires_internal = :SU2,
    sweep = (beta=[0.5, 1.0],),
    finite_N = 6,
    notes = "SU(2) invariance forces equal diagonal susceptibilities along all spin axes.",
)

# ── Σᵢ εᵢ = ⟨H⟩ ───────────────────────────────────────────────────────
# The local energy density is DEFINED by this sum rule — EnergyLocal's own
# docstring says "defined so that Σᵢ ε_i = ⟨H⟩_β" — so the check is the
# definition, asked of the implementation.
#
# It is a genuine cross-route check, which is the part worth stating.  On the
# free-fermion hubs the two sides are computed from DIFFERENT objects: ε_i sums
# local terms read off the thermal Majorana covariance matrix Σ, while
# Energy{:total} sums `-Λₙ/2 · tanh(βΛₙ/2)` over the BdG spectrum.  Neither is
# derived from the other.  MEASURED: the residual sits at rounding scale
# (1e-16 … 1e-14) and CHANGES SIGN with N, which is what an independent pair
# looks like; a shared routine would give exactly 0 at every N, the way
# LoschmidtRate does (see its MATERIALIZABLE_BUT_UNWIRED entry).
#
# NOT EQUALLY STRONG ON EVERY HUB, and the registry says which.  On the
# free-fermion hubs (TFIM, XYh1D — cost=:polynomial) the two sides come from
# different objects, as above.  On the dense-ED hubs (Heisenberg1D, XXZ1D,
# S1Heisenberg1D — cost=:exponential) BOTH sides come from `eigen(H)` of the same
# matrix: `Energy{:total}` is `Σₙ wₙ Eₙ`, `EnergyLocal` is `Σᵢ Σₙ wₙ ⟨n|bᵢ|n⟩`.
# There the check tests the LOCAL DECOMPOSITION — that the bond operators and the
# ½ splitting really sum to H — and nothing about the Hamiltonian or the spectrum,
# which are shared.  Still a real check, and it would catch a wrong bond operator;
# just a narrower one than the free-fermion hubs get.  Read `cost` alongside the
# pass count.
#
# And it DISCRIMINATES, which agreement alone does not establish.  Control,
# measured at N = 6: comparing Σᵢεᵢ(β=1) against E_total at the WRONG temperature
# (β=2) departs by 0.45 … 1.4 on the six hubs, against a matched residual of
# ≤ 3e-15 — a margin of 1e14 … 1e15.  So the agreement is a fact about the
# implementation, not an identity that holds whatever the numbers are.
#
# This is also the first edge whose participant is ARRAY-valued.  Nothing in
# `identity!` needed changing: the check lambda reduces the profile itself, so
# no "reduction relation" kind was required upstream — the narrower conclusion
# than the one #819 originally reached for this group.
@identity_edge(
    :local_energy_sum_rule,
    quantities = (ε=EnergyLocal, E=Energy{:total}),
    check = (v, p) -> (sum(v.ε), v.E),
    sweep = (beta=[0.5, 1.0, 2.0],),
    finite_N = 6,
    notes = "Σᵢ εᵢ = ⟨H⟩_β — the defining property of the local energy density, checked against a total energy computed from the spectrum rather than from the profile.",
)

# ── A quench relaxes to its GGE ───────────────────────────────────────
# lim_{t→∞} ⟨σˣ⟩(t)  =  ⟨σˣ⟩_GGE  for the infinite TFIM after h₀ → h_f.
#
# WHY THIS AND NOT A PROTOCOL AXIS.  #819 step (2) proposed that the quench
# quantities "need a protocol axis that does not exist yet in either package".
# The law that would justify such an axis turns out to be statable with the
# vocabulary already present: both quantities take `initial` and the quench one
# takes `t`, and `_sweep_points` carries a MODEL-valued sweep entry without
# changing anything.  That is the fourth §4(a) upstream change to shrink to
# nothing under measurement (after `SiteSupport`, the group law, and the
# reduction relation) — so the axis, if it is ever wanted, still needs its own
# first user.
#
# WHY IT IS INDEPENDENT, and exactly how far.  The two integrands are different
# expressions, not one routine called twice:
#
#     quench :  (1/π)∫₀^π [cos2θ_f·cos2Δθ + sin2θ_f·sin2Δθ·cos(2Λ_f t)]
#     GGE    :  (2/π)∫₀^π (h_f − J cos k)/Λ_f · (1 − 2 n_k),  n_k = sin²Δθ
#
# and they agree only via the non-trivial identities cos2θ_f = (h_f − J cos k)/Λ_f
# and cos2Δθ = 1 − 2n_k.  They also use SEPARATE helper families —
# `_tfim_two_theta`/`_tfim_lambda` versus `_tfim_gge_two_theta`/`_tfim_dispersion`.
# MEASURED: those pairs agree to every printed digit, i.e. they are DUPLICATE
# implementations of the same convention.  So this edge cannot detect an error in
# the Bogoliubov-angle convention itself — both sides would move together — and
# it is stated here rather than left for a reader to discover.
#
# MEASURED, at the hub's own default h_f = 1.0 (the CRITICAL point, and the
# gapless case where dephasing is slowest — so this is the unfavourable regime,
# not a flattering one).  |residual| at t = 800 against the wrong-h₀ margin:
#
#     h₀      0.2      0.5      0.8      1.5      2.0      3.0
#     resid   2.2e-7   1.1e-7   3.6e-8   6.5e-8   1.1e-7   1.6e-7
#     margin  1.9e-2   1.3e-2   5.8e-3   4.5e-3   3.4e-3   1.9e-3
#
# — a separation of ~10⁴, so `atol = 1e-5` sits 45× above the worst residual and
# 190× below the smallest margin.
#
# WHY t = 800 AND NOT LARGER.  The residual does NOT keep falling: it plateaus at
# ~1e-7 from t ≈ 800 and is unchanged at t = 3200, because the `cos(2Λ_f t)` term
# makes the integrand oscillate faster than the quadrature's `rtol = 1e-12` can
# follow.  Raising `t` buys nothing and eventually costs accuracy, which is the
# opposite of what "increase t for the long-time limit" would suggest.
@identity_edge(
    :quench_relaxes_to_gge,
    quantities = (m=QuenchLocalMagnetization{:x}, g=GGEValue{MagnetizationX}),
    check = (v, p) -> (v.m, v.g),
    sweep = (
        initial=[TFIM(1.0, 0.2), TFIM(1.0, 0.5), TFIM(1.0, 2.0), TFIM(1.0, 3.0)], t=[800.0]
    ),
    rtol = 1e-5,
    atol = 1e-5,
    notes = "lim_{t→∞} ⟨σˣ⟩(t) = ⟨σˣ⟩_GGE after a sudden h₀ → h_f quench of the infinite TFIM. Evaluated at t = 800, where the residual has plateaued; see the header for why larger t does not help.",
    references = ["CalabreseEsslerFagotti2012"],
)

# ── SU(2) isotropy of the magnetization family ────────────────────────
# m_x = m_y = m_z (all zero in a finite-N canonical ensemble, but the edge
# asserts the symmetry statement — equality — not the value).
@identity_edge(
    :su2_magnetization_isotropy,
    family = AbstractMagnetization,
    requires_internal = :SU2,
    sweep = (beta=[0.5, 1.0],),
    finite_N = 6,
    notes = "SU(2) invariance forces equal (vanishing) magnetization along all spin axes.",
)

# DELETION CRITERION for the hand-written duplicates (#698 migration): the
# harness identities in test/util/thermodynamic_identities.jl and the
# per-model files in test/identities/ stay until (a) the derivative-form
# identities (c_v, m_x = -∂f/∂h) are generatable, (b) the zero-value SU(2)
# checks (m_α = 0, not just pairwise equality) have a generated form, and
# (c) the parameter-conditional gates (XXZ1D at Δ ≈ 1) are expressible via
# profile `at` predicates.  Until then the overlap (Gibbs + pairwise
# isotropy on Heisenberg1D/S1Heisenberg1D) is accepted as the parity
# reference, NOT an oversight.
