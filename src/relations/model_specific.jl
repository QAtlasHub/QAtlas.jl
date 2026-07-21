# relations/model_specific.jl — relations that are TRUE but not UNIVERSAL (#730).
#
# AbstractQAtlas holds only what is model-independent; it purged these in ABQ#54
# ("keep only universal relations — move model-specific ones to QAtlas") because
# each of them assumes a particular model or band structure.  Their home is the
# implementing atlas, for the same reason the reference VALUES live here.
#
# They are declared with AbstractQAtlas's own `@relation` / `@inequality`, so
# they join the SAME registry and the same verbs (`residual` / `check` / `solve`
# / `relations_constraining`) as the universal ones — a caller cannot tell the
# two apart at the API, only by the domain tag.  Slots are type-keyed wherever
# the quantity genuinely exists in the shared vocabulary, which is what makes
# `relations_constraining(EdwardsAndersonParameter)` find these; slots that have
# no honest tag (a per-BOND energy, a relaxation time, the elementary charge)
# stay symbol-only rather than being forced onto an ill-fitting type.
#
# ── REGISTRY VISIBILITY (do not remove `__init__`) ─────────────────────────
# `@relation` performs its registry insertion as a LOAD-TIME SIDE EFFECT on
# AbstractQAtlas's `_RELATION_REGISTRY`.  Mutating another module's state does
# not survive QAtlas's precompilation — it is not serialized into our `.ji`, and
# a precompiled module body is never re-executed — so without re-registration
# these relations would exist as types with working `residual`/`solve` while
# being invisible to `all_relations()` / `relations_constraining()` in every
# fresh session.  This is not a theoretical worry: with `__init__` neutered and
# the package precompiled, a fresh session sees 123 relations and ALL SIX of
# these are missing; with it, 129 and all six present.  AbstractQAtlas documents
# exactly this for downstream packages ("re-register from your module `__init__`").  `MODEL_SPECIFIC_RELATIONS` below
# is the single list that both this file and `__init__` read, so a relation
# cannot be added here and silently forgotten there.

# ─── Spin glass (Edwards–Anderson / Nishimori / de Almeida–Thouless) ────────

"""
    EdwardsAndersonOrderParameter <: AbstractRelation

The Edwards–Anderson order parameter as the replica self-overlap
([EdwardsAnderson1975](@cite)),

`q_EA = [⟨s_i⟩²]`

(the disorder average of the squared thermal magnetization,
`(1/N) Σ_i [⟨s_i⟩²]`, the caller-supplied `overlap`).

Variables: `q_EA` ([`EdwardsAndersonParameter`](@ref)), `overlap`.
"""
@relation :spinglass EdwardsAndersonOrderParameter(
    q_EA::EdwardsAndersonParameter, overlap
) = q_EA - overlap

"""
    NishimoriEnergy <: AbstractRelation

The exact internal energy per BOND of the ±J Ising model on the Nishimori line
([Nishimori1981](@cite)),

`U = −J tanh(βJ)`,

a gauge-symmetry consequence that holds for any lattice and dimension — a rare
closed form in a disordered system, and a hard check on a disorder-averaged
simulation sitting on the Nishimori line.

`U` is per-bond, which is not one of the `Energy{G}` granularities, so it stays
symbol-only rather than being mis-tagged as a per-site energy.

Variables: `U`, `J`, `β` (or `T`).
"""
@relation :spinglass NishimoriEnergy(U, J, β) = U + J * tanh(β * J)

"""
    NishimoriMagnetizationOverlap <: AbstractRelation

The Nishimori gauge identity: on the Nishimori line the spin-glass order
parameter equals the ferromagnetic magnetization ([Nishimori1981](@cite)),

`q = m`,

so spin-glass and ferromagnetic order coincide — there is no
replica-symmetry-broken spin-glass phase below the Nishimori line.

Variables: `q` ([`EdwardsAndersonParameter`](@ref)), `m`.
"""
@relation :spinglass NishimoriMagnetizationOverlap(q::EdwardsAndersonParameter, m) = q - m

"""
    AlmeidaThoulessStability <: AbstractInequality

The de Almeida–Thouless replica-symmetric stability criterion
([AlmeidaThouless1978](@cite)): the replicon eigenvalue is non-negative in the
replica-symmetric phase,

`1 − (βJ)² [sech⁴(βh)] ≥ 0`

(the slack IS the replicon eigenvalue: `= 0` on the AT line, `< 0` in the
replica-symmetry-broken phase).  `sech4_avg = [sech⁴(βh)]` is the
disorder-averaged local-field factor.

Variables: `βJ`, `sech4_avg`.
"""
@inequality :spinglass AlmeidaThoulessStability(βJ, sech4_avg) = 1 - βJ^2 * sech4_avg

# ─── Transport (model-specific band / scattering assumptions) ───────────────

"""
    DrudeMobility <: AbstractRelation

Mobility in the Drude relaxation-time model ([Drude1900](@cite)),

`μ = e τ / m*`,

which assumes a single relaxation time τ and a parabolic band — unlike
`MobilityConductivity` (`σ = n e μ`), which is the *definition* of μ and stays
universal in AbstractQAtlas.

Variables: `μ` ([`Mobility`](@ref)), `e`, `τ`, `m` ([`EffectiveMass`](@ref)).
"""
@relation :transport DrudeMobility(μ::Mobility, e, τ, m::EffectiveMass) = μ - e * τ / m

"""
    SingleBandHall <: AbstractRelation

The single-band Hall coefficient ([Hall1879](@cite)),

`R_H = 1 / (n e)`,

valid only when one carrier type dominates — a two-band or compensated metal
violates it, which is precisely why it is not universal.

Variables: `R_H` ([`HallCoefficient`](@ref)), `n` ([`CarrierDensity`](@ref)), `e`.
"""
@relation :transport SingleBandHall(R_H::HallCoefficient, n::CarrierDensity, e) =
    R_H * n * e - 1

"""
    MODEL_SPECIFIC_RELATIONS

Every relation declared in this file, in one list.  `QAtlas.__init__` re-inserts
these into AbstractQAtlas's registry on each load — see the file header for why
the `@relation` insertion alone does not survive precompilation.  Adding a
relation above without adding it here makes it invisible to `all_relations()`;
`test/relations/test_model_specific.jl` fails on exactly that omission.
"""
const MODEL_SPECIFIC_RELATIONS = (
    EdwardsAndersonOrderParameter(),
    NishimoriEnergy(),
    NishimoriMagnetizationOverlap(),
    AlmeidaThoulessStability(),
    DrudeMobility(),
    SingleBandHall(),
)
