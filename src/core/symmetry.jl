# core/symmetry.jl — model symmetry attributes (#700, over the #697 kernel).
#
# Unlike the other constraint stores these are node ATTRIBUTES, not edges: a
# declarative record of each model family's symmetry data (internal group,
# translation invariance, time reversal, on-site spin) plus its declared
# spectral facts (gapped? ground-state degeneracy?).  Three consumers:
#
#   * LSM-type COHERENCE CHECKS (`check_lsm_consistency`, run in the C-suite):
#     the first checks that encode *theorems* over the registry rather than
#     bookkeeping invariants.  Declarations that contradict
#     Lieb–Schultz–Mattis are surfaced as `:error` — either a registry
#     mistake or a genuinely interesting claim; both are worth a loud finding.
#
#   * SYMMETRY-GATED IDENTITY GENERATION (core/identity.jl): an identity edge
#     with `requires_internal=:SU2` generates its checks exactly for the
#     models whose profile declares that internal symmetry — replacing the
#     hand-maintained `is_su2_symmetric` trait of the test harness with a
#     registry query.
#
#   * SPECTRAL CORROBORATION (`symmetry_checks`, the :symmetry kind of
#     `generated_checks`): a declared `gapped` fact is cross-checked against
#     the registered `MassGap` implementation, so the profile store and
#     REGISTRY cannot silently contradict each other.
#
# Profiles describe the model FAMILY at generic parameters.  A family whose
# symmetry is parameter-dependent (XXZ1D is U(1) generically, SU(2) only at
# Δ=1) declares the generic group; parameter-conditional enhancements stay in
# the per-model harness until profiles grow `at` predicates (the @realizes
# pattern).  Same for `gapped`: `nothing` means parameter-dependent/unknown —
# no claim, no check.

"""
    SymmetryProfile

The declared symmetry attributes of one model family — see [`@symmetry`](@ref).
`internal` is the on-site/internal symmetry group tag (`:SU2`, `:U1`, `:Z2`,
`:Z2xZ2`, `:none`, …); `site_spin` the on-site spin (`1//2`, `1`, …, or
`nothing` for non-spin models); `gapped` / `gs_degeneracy` the declared bulk
spectral facts at generic parameters (`nothing` = parameter-dependent or not
declared).
"""
struct SymmetryProfile
    model::Type
    internal::Symbol
    translation::Bool
    time_reversal::Bool
    site_spin::Union{Rational{Int},Nothing}
    gapped::Union{Bool,Nothing}
    gs_degeneracy::Union{Int,Nothing}
    notes::String
    references::Vector{String}
end

"""
    SYMMETRY_PROFILES :: Vector{SymmetryProfile}

The model-symmetry attribute store, populated at include-time by
[`symmetry!`](@ref) / [`@symmetry`](@ref) (one profile per model family).
Query with [`symmetry_profile`](@ref) / [`models_with_symmetry`](@ref).
"""
const SYMMETRY_PROFILES = SymmetryProfile[]

"""
    symmetry!(model_T; internal, translation=false, time_reversal=false,
              site_spin=nothing, gapped=nothing, gs_degeneracy=nothing,
              notes="", references=String[])

Record `model_T`'s symmetry profile.  Invariants: one profile per model;
`gs_degeneracy` is only meaningful for a declared `gapped=true` family and
must be `≥ 1`.
"""
function symmetry!(
    model_T::Type;
    internal::Symbol,
    translation::Bool=false,
    time_reversal::Bool=false,
    site_spin::Union{Rational{Int},Integer,Nothing}=nothing,
    gapped::Union{Bool,Nothing}=nothing,
    gs_degeneracy::Union{Int,Nothing}=nothing,
    notes::AbstractString="",
    references::AbstractVector{<:AbstractString}=String[],
)
    any(p -> p.model === model_T, SYMMETRY_PROFILES) && throw(
        ArgumentError("symmetry!: $(model_T) already has a profile (one per model family)"),
    )
    if gs_degeneracy !== nothing
        gapped === true || throw(
            ArgumentError(
                "symmetry!: gs_degeneracy is the degeneracy OF the gapped ground " *
                "state — declare gapped=true (got gapped=$(repr(gapped)))",
            ),
        )
        gs_degeneracy ≥ 1 || throw(ArgumentError("symmetry!: gs_degeneracy must be ≥ 1"))
    end
    spin = site_spin === nothing ? nothing : Rational{Int}(site_spin)
    push!(
        SYMMETRY_PROFILES,
        SymmetryProfile(
            model_T,
            internal,
            translation,
            time_reversal,
            spin,
            gapped,
            gs_degeneracy,
            String(notes),
            String[r for r in references],
        ),
    )
    return nothing
end

"""
    @symmetry Model internal=:SU2 translation=true site_spin=1//2 …

Macro sugar around [`symmetry!`](@ref): the positional `Model` is spliced as a
type; the remaining `key=value` pairs are forwarded as keyword arguments.

```julia
@symmetry Heisenberg1D internal=:SU2 translation=true time_reversal=true site_spin=1//2 gapped=false
```
"""
macro symmetry(model_T, kwargs...)
    return _forward_kw_macro(symmetry!, :symmetry, (model_T,), kwargs)
end

"""
    symmetry_profile(model) -> Union{SymmetryProfile,Nothing}

The declared symmetry profile of `model` (instance or type), or `nothing` if
the model has no [`@symmetry`](@ref) declaration yet.
"""
function symmetry_profile(model)
    m_T = _as_type(model)
    for p in SYMMETRY_PROFILES
        p.model === m_T && return p
    end
    return nothing
end

"""
    models_with_symmetry(internal::Symbol) -> Vector{Type}

The model families whose profile declares the given `internal` symmetry group
— the registry query behind symmetry-gated identity generation.
"""
function models_with_symmetry(internal::Symbol)
    return Type[p.model for p in SYMMETRY_PROFILES if p.internal === internal]
end

# ──────────────────────────────────────────────────────────────────────
# C10 — LSM-type consistency (theorems as static coherence checks)
# ──────────────────────────────────────────────────────────────────────

# Internal symmetry groups that carry the Lieb–Schultz–Mattis obstruction for
# a translation-invariant chain of half-odd-integer on-site spin: the ground
# state can be neither gapped-and-unique nor symmetric — it is gapless or
# degenerate (LSM 1961; Affleck–Lieb 1986; Oshikawa 2000; Hastings 2004).
const _LSM_INTERNAL = (:SU2, :U1)

"""
    check_lsm_consistency() -> Vector{CoherenceFinding}

C10: Lieb–Schultz–Mattis coherence over the [`@symmetry`](@ref) store — no
test execution, declarations only.  For every profile where the LSM theorem
applies (internal ⊇ U(1) spin rotation, translation invariance, half-odd-
integer `site_spin`):

  * declared `gapped=true` with `gs_degeneracy == 1` is a `:error` — the
    declaration contradicts the theorem (registry mistake, or a claim that
    needs extraordinary evidence);
  * declared `gapped=true` with no `gs_degeneracy` is a `:gap` — the
    profile owes the degeneracy that reconciles it with LSM.
"""
function check_lsm_consistency()
    out = CoherenceFinding[]
    for p in SYMMETRY_PROFILES
        p.translation || continue
        p.internal in _LSM_INTERNAL || continue
        p.site_spin === nothing && continue
        isinteger(p.site_spin) && continue   # integer spin: Haldane, no obstruction
        p.gapped === true || continue
        if p.gs_degeneracy === nothing
            push!(
                out,
                CoherenceFinding(
                    :lsm_consistency,
                    :gap,
                    "$(_kgshort(p.model)) declares a gapped translation-invariant " *
                    "spin-$(p.site_spin) chain with internal :$(p.internal) but no " *
                    "gs_degeneracy — LSM requires degenerate ground states; declare " *
                    "gs_degeneracy (≥ 2) or revisit gapped=true",
                ),
            )
        elseif p.gs_degeneracy == 1
            push!(
                out,
                CoherenceFinding(
                    :lsm_consistency,
                    :error,
                    "$(_kgshort(p.model)) declares gapped=true with a UNIQUE ground " *
                    "state on a translation-invariant spin-$(p.site_spin) chain with " *
                    "internal :$(p.internal) — this contradicts Lieb–Schultz–Mattis " *
                    "(gapless or degenerate); the profile (or the model's spectral " *
                    "claim) is wrong",
                ),
            )
        end
    end
    return out
end

register_edge_store!(
    :symmetry, SYMMETRY_PROFILES; location_of=p -> "symmetry $(_kgshort(p.model))"
)

# ──────────────────────────────────────────────────────────────────────
# Generator — the :symmetry kind of generated_checks()
# ──────────────────────────────────────────────────────────────────────

# Threshold separating "vanishing" from "open" gaps in the corroboration
# checks: closed-form / literature gaps are exact, so the floor only absorbs
# round-off — a declared-gapless model must fetch |Δ| ≤ this, a
# declared-gapped one Δ > this.
const _GAP_ATOL = 1e-10

"""
    symmetry_checks() -> Vector{GeneratedCheck}

Cross-store corroboration of the [`@symmetry`](@ref) spectral declarations:
for every profile that declares `gapped` (true or false) AND whose model has
a canonical, independent `MassGap` row at `Infinite`, emit a check comparing
the declaration against the fetched gap.  This closes the
two-sources-of-truth seam between the profile store and `REGISTRY`: a
`MassGap` implementation change that contradicts the declared profile (or a
wrong profile) fails loudly instead of drifting silently.  Profiles with
`gapped=nothing` (parameter-dependent families) and models without an
independent `MassGap` row emit nothing — the declaration carries no claim to
corroborate, or no second implementation exists to corroborate it against.
"""
function symmetry_checks()
    out = GeneratedCheck[]
    for p in SYMMETRY_PROFILES
        p.gapped === nothing && continue
        row = _canonical_row(p.model, MassGap, Infinite)
        (row !== nothing && _is_independent_row(row)) || continue
        model_T = p.model
        gapped = p.gapped
        id = string("symmetry/gapped/", _kgshort(model_T), "/Infinite")
        runner = function ()
            gap = Float64(fetch(model_T(), MassGap(), Infinite()))
            if gapped
                return CheckOutcome(
                    gap > _GAP_ATOL ? :pass : :fail,
                    gap,
                    _GAP_ATOL,
                    NaN,
                    NaN,
                    "profile declares gapped=true: fetched MassGap must exceed " *
                    "the $(_GAP_ATOL) floor",
                )
            else
                return _outcome(
                    gap,
                    0.0;
                    rtol=0.0,
                    atol=_GAP_ATOL,
                    detail="profile declares gapped=false: fetched MassGap must vanish",
                )
            end
        end
        push!(
            out,
            GeneratedCheck(
                :symmetry,
                id,
                string(
                    "symmetry profile of ",
                    _kgshort(model_T),
                    " declares gapped=",
                    gapped,
                    " — corroborated against the registered MassGap at Infinite",
                ),
                runner,
            ),
        )
    end
    return out
end

register_check_generator!(:symmetry, symmetry_checks)
