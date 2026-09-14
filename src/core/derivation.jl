# core/derivation.jl — the AbstractQAtlas relation NETWORK as a check generator.
#
# The other constraint edge types (identity.jl, bound.jl, duality.jl, limits.jl)
# each own a QAtlas-side STORE of hand-written edges.  This one does not: its
# store is AbstractQAtlas's relation registry, walked by `consistency_report`,
# which holds out each variable of a data set in turn and solves for it by every
# relation that reaches it from the rest.  What QAtlas supplies is the data and
# the scope; what the network supplies is which laws close over it.
#
# Scope of this file: the SCALING plane — a hub that fetches `CriticalExponents`
# is cross-checked against the `:scaling` domain (Rushbrooke, Widom, Fisher,
# Josephson).  `derivation_reach` below measures what the other planes would
# reach, and the measurement is why they are not here: 27 of 112 hubs close on
# some relation, 22 of those on `FreeEnergyLegendre` alone, which the
# `@identity_edge :gibbs` edge already covers.  Pinned by
# test/lint/test_derivation_reach.jl.
#
# What a passing row does and does not say.  On the scaling plane this is
# unusually sharp and worth stating exactly.  The `:scaling` algebra is four
# independent relations on seven numbers, so its solution set is
# three-dimensional, and those three are the renormalization parameters
# `(y_t, y_h, d)`: `ν = 1/y_t`, `β/ν = d − y_h`, `γ/ν = 2y_h − d`,
# `δ = y_h/(d − y_h)`, `η = d + 2 − 2y_h`, `α = 2 − dν`.  So a green row says the
# six exponents came from ONE fixed point's two eigenvalues, and says nothing
# about whether those eigenvalues are right.  What it catches is a value spliced
# in from a different source, which is what it found in the 3D percolation table.
# Pinned in test/core/test_derivation.jl, both directions.
#
# It is INTERNAL CONSISTENCY, the claim identity.jl makes, not corroboration
# against the literature; the literature plane is the `verify` cards with
# `route = :literature_value`.  Three ways the claim can be overstated are
# mechanised rather than left to prose:
#
#   * a value the table OBTAINED from one of these relations cannot check it.
#     `derived_from` names those routes and each is emitted as a visible :skip.
#   * a route absent for want of an INPUT is not a route excluded for want of
#     APPLICABILITY, and the network cannot tell them apart.  A hub whose
#     `dimension` is `nothing` therefore emits an explicit :skip for the
#     hyperscaling routes rather than letting them drop out silently.
#   * `consistency_report` on a NamedTuple keys on formula LETTERS, so a table
#     whose `α` is not the specific-heat exponent must be refused, not fed in.
#     `RefusedExponents` is that declaration; KPZ is the measured case.
#
# Tolerance is not a guessed rtol.  Each table carries its own quoted errors
# (`α_err`, …); a route's value carries those errors PROPAGATED through the
# relation, and the pass criterion is `|value − held| ≤ k·√(σ_held² + σ_route²)`.
# A table with an inexact participant and no quoted error for it is refused,
# never judged against a tolerance it did not state.

"""
    AbstractExponentSweep

A declared hub of the scaling plane: a `(model, bc)` whose `CriticalExponents`
row either IS handed to the network ([`SweptExponents`](@ref)) or is deliberately
kept out of it ([`RefusedExponents`](@ref)).

Two types rather than one with a `refused` field, for the reason
[`AbstractIdentityEdge`](@ref) gives: a refused hub has no sweep, no dimension and
no derived-from list, so one struct would carry a mode tag and a nothing-filled
half, and `_scaling_checks` would branch where it can dispatch.
"""
abstract type AbstractExponentSweep end

"""
    SweptExponents

A hub whose exponent table is cross-checked: the `(model, bc)`, the `sweep` of
fetch kwargs to run it at, the spatial `dimension` the hyperscaling relations
need, and the routes that may not judge it.

`dimension` is `nothing` where hyperscaling does not apply, `:sweep` where the
sweep's own `d` is the dimension, or a number where the hub sits at a fixed `d`
its `fetch` does not take as a kwarg.

`derived_from` maps a target exponent to the relation NAMES the shipped value was
obtained from; those routes cannot check it and are skipped with a reason.  `k`
is the sigma multiplier of the pass criterion.
"""
struct SweptExponents <: AbstractExponentSweep
    model::Type
    bc::Type
    sweep::NamedTuple
    dimension::Union{Nothing,Symbol,Real}
    derived_from::Vector{Pair{Symbol,Vector{Symbol}}}
    k::Float64
    notes::String
    references::Vector{String}
end

"""
    RefusedExponents

A hub whose `CriticalExponents` are kept OUT of the network, carrying `reason`.

It exists because the name-keyed `consistency_report` cannot tell two quantities
wearing one letter apart, so a table whose `α` is not the specific-heat exponent
has to be excluded rather than fed in and reported as a contradiction.  The
refusal is emitted as a visible skip, because an absence and an exclusion read
the same in a pass count.
"""
struct RefusedExponents <: AbstractExponentSweep
    model::Type
    bc::Type
    reason::String
    references::Vector{String}
end

"""
    EXPONENT_SWEEPS :: Vector{AbstractExponentSweep}

The declared scaling-plane hubs, populated at include-time by
[`exponent_sweep!`](@ref) / [`refuse_exponents!`](@ref).  Every hub with a
`CriticalExponents` method must appear here, or
[`check_derivation_coverage`](@ref) reports it: an exponent table nothing
declares is a table the network never sees.
"""
const EXPONENT_SWEEPS = AbstractExponentSweep[]

# The relation names the `:scaling` domain actually offers, read from the
# registry rather than listed.  `derived_from` is a second copy of a provenance
# header, and an unvalidated second copy fails in the direction that looks fine:
# a mistyped `:widom` never matches, the route it was meant to suppress runs, and
# it PASSES, because the shipped value is what that relation predicts. The
# anti-tautology guarantee would then be a tautology.
function _scaling_relation_names()
    return Set(
        Symbol(nameof(typeof(r))) for r in AbstractQAtlas.all_relations(; domain=:scaling)
    )
end

"""
    exponent_sweep!(model, bc; sweep=(;), dimension=:sweep, derived_from=[],
                    k=3.0, notes="", references=String[])

Declare a [`SweptExponents`](@ref) hub.  See `src/derivation_registry.jl` for the
catalog.

Refuses a duplicate `(model, bc, sweep)`, a non-positive `k`, a `dimension` that
is neither `nothing`, `:sweep` nor a real, a `:sweep` dimension whose sweep
carries no real `d`, and a `derived_from` naming an exponent or a relation that
does not exist.  `dimension=nothing` needs `notes`: it suppresses the
hyperscaling routes, and an unexplained suppression reads as coverage.
"""
function exponent_sweep!(
    model::Type,
    bc::Type;
    sweep::NamedTuple=NamedTuple(),
    dimension::Union{Nothing,Symbol,Real}=:sweep,
    derived_from::AbstractVector=Pair{Symbol,Vector{Symbol}}[],
    k::Real=3.0,
    notes::AbstractString="",
    references::AbstractVector{<:AbstractString}=String[],
)
    _reject_duplicate_hub(model, bc, sweep)
    k > 0 || throw(ArgumentError("exponent_sweep!: k must be > 0; got $(k)"))
    dimension isa Symbol &&
        dimension !== :sweep &&
        throw(
            ArgumentError(
                "exponent_sweep!: dimension is `nothing`, `:sweep` or a number; " *
                "got :$(dimension)",
            ),
        )
    if dimension === :sweep
        haskey(sweep, :d) || throw(
            ArgumentError(
                "exponent_sweep!: dimension=:sweep needs the sweep to carry `d` " *
                "(the algebra's d must be the d the exponents were fetched at)",
            ),
        )
        # Checked here rather than left to the solver: `d` is the one value that
        # reaches the relations without passing the `isa Real` test the fetched
        # exponents get, so a non-numeric one surfaces as a MethodError from inside
        # AbstractQAtlas instead of naming the declaration that caused it.
        all(x -> x isa Real, sweep.d) || throw(
            ArgumentError(
                "exponent_sweep!: every swept `d` must be a real number; got $(sweep.d)"
            ),
        )
    end
    dimension === nothing &&
        isempty(notes) &&
        throw(
            ArgumentError(
                "exponent_sweep!: dimension=nothing suppresses the hyperscaling " *
                "routes, so it needs `notes` saying why they do not apply",
            ),
        )
    known = _scaling_relation_names()
    df = Pair{Symbol,Vector{Symbol}}[]
    for p in derived_from
        target = Symbol(first(p))
        target in _SCALING_EXPONENTS || throw(
            ArgumentError(
                "exponent_sweep!: derived_from names :$(target), which is not one of " *
                "$(_SCALING_EXPONENTS)",
            ),
        )
        rels = Symbol[last(p)...]
        for r in rels
            r in known || throw(
                ArgumentError(
                    "exponent_sweep!: derived_from names the relation :$(r), which is " *
                    "not in the :scaling domain — a name that matches nothing " *
                    "suppresses nothing, and the route it meant to exclude would pass",
                ),
            )
        end
        push!(df, Pair{Symbol,Vector{Symbol}}(target, rels))
    end
    push!(
        EXPONENT_SWEEPS,
        SweptExponents(
            model,
            bc,
            sweep,
            dimension,
            df,
            Float64(k),
            String(notes),
            String[r for r in references],
        ),
    )
    return nothing
end

"""
    refuse_exponents!(model, bc; reason, references=String[])

Declare a [`RefusedExponents`](@ref) hub: its `CriticalExponents` are kept out of
the network and `reason` says why.

A refusal is keyed on `(model, bc)` ALONE, because its generated check carries no
sweep point in its id; two refusals for one hub would pass a sweep-aware
duplicate test and then collide on that id, which surfaces three layers away as
`generated_checks` calling the generator non-deterministic.
"""
function refuse_exponents!(
    model::Type,
    bc::Type;
    reason::AbstractString,
    references::AbstractVector{<:AbstractString}=String[],
)
    isempty(reason) &&
        throw(ArgumentError("refuse_exponents!: a refusal carries its reason"))
    _reject_duplicate_hub(model, bc, nothing)
    push!(
        EXPONENT_SWEEPS,
        RefusedExponents(model, bc, String(reason), String[r for r in references]),
    )
    return nothing
end

# A swept hub is identified by `(model, bc, sweep)` — Ising is legitimately
# declared twice, at different `d`, and the sweep is what separates their check
# ids.  A refusal (`sweep === nothing`) is identified by `(model, bc)`, and
# collides with any other declaration of the same hub in either direction.
function _reject_duplicate_hub(model::Type, bc::Type, sweep)
    for s in EXPONENT_SWEEPS
        (s.model === model && s.bc === bc) || continue
        (sweep !== nothing && s isa SweptExponents && s.sweep != sweep) && continue
        throw(
            ArgumentError(
                "declaring $(_kgshort(model))/$(_kgshort(bc)): already declared as " *
                "$(nameof(typeof(s)))",
            ),
        )
    end
    return nothing
end

"""
    @exponent_sweep Model BC key=value …

Macro sugar around [`exponent_sweep!`](@ref), matching `@identity_edge`'s shape:
the two positional arguments are the hub, the rest are forwarded as keywords.
"""
macro exponent_sweep(model, bc, kwargs...)
    return _forward_kw_macro(exponent_sweep!, :exponent_sweep, (model, bc), kwargs)
end

"""
    @refuse_exponents Model BC reason=…

Macro sugar around [`refuse_exponents!`](@ref).
"""
macro refuse_exponents(model, bc, kwargs...)
    return _forward_kw_macro(refuse_exponents!, :refuse_exponents, (model, bc), kwargs)
end

# The six exponents the `:scaling` domain is written on.  `c`, `ψ`, `x_m` and the
# `_err` companions a table also carries are not variables of those relations.
const _SCALING_EXPONENTS = (:α, :β, :γ, :δ, :ν, :η)

# The routes that need a spatial dimension.  Named rather than detected, so a hub
# that suppresses them says so and the skip is visible; a new hyperscaling law
# upstream is added here with its reason.
const _HYPERSCALING_RELATIONS = (:Josephson, :QuantumHyperscaling)

# `(values, errors)` for one fetched exponent table, or a refusal string.
# An exact participant (Rational/Integer) needs no quoted error; an inexact one
# does, because its route's tolerance is that error propagated and there is no
# honest default for a number whose precision the table declined to state.
function _scaling_data(nt::NamedTuple, dim)
    ks = Tuple(k for k in _SCALING_EXPONENTS if haskey(nt, k))
    length(ks) ≥ 3 ||
        return "table carries $(length(ks)) of the six scaling exponents, too few to close"
    vals = Any[]
    errs = Any[]
    for k in ks
        v = nt[k]
        v isa Real || return "$(k) is $(typeof(v)), not a real number"
        push!(vals, v)
        if v isa Union{Rational,Integer}
            push!(errs, 0.0)
        else
            e = get(nt, Symbol(k, "_err"), nothing)
            e === nothing && return "$(k) is inexact and the table quotes no $(k)_err"
            push!(errs, Float64(e))
        end
    end
    names = ks
    if dim !== nothing
        names = (names..., :d)
        push!(vals, dim)
        push!(errs, 0.0)
    end
    return (NamedTuple{names}(Tuple(vals)), NamedTuple{names}(Tuple(errs)))
end

function _solve_step(step, data::NamedTuple)
    return solve(step.relation, Val(step.output); (v => data[v] for v in step.inputs)...)
end

# The route value's own uncertainty: each input's quoted error pushed through the
# relation, added in quadrature.  Linear propagation by finite difference AT the
# quoted error, exact for a relation affine in that input, which the `:scaling`
# relations are.
#
# A non-finite contribution returns `NaN` rather than entering the sum.  These
# routes divide by a variable (Widom by `β` and by `δ−1`, Fisher by `ν` and by
# `2−η`, Josephson by `ν`), so a quoted error that straddles a pole would return
# an infinite sigma, and an infinite sigma is an infinite tolerance that accepts
# every disagreement.  `_sigma_outcome` refuses the `NaN` instead of passing it.
# Upstream's own `_families_satisfied` drops non-finite scales for this reason.
function _route_sigma(step, data::NamedTuple, errs::NamedTuple, base::Real)
    b = Float64(base)
    isfinite(b) || return NaN
    s2 = 0.0
    for v in step.inputs
        e = getfield(errs, v)
        iszero(e) && continue
        bumped = merge(data, NamedTuple{(v,)}((data[v] + e,)))
        δ = Float64(_solve_step(step, bumped)) - b
        isfinite(δ) || return NaN
        s2 += abs2(δ)
    end
    return sqrt(s2)
end

# Round-off floor, so an exact table is not asked for bit equality through a
# solver that may have promoted its rationals.  Machine precision scaled to the
# compared magnitudes; never a physics tolerance.
function _roundoff_floor(a::Real, b::Real)
    return 8 * eps(Float64) * max(abs(Float64(a)), abs(Float64(b)), 1.0)
end

"""
    _sigma_outcome(held, value, sigma_held, sigma_route; k, detail="") -> CheckOutcome

The scaling plane's pass criterion: `|value − held| ≤ k·√(σ_held² + σ_route²)`,
floored at round-off.  Reported as `lhs = held`, `rhs = value`, with `rel_err`
carrying the DEVIATION IN SIGMA rather than a relative difference — what a reader
of a failing row needs is how far outside the stated errors it is.

A non-finite input is `:error`, never `:pass`.  The tolerance is built FROM the
sigmas, so an infinite one accepts every disagreement; refusing is the only
answer that does not turn a numerical breakdown into agreement.
"""
function _sigma_outcome(
    held::Real, value::Real, sigma_held::Real, sigma_route::Real; k::Real, detail::String=""
)
    l, r = Float64(held), Float64(value)
    sh, sr = Float64(sigma_held), Float64(sigma_route)
    if !(isfinite(l) && isfinite(r) && isfinite(sh) && isfinite(sr))
        return CheckOutcome(
            :error,
            l,
            r,
            NaN,
            NaN,
            "no tolerance can be formed: held=$(l) value=$(r) σ_held=$(sh) σ_route=$(sr)",
        )
    end
    abs_err = abs(l - r)
    sigma = hypot(sh, sr)
    tol = max(k * sigma, _roundoff_floor(l, r))
    status = abs_err ≤ tol ? :pass : :fail
    n_sigma = sigma > 0 ? abs_err / sigma : (status === :pass ? 0.0 : Inf)
    return CheckOutcome(status, l, r, abs_err, n_sigma, detail)
end

# The routes a sweep declares circular for `target`, plus the hyperscaling ones
# when the hub has no dimension.  `name => reason`, so the skip says which.
function _suppressed_routes(s::SweptExponents, target::Symbol)
    out = Pair{Symbol,String}[]
    for (t, rels) in s.derived_from
        t === target || continue
        for r in rels
            push!(
                out,
                r =>
                    "the shipped $(target) was obtained from $(r); a route cannot " *
                    "check the value it produced",
            )
        end
    end
    return out
end

# The hyperscaling routes of a hub that has no dimension, as STANDALONE skips.
#
# They cannot be suppressed route-by-route the way `derived_from` is: with no
# dimension, `d` never enters the data, so `consistency_report` never reaches a
# relation that needs it and there is no route to mark.  That silent absence is
# the very thing upstream warns about — a route missing for want of an input
# reads exactly like one excluded for want of applicability — so the exclusion is
# emitted on its own rather than inferred from a gap.
function _push_hyperscaling_skips!(out, hub::AbstractString, s::SweptExponents)
    s.dimension === nothing || return out
    for r in _HYPERSCALING_RELATIONS
        _push_excluded_check!(
            out,
            :derivation,
            string(hub, "/", r),
            "hyperscaling does not apply here, so $(r) is not a route: $(s.notes)",
        )
    end
    return out
end

# ──────────────────────────────────────────────────────────────────────
# Generator — the :derivation kind of generated_checks()
# ──────────────────────────────────────────────────────────────────────

"""
    ScalingRefusal(status, reason)

Why a hub produced no routes at a sweep point.  `status` is `:skip` where the
table is DECLARED not to close (too few exponents, no quoted error, no relation
reaches it) and `:error` where a call THREW.

The split is the point.  A renamed `fetch` keyword and a BKT table carrying η
alone both end a sweep point, and reported the same way the first disappears into
the second: a suite that treats every non-route as a declared skip goes green
while a whole hub's cross-checks stop existing.  `:error` fails the suite,
matching [`CheckOutcome`](@ref)'s own split of a config bug from a contradiction.
"""
struct ScalingRefusal
    status::Symbol
    reason::String
end

# Fetch the table, shape it for the network, and ask the network which routes
# reach each exponent.  Returns `(data, errs, rows)` or a [`ScalingRefusal`](@ref).
# The fetch happens HERE rather than inside each check's runner: one table feeds
# every route of the hub, and re-fetching per route would let a stale-value bug
# hide behind rows that each fetched their own copy.
function _prepare_scaling(s::SweptExponents, point::NamedTuple, dim)
    table = try
        fetch(s.model(), CriticalExponents(), _bc_instance(s.bc; finite_N=8); point...)
    catch err
        return ScalingRefusal(:error, "fetch threw: $(sprint(showerror, err))")
    end
    table isa NamedTuple ||
        return ScalingRefusal(:error, "fetch returned $(typeof(table)), not a NamedTuple")
    prepared = _scaling_data(table, dim)
    prepared isa String && return ScalingRefusal(:skip, prepared)
    data, errs = prepared
    rows = try
        AbstractQAtlas.consistency_report(data; domain=:scaling, atol=0, rtol=0)
    catch err
        return ScalingRefusal(:error, "consistency_report threw: $(sprint(showerror, err))")
    end
    # `agree` is deliberately not read: that verdict is one tolerance for a whole
    # row, and this plane judges each route against ITS OWN propagated error.
    # What the report is used for here is the enumeration of routes.
    isempty(rows) && return ScalingRefusal(
        :skip, "no :scaling relation reaches any exponent of this table"
    )
    return (data, errs, rows)
end

# A refusal is reported as a skip, never as an absence: `_push_excluded_check!`
# for the declared kind, and a check that reports `:error` where a call THREW,
# so a renamed `fetch` keyword fails the suite instead of joining the several
# skips this table legitimately carries.
function _push_refusal_check!(out, id::AbstractString, r::ScalingRefusal)
    r.status === :skip && return _push_excluded_check!(out, :derivation, id, r.reason)
    push!(
        out,
        GeneratedCheck(
            :derivation,
            id,
            "BROKEN: $(r.reason)",
            () -> CheckOutcome(:error, NaN, NaN, NaN, NaN, r.reason),
        ),
    )
    return out
end

function _hub_id(s::AbstractExponentSweep)
    return string("derivation/scaling/", _kgshort(s.model), "/", _kgshort(s.bc))
end

function _scaling_checks(s::RefusedExponents)
    out = GeneratedCheck[]
    _push_excluded_check!(out, :derivation, _hub_id(s), s.reason)
    return out
end

function _scaling_checks(s::SweptExponents)
    out = GeneratedCheck[]
    hub = _hub_id(s)
    _push_hyperscaling_skips!(out, hub, s)
    for point in _sweep_points(s.sweep)
        pid = _point_suffix(point)
        dim = s.dimension === :sweep ? getfield(point, :d) : s.dimension
        prepared = _prepare_scaling(s, point, dim)
        if prepared isa ScalingRefusal
            _push_refusal_check!(out, hub * pid, prepared)
            continue
        end
        data, errs, rows = prepared
        for row in rows
            target = row.target
            suppressed = _suppressed_routes(s, target)
            for (step, value) in zip(row.steps, row.values)
                rname = Symbol(nameof(typeof(step.relation)))
                id = string(hub, pid, "/", target, "/", rname)
                hit = findfirst(p -> first(p) === rname, suppressed)
                if hit !== nothing
                    _push_excluded_check!(out, :derivation, id, last(suppressed[hit]))
                    continue
                end
                held = row.held_out
                sigma_held = getfield(errs, target)
                kk = s.k
                desc = string(
                    ":",
                    target,
                    " of ",
                    _kgshort(s.model),
                    isempty(pid) ? "" : " ($(_point_id(point)))",
                    " against ",
                    rname,
                    " on the rest of the table",
                )
                push!(
                    out,
                    GeneratedCheck(
                        :derivation,
                        id,
                        desc,
                        # The propagation runs INSIDE the runner. Computed while
                        # generating, a throwing `solve` would escape
                        # `run_generated_check` and abort `generated_checks()` for
                        # every kind, not just this route.
                        function ()
                            sigma_route = _route_sigma(step, data, errs, value)
                            return _sigma_outcome(
                                held, value, sigma_held, sigma_route; k=kk
                            )
                        end,
                    ),
                )
            end
        end
    end
    return out
end

function derivation_checks()
    out = GeneratedCheck[]
    for s in EXPONENT_SWEEPS
        append!(out, _scaling_checks(s))
    end
    return out
end

register_check_generator!(:derivation, derivation_checks)

# ──────────────────────────────────────────────────────────────────────
# Coverage — an undeclared exponent table is a table the network never sees
# ──────────────────────────────────────────────────────────────────────

"""
    exponent_hubs() -> Vector{Type}

Every model type with a `fetch(model, ::CriticalExponents, …)` method, found by
scanning `methods(fetch)`.

Reflection rather than [`REGISTRY`](@ref), because most universality classes
carry no registry row: of the fifteen types with such a method, five are
registered (`Universality{:MeanField}` plus the four delegating models).  A
coverage guard keyed on the registry alone would be structurally unable to see
the ten classes that carry the exponent tables.
"""
function exponent_hubs()
    out = Type[]
    for (M, _) in _exponent_methods()
        M === nothing && continue
        M in out || push!(out, M)
    end
    return sort!(out; by=_kgshort)
end

# Every `CriticalExponents` method as `(hub_or_nothing, method)`.  `nothing` where
# the model slot is a TypeVar: a `fetch(m::M, ::CriticalExponents, …) where {M<:…}`
# has no hub to name, so it cannot be declared and cannot be exempted, and
# dropping it would let it read as "no gap" instead of "a gap nothing can close".
function _exponent_methods()
    out = Tuple{Union{Nothing,Type},Method}[]
    for m in methods(fetch)
        p = Base.unwrap_unionall(m.sig).parameters
        length(p) ≥ 3 || continue
        M, Q = p[2], p[3]
        (Q isa Type && Q <: CriticalExponents) || continue
        hub = if M isa Type && isconcretetype(M) && M <: AbstractQAtlasModel
            M
        else
            nothing
        end
        push!(out, (hub, m))
    end
    return out
end

"""
    check_derivation_coverage() -> Vector{CoherenceFinding}

Every hub that can fetch `CriticalExponents` is declared — swept or refused —
and every declaration generates at least one check.

Both directions are needed and they fail differently: an undeclared hub is a
table nothing cross-checks and nothing says so, and a declaration that judges
nothing reads as coverage while constraining nothing.  Delegation rows of
[`REGISTRY`](@ref) are exempt: a delegated table is the same numbers as the hub
it routes to, so sweeping it would run one check several times under different
names.
"""
function check_derivation_coverage()
    out = CoherenceFinding[]
    declared = Set(s.model for s in EXPONENT_SWEEPS)
    delegating = Set(
        e.model for
        e in REGISTRY if e.quantity === CriticalExponents && _is_delegation(e.method)
    )
    for (M, m) in _exponent_methods()
        if M === nothing
            push!(
                out,
                CoherenceFinding(
                    :derivation_coverage,
                    :gap,
                    "a CriticalExponents method with a non-concrete model slot " *
                    "($(m.file):$(m.line)) names no hub, so it can be neither " *
                    "declared nor exempted — the scaling algebra cannot reach it",
                ),
            )
            continue
        end
        (M in declared || M in delegating) && continue
        push!(
            out,
            CoherenceFinding(
                :derivation_coverage,
                :gap,
                "$(_kgshort(M)) fetches CriticalExponents but declares no " *
                "exponent_sweep! — the scaling algebra never sees it",
            ),
        )
    end
    for s in EXPONENT_SWEEPS
        _emits_only_skips(s) && push!(
            out,
            CoherenceFinding(
                :derivation_coverage,
                :gap,
                "exponent_sweep! $(_kgshort(s.model))/$(_kgshort(s.bc)) at $(s.sweep) " *
                "emits only skips — it judges nothing",
            ),
        )
    end
    return out
end

"""
    _emits_only_skips(s::AbstractExponentSweep) -> Bool

Whether a SWEPT declaration produces no check that can return a verdict.

Emptiness is the test that reads naturally here and it cannot fire: a point
`_scaling_checks` cannot prepare still emits a refusal check, so the vector is
never empty and a guard on `isempty` would be a guard unable to fail.  What can
happen, and is what this asks, is a declaration every one of whose points was
skipped.  An `:error` counts as content: it fails the suite, so a hub that is
merely broken is loud already and does not also need a coverage finding.  A
refused hub is exempt by construction; its skip IS its content.
"""
_emits_only_skips(::RefusedExponents) = false
function _emits_only_skips(s::SweptExponents)
    return !any(c -> run_generated_check(c).status !== :skip, _scaling_checks(s))
end

# ──────────────────────────────────────────────────────────────────────
# Reach — what the OTHER planes would cross-check, measured not guessed
# ──────────────────────────────────────────────────────────────────────

"""
    DerivationReach

One `(model, bc)` hub and the relation names whose type-keyed derivation step is
CLOSED on it: the hub fetches the step's output and every quantity the step
needs.  `quantities` is how many distinct quantity families the hub implements.
"""
struct DerivationReach
    model::Type
    bc::Type
    quantities::Int
    relations::Vector{Symbol}
end

# An `InverseTemperature` slot is a fetch kwarg here, not a second fetched value,
# so a step whose only other typed slot is the temperature relates one atlas
# number to itself and cross-checks nothing.
function _is_supplied_slot(@nospecialize(T::Type))
    return Base.typename(T).wrapper in
           (AbstractQAtlas.InverseTemperature, AbstractQAtlas.Temperature)
end

"""
    derivation_reach() -> Vector{DerivationReach}

Which hubs the AbstractQAtlas relation network could cross-check today: for each
`(model, bc)` in [`REGISTRY`](@ref), the relations whose type-keyed derivation
step has its output and all of its quantity inputs implemented there.

A genuine cross-check needs the hub to fetch at least TWO of the step's
quantities, so a step whose only other typed slot is a supplied temperature does
not count.

Structural, and an upper bound on what fires, in two ways worth keeping apart. A
step's UNTYPED slots are not consulted here at all, so a relation can close on a
hub's quantities and still have no route because a supplied value it needs does
not exist — four of the reachable relations sit in `test_abq_conformance.jl`'s
`MATERIALIZABLE_BUT_UNWIRED` for exactly that. And whether the solve computes at
all is [`generated_checks`](@ref)'s business, not this one's.

Name-sorted and deterministic.  Pinned by `test/lint/test_derivation_reach.jl`,
so the reachable set cannot shrink unnoticed.
"""
function derivation_reach()
    steps = AbstractQAtlas.typed_derivation_steps()
    hubs = Dict{Tuple{Type,Type},Set{Any}}()
    for e in REGISTRY
        e.method === :not_implemented && continue
        push!(get!(hubs, (e.model, e.bc), Set{Any}()), Base.typename(e.quantity).wrapper)
    end
    out = DerivationReach[]
    for ((M, BC), qs) in hubs
        names = Symbol[]
        for st in steps
            Base.typename(st.output.type).wrapper in qs || continue
            all(
                i -> _is_supplied_slot(i.type) || Base.typename(i.type).wrapper in qs,
                st.inputs,
            ) || continue
            count(i -> !_is_supplied_slot(i.type), st.inputs) ≥ 1 || continue
            push!(names, Symbol(nameof(typeof(st.relation))))
        end
        isempty(names) && continue
        push!(out, DerivationReach(M, BC, length(qs), sort!(unique!(names))))
    end
    sort!(out; by=r -> (_kgshort(r.model), _kgshort(r.bc)))
    return out
end
