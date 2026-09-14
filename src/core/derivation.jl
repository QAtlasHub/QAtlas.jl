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
# reach; the number is small and is pinned by test/lint/test_derivation_reach.jl.
#
# What a passing row does and does not say.  This is INTERNAL CONSISTENCY, the
# claim identity.jl makes, not corroboration against the literature: the routes
# and the held-out value come from one table.  Three ways that claim can be
# overstated are mechanised rather than left to prose:
#
#   * a value the table OBTAINED from one of these relations cannot check it.
#     `derived_from` names those routes and each is emitted as a visible :skip.
#   * a route absent for want of an INPUT is not a route excluded for want of
#     APPLICABILITY, and the network cannot tell them apart.  A hub whose
#     `dimension` is `nothing` therefore emits an explicit :skip for the
#     hyperscaling routes rather than letting them drop out silently.
#   * `consistency_report` on a NamedTuple keys on formula LETTERS, so a table
#     whose `α` is not the specific-heat exponent must be refused, not fed in.
#     `refused` is that declaration; KPZ is the measured case.
#
# Tolerance is not a guessed rtol.  Each table carries its own quoted errors
# (`α_err`, …); a route's value carries those errors PROPAGATED through the
# relation, and the pass criterion is `|value − held| ≤ k·√(σ_held² + σ_route²)`.
# A table with an inexact participant and no quoted error for it is refused,
# never judged against a tolerance it did not state.

"""
    ExponentSweep

One declared hub of the scaling plane: the `(model, bc)` whose
`CriticalExponents` row is cross-checked, the `sweep` of fetch kwargs to run it
at, the spatial `dimension` the hyperscaling relations need, and the routes that
may not judge it.

`dimension` is `nothing` where hyperscaling does not apply, `:sweep` where the
sweep's own `d` is the dimension, or a number where the hub sits at a fixed `d`
its `fetch` does not take as a kwarg.

`derived_from` maps a target exponent to the relation NAMES the shipped value
was obtained from; those routes cannot check it and are skipped with a reason.
`refused`, when set, makes the whole hub one visible skip carrying that reason.
`k` is the sigma multiplier of the pass criterion.
"""
struct ExponentSweep
    model::Type
    bc::Type
    sweep::NamedTuple
    dimension::Union{Nothing,Symbol,Real}
    derived_from::Vector{Pair{Symbol,Vector{Symbol}}}
    refused::Union{Nothing,String}
    k::Float64
    notes::String
    references::Vector{String}
end

"""
    EXPONENT_SWEEPS :: Vector{ExponentSweep}

The declared scaling-plane hubs, populated at include-time by
[`exponent_sweep!`](@ref).  Every hub with a `CriticalExponents` method must
appear here — swept or refused — or [`check_derivation_coverage`](@ref) reports
it: an exponent table nothing declares is a table the network never sees.
"""
const EXPONENT_SWEEPS = ExponentSweep[]

"""
    exponent_sweep!(model, bc; sweep=(;), dimension=:sweep, derived_from=[],
                    refused=nothing, k=3.0, notes="", references=String[])

Declare one scaling-plane hub.  See [`ExponentSweep`](@ref) for the fields and
`src/derivation_registry.jl` for the catalog.

Refuses a duplicate `(model, bc, sweep)`, a non-positive `k`, a `dimension` that
is neither `nothing`, `:sweep` nor a real, and a `:sweep` dimension whose sweep
carries no `d` (the algebra's `d` must be the `d` the exponents were fetched at).
`dimension=nothing` and `refused` each require `notes`, because both suppress
checks and an unexplained suppression reads as coverage.
"""
function exponent_sweep!(
    model::Type,
    bc::Type;
    sweep::NamedTuple=NamedTuple(),
    dimension::Union{Nothing,Symbol,Real}=:sweep,
    derived_from::AbstractVector=Pair{Symbol,Vector{Symbol}}[],
    refused::Union{Nothing,AbstractString}=nothing,
    k::Real=3.0,
    notes::AbstractString="",
    references::AbstractVector{<:AbstractString}=String[],
)
    any(s -> s.model === model && s.bc === bc && s.sweep == sweep, EXPONENT_SWEEPS) &&
        throw(
            ArgumentError(
                "exponent_sweep!: $(_kgshort(model))/$(_kgshort(bc)) at $(sweep) is " *
                "already declared",
            ),
        )
    k > 0 || throw(ArgumentError("exponent_sweep!: k must be > 0; got $(k)"))
    if refused !== nothing
        isempty(notes) || throw(
            ArgumentError(
                "exponent_sweep!: a refused hub carries its reason in " *
                "`refused`; drop `notes`",
            ),
        )
        push!(
            EXPONENT_SWEEPS,
            ExponentSweep(
                model,
                bc,
                sweep,
                nothing,
                Pair{Symbol,Vector{Symbol}}[],
                String(refused),
                Float64(k),
                "",
                String[r for r in references],
            ),
        )
        return nothing
    end
    dimension isa Symbol &&
        dimension !== :sweep &&
        throw(
            ArgumentError(
                "exponent_sweep!: dimension is `nothing`, `:sweep` or a number; " *
                "got :$(dimension)",
            ),
        )
    dimension === :sweep &&
        !haskey(sweep, :d) &&
        throw(
            ArgumentError(
                "exponent_sweep!: dimension=:sweep needs the sweep to carry `d` " *
                "(the algebra's d must be the d the exponents were fetched at)",
            ),
        )
    dimension === nothing &&
        isempty(notes) &&
        throw(
            ArgumentError(
                "exponent_sweep!: dimension=nothing suppresses the hyperscaling " *
                "routes, so it needs `notes` saying why they do not apply",
            ),
        )
    df = Pair{Symbol,Vector{Symbol}}[]
    for p in derived_from
        push!(df, Pair{Symbol,Vector{Symbol}}(Symbol(first(p)), Symbol[last(p)...]))
    end
    push!(
        EXPONENT_SWEEPS,
        ExponentSweep(
            model,
            bc,
            sweep,
            dimension,
            df,
            nothing,
            Float64(k),
            String(notes),
            String[r for r in references],
        ),
    )
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
function _route_sigma(step, data::NamedTuple, errs::NamedTuple, base::Real)
    s2 = 0.0
    for v in step.inputs
        e = getfield(errs, v)
        iszero(e) && continue
        bumped = merge(data, NamedTuple{(v,)}((data[v] + e,)))
        s2 += abs2(Float64(_solve_step(step, bumped)) - Float64(base))
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
"""
function _sigma_outcome(
    held::Real, value::Real, sigma_held::Real, sigma_route::Real; k::Real, detail::String=""
)
    l, r = Float64(held), Float64(value)
    abs_err = abs(l - r)
    sigma = sqrt(abs2(Float64(sigma_held)) + abs2(Float64(sigma_route)))
    tol = max(k * sigma, _roundoff_floor(l, r))
    status = abs_err ≤ tol ? :pass : :fail
    n_sigma = sigma > 0 ? abs_err / sigma : (status === :pass ? 0.0 : Inf)
    return CheckOutcome(status, l, r, abs_err, n_sigma, detail)
end

# The routes a sweep declares circular for `target`, plus the hyperscaling ones
# when the hub has no dimension.  `name => reason`, so the skip says which.
function _suppressed_routes(s::ExponentSweep, target::Symbol)
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
    if s.dimension === nothing
        for r in _HYPERSCALING_RELATIONS
            push!(out, r => "hyperscaling does not apply here: $(s.notes)")
        end
    end
    return out
end

# ──────────────────────────────────────────────────────────────────────
# Generator — the :derivation kind of generated_checks()
# ──────────────────────────────────────────────────────────────────────

# Fetch the table, shape it for the network, and ask the network which routes
# reach each exponent.  Returns `(data, errs, rows)` or a refusal string.
# The fetch happens HERE rather than inside each check's runner: one table feeds
# every route of the hub, and re-fetching per route would let a stale-value bug
# hide behind rows that each fetched their own copy.
function _prepare_scaling(s::ExponentSweep, point::NamedTuple, dim)
    table = try
        fetch(s.model(), CriticalExponents(), _bc_instance(s.bc; finite_N=8); point...)
    catch err
        return "fetch threw: $(sprint(showerror, err))"
    end
    table isa NamedTuple || return "fetch returned $(typeof(table)), not a NamedTuple"
    prepared = _scaling_data(table, dim)
    prepared isa String && return prepared
    data, errs = prepared
    rows = try
        AbstractQAtlas.consistency_report(data; domain=:scaling, atol=0, rtol=0)
    catch err
        return "consistency_report threw: $(sprint(showerror, err))"
    end
    # `agree` is deliberately not read: that verdict is one tolerance for a whole
    # row, and this plane judges each route against ITS OWN propagated error.
    # What the report is used for here is the enumeration of routes.
    isempty(rows) && return "no :scaling relation reaches any exponent of this table"
    return (data, errs, rows)
end

function _scaling_checks(s::ExponentSweep)
    out = GeneratedCheck[]
    hub = string("derivation/scaling/", _kgshort(s.model), "/", _kgshort(s.bc))
    if s.refused !== nothing
        _push_excluded_check!(out, :derivation, hub, s.refused)
        return out
    end
    for point in _sweep_points(s.sweep)
        pid = isempty(keys(point)) ? "" : "/" * _point_id(point)
        dim = s.dimension === :sweep ? getfield(point, :d) : s.dimension
        prepared = _prepare_scaling(s, point, dim)
        if prepared isa String
            _push_excluded_check!(out, :derivation, hub * pid, prepared)
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
                sigma_route = _route_sigma(step, data, errs, value)
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
                        () -> _sigma_outcome(held, value, sigma_held, sigma_route; k=kk),
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
    for m in methods(fetch)
        p = Base.unwrap_unionall(m.sig).parameters
        length(p) ≥ 3 || continue
        M, Q = p[2], p[3]
        (M isa Type && Q isa Type) || continue
        isconcretetype(M) && M <: AbstractQAtlasModel || continue
        Q <: CriticalExponents || continue
        M in out || push!(out, M)
    end
    return sort!(out; by=_kgshort)
end

"""
    check_derivation_coverage() -> Vector{CoherenceFinding}

Every hub that can fetch `CriticalExponents` is declared — swept or refused —
and every declaration generates at least one check.

Both directions are needed and they fail differently: an undeclared hub is a
table nothing cross-checks and nothing says so, and a declaration that generates
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
    for M in exponent_hubs()
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
        isempty(_scaling_checks(s)) && push!(
            out,
            CoherenceFinding(
                :derivation_coverage,
                :gap,
                "exponent_sweep! $(_kgshort(s.model))/$(_kgshort(s.bc)) at $(s.sweep) " *
                "generates no checks — it constrains nothing",
            ),
        )
    end
    return out
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
not count.  Structural, and an upper bound on what fires: whether the solve
actually computes is [`generated_checks`](@ref)'s business, not this one's.

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
