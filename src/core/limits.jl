# core/limits.jl — asymptotic limit edges (#701, over the #697 kernel).
#
# A limit edge declares that the `source` family approaches the `target`
# family as one source parameter is driven to a limit — distinct from
# `@reduces`, which records EXACT delegation at a point in parameter space.
# `XXZ1D` at Δ = 1 *reduces to* Heisenberg1D (equality); `XXZ1D` as Δ → 1⁺
# *limits to* it (convergence along a sequence).  Asymptotic edges document,
# in the graph, which models are strong-coupling / isotropic / classical
# shadows of which, and extend selective CI's blast radius accordingly.
#
# The generated-test shape (hand-written precedent: the Δ → 1⁺ two-spinon DSF
# check of #691) is a CONVERGENCE-SEQUENCE check: evaluate the quantity on the
# source at each point of the declared `approach` sequence, compare against
# the target value, and require the error to (weakly) shrink along the
# sequence and to land below `final_atol` at its end.  Both endpoints must be
# independent (non-delegating) implementations — a delegating endpoint would
# make the sequence converge to itself by construction (#701's circularity
# rule, enforced via the kernel's independence filter).

"""
    LimitQuantitySpec

One quantity a [`LimitEdge`](@ref)'s generated convergence check runs on:
compared at boundary condition `bc` on the fetch-kwargs grid `sweep`, with the
sequence's terminal error required below `final_atol` (set from the measured
convergence at declaration time; the looser early-sequence behaviour is
covered by the monotonicity requirement).
"""
struct LimitQuantitySpec
    quantity::Type
    bc::Type
    sweep::NamedTuple
    final_atol::Float64
end

"""
    LimitEdge

One asymptotic model→model limit — see [`@limits_to`](@ref).  `param` is the
driven source field, `approach` the strictly-monotone parameter sequence
(ordered toward the limit), `rate` optional human-readable convergence-rate
metadata, `quantities` the per-quantity convergence specs, and `mono_slack`
the relative slack of the error-shrinkage requirement (declare a looser value
for slowly/non-uniformly converging limits, e.g. logarithmic rates).
"""
struct LimitEdge
    name::Symbol
    source::Type
    target::Type
    param::Symbol
    approach::Vector{Float64}
    regime::String
    rate::Union{String,Nothing}
    quantities::Vector{LimitQuantitySpec}
    finite_N::Int
    mono_slack::Float64
    notes::String
    references::Vector{String}
end

"""
    LIMIT_EDGES :: Vector{LimitEdge}

The asymptotic-limit store, populated at include-time by [`limits_to!`](@ref)
/ [`@limits_to`](@ref).  Query with [`limits_from`](@ref) /
[`limits_into`](@ref); the generated convergence-sequence checks are the
`:limit` kind of [`generated_checks`](@ref).
"""
const LIMIT_EDGES = LimitEdge[]

"""
    limits_to!(name, source_T, target_T; param, approach, regime,
               quantities, rate=nothing, finite_N=8, mono_slack=0.1,
               notes="", references=String[])

Record an asymptotic limit edge.  `approach` must be strictly monotone (its
direction encodes the side the limit is taken from); `quantities` is an
iterable of NamedTuples `(quantity=Q, bc=BC, final_atol=ε[, sweep=(…)])`;
`mono_slack` is the relative tolerance of the error-shrinkage requirement
(each error may exceed its predecessor by at most this fraction).
"""
function limits_to!(
    name::Symbol,
    source_T::Type,
    target_T::Type;
    param::Symbol,
    approach::AbstractVector{<:Real},
    regime::AbstractString,
    quantities,
    rate::Union{AbstractString,Nothing}=nothing,
    finite_N::Int=8,
    mono_slack::Real=0.1,
    notes::AbstractString="",
    references::AbstractVector{<:AbstractString}=String[],
)
    mono_slack ≥ 0 || throw(ArgumentError("limits_to!: mono_slack must be ≥ 0"))
    any(l -> l.name === name, LIMIT_EDGES) &&
        throw(ArgumentError("limits_to!: :$(name) already declared"))
    length(approach) ≥ 2 ||
        throw(ArgumentError("limits_to!: approach needs ≥ 2 points to witness convergence"))
    diffs = diff(collect(Float64, approach))
    (all(>(0), diffs) || all(<(0), diffs)) ||
        throw(ArgumentError("limits_to!: approach must be strictly monotone"))
    isempty(strip(regime)) &&
        throw(ArgumentError("limits_to!: regime must describe the limit (e.g. \"Δ → 1⁺\")"))
    specs = LimitQuantitySpec[]
    for spec in quantities
        Q = spec.quantity
        Q isa Type && Q <: AbstractQuantity ||
            throw(ArgumentError("limits_to!: quantity $(Q) is not a quantity type"))
        _quantity_instance(Q)
        BC = spec.bc
        BC isa Type && BC <: BoundaryCondition ||
            throw(ArgumentError("limits_to!: bc $(BC) is not a boundary-condition type"))
        spec.final_atol ≥ 0 || throw(
            ArgumentError("limits_to!: final_atol must be ≥ 0; got $(spec.final_atol)")
        )
        push!(
            specs,
            LimitQuantitySpec(
                Q,
                BC,
                haskey(spec, :sweep) ? spec.sweep : NamedTuple(),
                Float64(spec.final_atol),
            ),
        )
    end
    isempty(specs) &&
        throw(ArgumentError("limits_to!: the quantity list must be non-empty"))
    push!(
        LIMIT_EDGES,
        LimitEdge(
            name,
            source_T,
            target_T,
            param,
            collect(Float64, approach),
            String(regime),
            rate === nothing ? nothing : String(rate),
            specs,
            finite_N,
            Float64(mono_slack),
            String(notes),
            String[r for r in references],
        ),
    )
    return nothing
end

"""
    @limits_to :name Source Target param=:Δ approach=[…] regime="…" quantities=[…] …

Macro sugar around [`limits_to!`](@ref): `:name` and the `Source`/`Target`
model types are positional; the remaining `key=value` pairs are forwarded as
keyword arguments.  See `src/limits_registry.jl` for the declared catalog.
"""
macro limits_to(name, source_T, target_T, kwargs...)
    return _forward_kw_macro(limits_to!, :limits_to, (name, source_T, target_T), kwargs)
end

"""
    limits_from(model) -> Vector{NamedTuple}

The asymptotic limits of `model` (as the driven source):
`(name, target, param, regime, rate, references)` rows.
"""
function limits_from(model)
    m_T = _as_type(model)
    return [
        (
            name=l.name,
            target=l.target,
            param=l.param,
            regime=l.regime,
            rate=l.rate,
            references=l.references,
        ) for l in LIMIT_EDGES if l.source === m_T
    ]
end

"""
    limits_into(model) -> Vector{NamedTuple}

The models that asymptotically approach `model` — the inverse of
[`limits_from`](@ref): `(name, source, param, regime, rate, references)` rows.
"""
function limits_into(model)
    m_T = _as_type(model)
    return [
        (
            name=l.name,
            source=l.source,
            param=l.param,
            regime=l.regime,
            rate=l.rate,
            references=l.references,
        ) for l in LIMIT_EDGES if l.target === m_T
    ]
end

# ──────────────────────────────────────────────────────────────────────
# C13 — static limit coherence (no fetch execution)
# ──────────────────────────────────────────────────────────────────────

"""
    check_limit_edges() -> Vector{CoherenceFinding}

C13: static sanity of every limit edge:

  * the driven `param` is a field of the source model (`:error` — checked by
    reconstructing the default instance at the first approach point);
  * each quantity spec has canonical, independent (non-delegating) registry
    rows on BOTH endpoints at its `bc` — missing or delegation-backed rows
    are a `:gap` (the convergence check would be absent or circular,
    #701's rule).
"""
function check_limit_edges()
    out = CoherenceFinding[]
    for l in LIMIT_EDGES
        try
            _with_param(l.source(), l.param, first(l.approach))
        catch err
            push!(
                out,
                CoherenceFinding(
                    :limit_edge,
                    :error,
                    "limits_to :$(l.name): cannot drive $(_kgshort(l.source)).$(l.param)" *
                    " ($(typeof(err)))",
                ),
            )
        end
        for spec in l.quantities
            _check_endpoint_rows!(
                out,
                l.source,
                l.target,
                spec.quantity,
                spec.bc,
                :limit_edge,
                "limits_to :$(l.name)",
            )
        end
    end
    return out
end

# ──────────────────────────────────────────────────────────────────────
# Generator — the :limit kind of generated_checks()
# ──────────────────────────────────────────────────────────────────────

# Absolute floor of the error-shrinkage comparison: two consecutive errors
# both at machine precision must not fail monotonicity on round-off alone.
const _MONO_ATOL_FLOOR = 1e-14

function limit_checks()
    out = GeneratedCheck[]
    for l in LIMIT_EDGES, spec in l.quantities
        _both_endpoints_independent(l.source, l.target, spec.quantity, spec.bc) || continue
        for point in _sweep_points(spec.sweep)
            id = string(
                "limit/",
                l.name,
                "/",
                _kgshort(spec.quantity),
                "/",
                _kgshort(spec.bc),
                _point_suffix(point),
            )
            runner = function ()
                bc = _bc_instance(spec.bc; finite_N=l.finite_N)
                q = _quantity_instance(spec.quantity)
                t = fetch(l.target(), q, bc; point...)
                errs = Float64[]
                for x in l.approach
                    v = fetch(_with_param(l.source(), l.param, x), q, bc; point...)
                    push!(errs, abs(Float64(v) - Float64(t)))
                end
                # weak monotonicity along the sequence: the declared
                # mono_slack absorbs the limit's own convergence wiggle, the
                # absolute floor the machine-precision noise of an
                # already-converged tail …
                shrinking = all(
                    errs[i + 1] ≤ (1 + l.mono_slack) * errs[i] + _MONO_ATOL_FLOOR for
                    i in 1:(length(errs) - 1)
                )
                # … and the terminal error must land below the declared atol.
                detail = string(
                    l.regime,
                    "; |source−target| along ",
                    l.param,
                    "=",
                    l.approach,
                    " → ",
                    errs,
                )
                shrinking || return CheckOutcome(
                    :fail,
                    errs[end],
                    0.0,
                    errs[end],
                    errs[end] / max(errs[end], eps()),
                    "error sequence not shrinking: " * detail,
                )
                return _outcome(
                    errs[end], 0.0; rtol=0.0, atol=spec.final_atol, detail=detail
                )
            end
            push!(
                out,
                GeneratedCheck(
                    :limit,
                    id,
                    string(
                        "limit :",
                        l.name,
                        " (",
                        l.regime,
                        "): ",
                        _kgshort(spec.quantity),
                        " convergence ",
                        _kgshort(l.source),
                        " → ",
                        _kgshort(l.target),
                    ),
                    runner,
                ),
            )
        end
    end
    return out
end

register_check_generator!(:limit, limit_checks)
register_edge_store!(:limits_to, LIMIT_EDGES; location_of=l -> "limits_to :$(l.name)")
