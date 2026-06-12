# core/duality.jl — parameter-mapped duality edges (#699, over the #697
# kernel).
#
# A duality edge declares a bijective, parameter-mapped equivalence between
# two model families: `target(param_map(m))` computes the same physics as the
# source instance `m`, possibly up to a per-quantity renormalisation
# (`value_map` — operator conventions, additive constants).  Distinct from
# `@reduces` (exact delegation at a point in parameter space): a duality holds
# on the whole parameter manifold and maps it nontrivially.
#
# WHY this store exists — the SU(2)₁ boundary-entropy incident: a wrong
# registered value survived review because its test parroted the source's own
# claimed numbers (circular verification).  The bug was caught by a
# hand-derived cross-implementation argument.  A declared duality edge
# *generates* that argument for every quantity both sides implement: fetch on
# side A, map parameters, fetch on side B, compare.  These are independent-
# implementation cross-checks — the strongest verification class available
# without external data — and they are precisely what delegation-backed rows
# cannot provide, so the generator enforces the kernel's independence filter.
#
# The per-quantity comparison spec is an explicit ALLOWLIST: a duality maps
# the spectrum, but not every observable is self-dual (the TFIM
# Kramers–Wannier map exchanges order and disorder parameters, so
# MagnetizationX is deliberately NOT on the edge's quantity list while the
# thermal potentials and the gap are).  Explicit beats inferred here — a
# wrongly-inferred quantity comparison would manufacture false :error
# findings out of correct physics.

"""
    DualityQuantitySpec

One quantity both sides of a [`Duality`](@ref) must agree on: compared at
boundary condition `bc`, on the fetch-kwargs grid `sweep`, after mapping the
target-side value through `value_map(value, source_instance)` (identity by
default; carries operator renormalisations and additive constants, e.g. the
Jordan–Wigner `−h` energy-density offset between TFIM and the Kitaev wire).
"""
struct DualityQuantitySpec
    quantity::Type
    bc::Type
    sweep::NamedTuple
    value_map::Function
end

"""
    Duality

One parameter-mapped model↔model equivalence — see [`@dual`](@ref).
`param_map` sends a source instance to the equivalent target instance;
`examples` are the source instances the generated cross-checks run at (chosen
off any self-dual locus so the map is exercised nontrivially); `involution`
asserts `param_map ∘ param_map ≈ id` (checked statically on the examples).
"""
struct Duality
    name::Symbol
    source::Type
    target::Type
    param_map::Function
    kind::Symbol
    involution::Bool
    examples::Vector{Any}
    quantities::Vector{DualityQuantitySpec}
    finite_N::Int
    rtol::Float64
    atol::Float64
    regime::String
    notes::String
    references::Vector{String}
end

"""
    DUALITIES :: Vector{Duality}

The duality-edge store, populated at include-time by [`dual!`](@ref) /
[`@dual`](@ref).  Query with [`dualities`](@ref); the generated
cross-implementation checks are the `:dual` kind of
[`generated_checks`](@ref).
"""
const DUALITIES = Duality[]

"""
    dual!(name, source_T, target_T; param_map, kind, quantities,
          examples, involution=false, finite_N=8, rtol=1e-8, atol=1e-10,
          regime="", notes="", references=String[])

Record a duality edge.  `quantities` is an iterable of NamedTuples
`(quantity=Q, bc=BC[, sweep=(…)][, value_map=f])` — the explicit allowlist of
observables the duality maps; `examples` must be non-empty source instances.
"""
function dual!(
    name::Symbol,
    source_T::Type,
    target_T::Type;
    param_map::Function,
    kind::Symbol,
    quantities,
    examples::AbstractVector,
    involution::Bool=false,
    finite_N::Int=8,
    rtol::Real=1e-8,
    atol::Real=1e-10,
    regime::AbstractString="",
    notes::AbstractString="",
    references::AbstractVector{<:AbstractString}=String[],
)
    any(d -> d.name === name, DUALITIES) &&
        throw(ArgumentError("dual!: :$(name) already declared"))
    isempty(examples) && throw(
        ArgumentError(
            "dual!: a duality needs example source instances for its generated " *
            "cross-checks (pick points off any self-dual locus)",
        ),
    )
    # An involution maps a parameter manifold to ITSELF; with distinct
    # endpoint families param_map² is not even type-stable, so the C12
    # involution probe would feed param_map a foreign type.
    involution &&
        source_T !== target_T &&
        throw(
            ArgumentError(
                "dual!: involution=true requires source === target (a self-duality); " *
                "got $(source_T) ↔ $(target_T)",
            ),
        )
    specs = DualityQuantitySpec[]
    for spec in quantities
        Q = spec.quantity
        Q isa Type && Q <: AbstractQuantity ||
            throw(ArgumentError("dual!: quantity $(Q) is not a quantity type"))
        _quantity_instance(Q)
        BC = spec.bc
        BC isa Type && BC <: BoundaryCondition ||
            throw(ArgumentError("dual!: bc $(BC) is not a boundary-condition type"))
        push!(
            specs,
            DualityQuantitySpec(
                Q,
                BC,
                haskey(spec, :sweep) ? spec.sweep : NamedTuple(),
                haskey(spec, :value_map) ? spec.value_map : (v, _) -> v,
            ),
        )
    end
    isempty(specs) &&
        throw(ArgumentError("dual!: the quantity allowlist must be non-empty"))
    push!(
        DUALITIES,
        Duality(
            name,
            source_T,
            target_T,
            param_map,
            kind,
            involution,
            Any[m for m in examples],
            specs,
            finite_N,
            Float64(rtol),
            Float64(atol),
            String(regime),
            String(notes),
            String[r for r in references],
        ),
    )
    return nothing
end

"""
    @dual :name Source Target param_map=… kind=… quantities=[…] examples=[…] …

Macro sugar around [`dual!`](@ref): `:name` and the `Source`/`Target` model
types are positional; the remaining `key=value` pairs are forwarded as
keyword arguments.  See `src/duality_registry.jl` for the declared catalog.
"""
macro dual(name, source_T, target_T, kwargs...)
    return _forward_kw_macro(dual!, :dual, (name, source_T, target_T), kwargs)
end

"""
    dualities(model) -> Vector{NamedTuple}

The duality edges touching `model` (as source or target):
`(name, source, target, kind, regime, references)` rows.
"""
function dualities(model)
    m_T = _as_type(model)
    return [
        (
            name=d.name,
            source=d.source,
            target=d.target,
            kind=d.kind,
            regime=d.regime,
            references=d.references,
        ) for d in DUALITIES if d.source === m_T || d.target === m_T
    ]
end

# ──────────────────────────────────────────────────────────────────────
# C12 — static duality coherence (param_map sanity; no fetch execution)
# ──────────────────────────────────────────────────────────────────────

# Field-wise approximate equality of two model instances (for the involution
# check); models are plain concrete structs of Real fields.
function _instances_approx(a, b; rtol=1e-12)
    typeof(a) === typeof(b) || return false
    return all(
        isapprox(getfield(a, f), getfield(b, f); rtol=rtol) for f in propertynames(a)
    )
end

"""
    check_duality_maps() -> Vector{CoherenceFinding}

C12: static sanity of every duality edge, evaluated on its registered
examples only (the C8 pattern — predicates run, `fetch` does not):

  * each example is a `source` instance and `param_map(example)` is a
    `target` instance (`:error` otherwise — the map is malformed);
  * an `involution=true` edge satisfies `param_map(param_map(x)) ≈ x` on
    every example (`:error`);
  * each quantity spec has canonical, independent (non-delegating) registry
    rows on BOTH endpoints at its `bc` — a missing or delegation-backed row
    is a `:gap`: the cross-check the edge promises cannot be generated yet.
"""
function check_duality_maps()
    out = CoherenceFinding[]
    for d in DUALITIES
        for (i, ex) in enumerate(d.examples)
            if !(ex isa d.source)
                push!(
                    out,
                    CoherenceFinding(
                        :duality_map,
                        :error,
                        "dual :$(d.name) example #$(i) is a $(typeof(ex)), not a " *
                        "$(_kgshort(d.source)) instance",
                    ),
                )
                continue
            end
            mapped = try
                d.param_map(ex)
            catch err
                push!(
                    out,
                    CoherenceFinding(
                        :duality_map,
                        :error,
                        "dual :$(d.name) param_map threw on example #$(i) " *
                        "($(typeof(err)))",
                    ),
                )
                continue
            end
            if !(mapped isa d.target)
                push!(
                    out,
                    CoherenceFinding(
                        :duality_map,
                        :error,
                        "dual :$(d.name) param_map(example #$(i)) is a " *
                        "$(typeof(mapped)), not a $(_kgshort(d.target)) instance",
                    ),
                )
                continue   # a malformed image must not feed the involution probe
            end
            if d.involution   # dual! guarantees source === target here
                back = try
                    d.param_map(mapped)
                catch err
                    push!(
                        out,
                        CoherenceFinding(
                            :duality_map,
                            :error,
                            "dual :$(d.name) param_map threw on its own image of " *
                            "example #$(i) ($(typeof(err))) — not an involution",
                        ),
                    )
                    continue
                end
                _instances_approx(back, ex) || push!(
                    out,
                    CoherenceFinding(
                        :duality_map,
                        :error,
                        "dual :$(d.name) declares involution=true but " *
                        "param_map²(example #$(i)) ≠ example",
                    ),
                )
            end
        end
        for spec in d.quantities
            for (side, m_T) in ((:source, d.source), (:target, d.target))
                row = _canonical_row(m_T, spec.quantity, spec.bc)
                if row === nothing
                    push!(
                        out,
                        CoherenceFinding(
                            :duality_map,
                            :gap,
                            "dual :$(d.name) lists $(_kgshort(spec.quantity)) at " *
                            "$(_kgshort(spec.bc)) but the $(side) side " *
                            "($(_kgshort(m_T))) has no canonical row — the promised " *
                            "cross-check cannot be generated",
                        ),
                    )
                elseif !_is_independent_row(row)
                    push!(
                        out,
                        CoherenceFinding(
                            :duality_map,
                            :gap,
                            "dual :$(d.name): the $(side) row for " *
                            "$(_kgshort(spec.quantity)) at $(_kgshort(spec.bc)) is " *
                            "delegation-backed — the cross-check would be circular " *
                            "and is not generated",
                        ),
                    )
                end
            end
        end
    end
    return out
end

# ──────────────────────────────────────────────────────────────────────
# Generator — the :dual kind of generated_checks()
# ──────────────────────────────────────────────────────────────────────

# A duality cross-check is only emitted when both endpoint rows exist and are
# independent (non-delegating) — the structural condition that makes it a
# genuine two-implementation comparison.  Missing/delegating rows are visible
# as C12 :gap findings, not silently skipped here.
function _dual_rows_independent(d::Duality, spec::DualityQuantitySpec)
    for m_T in (d.source, d.target)
        row = _canonical_row(m_T, spec.quantity, spec.bc)
        (row !== nothing && _is_independent_row(row)) || return false
    end
    return true
end

function dual_checks()
    out = GeneratedCheck[]
    for d in DUALITIES, spec in d.quantities
        _dual_rows_independent(d, spec) || continue
        for (i, ex) in enumerate(d.examples), point in _sweep_points(spec.sweep)
            id = string(
                "dual/",
                d.name,
                "/",
                _kgshort(spec.quantity),
                "/",
                _kgshort(spec.bc),
                "/ex",
                i,
                isempty(keys(point)) ? "" : "/" * _point_id(point),
            )
            runner = function ()
                bc = _bc_instance(spec.bc; finite_N=d.finite_N)
                q = _quantity_instance(spec.quantity)
                a = fetch(ex, q, bc; point...)
                b = fetch(d.param_map(ex), q, bc; point...)
                return _outcome(a, spec.value_map(b, ex); rtol=d.rtol, atol=d.atol)
            end
            push!(
                out,
                GeneratedCheck(
                    :dual,
                    id,
                    string(
                        "dual :",
                        d.name,
                        " (",
                        d.kind,
                        "): ",
                        _kgshort(spec.quantity),
                        " of ",
                        ex,
                        " vs ",
                        _kgshort(d.target),
                        " image",
                    ),
                    runner,
                ),
            )
        end
    end
    return out
end

register_check_generator!(:dual, dual_checks)
register_edge_store!(:dual, DUALITIES; location_of=d -> "dual :$(d.name)")
