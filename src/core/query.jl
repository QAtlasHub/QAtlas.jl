# core/query.jl — availability search over the atlas (the "use" face)
#
# `implementation_status()` returns every registry row, but there was no way to ask "does the atlas
# HAVE X?" and get a yes/no. `search` answers that — by model, quantity, boundary condition, and/or
# a coarse `regime` facet (ground-state / finite-temperature / dynamics / universality) — and
# `search_jsonl` streams the answer as JSONL: a summary line `{available, query, count, verified}`
# then one object per hit. A consumer (an LLM, or a shell pipeline) can ask "is there
# finite-temperature TFIM data?" and branch on the first line alone.
#
# Each hit carries its provenance (status, reliability, references) and whether the model is tied
# into QAtlas's executable cross-checks (`cross_checked`) — so "available" is qualified by "and how
# trustworthy", not a bare bool. The live pass/fail verdict (running the numerics) is a deeper,
# on-demand query, deliberately NOT run here so search stays fast. JSON is hand-rolled (the fields
# are clean identifiers / bibkeys) to keep the core dependency-free.

"""
    REGIMES

Coarse capability facets, each a predicate over a registry row. Declarative and extensible — add a
regime by adding a predicate. `:dynamics` is intentionally sparse today (QAtlas has velocities but
no spectral / real-time data yet); a search for it honestly reports the gap.
"""
const REGIMES = (
    ground_state=impl -> (
        impl.quantity <: AbstractGap ||
        impl.quantity <: AbstractEntanglementMeasure ||
        impl.quantity === GroundStateEnergyDensity
    ),
    finite_temperature=impl -> impl.quantity <: AbstractThermalPotential,
    universality=impl -> impl.status === :universal,
    dynamics=impl -> impl.quantity <: AbstractVelocity,
)

"""
    Hit

One available (model, quantity, bc) hub, with its provenance: `method` (how the value is obtained),
`status`, `reliability`, `references`, and `cross_checked` — `true` when the model participates in
an executable cross-check (a duality, limit, or LSM-symmetry edge), the QAtlas-specific trust signal
(model-level in this version).
"""
struct Hit
    model::String
    quantity::String
    bc::String
    method::Symbol        # how the value is obtained (:bdg, :bethe_ansatz, :analytic, …)
    status::Symbol        # :exact / :bound / :approx / :universal
    reliability::Symbol   # :high / :medium / :low
    cross_checked::Bool
    references::Vector{String}
end

"""
    QueryResult

The result of a [`search`](@ref): the echoed `query`, an `available` flag, and the matching `hits`.
"""
struct QueryResult
    query::NamedTuple
    available::Bool
    hits::Vector{Hit}
end

_label(@nospecialize T) = replace(string(T), "QAtlas." => "")
_norm(s) = lowercase(replace(string(s), "_" => ""))

# Facet matching: `nothing` matches anything; a Symbol/String matches by case- and
# underscore-insensitive substring on the type name (so `:specific_heat` finds `SpecificHeat`,
# `:ising` finds every Ising* model); a Type matches as a supertype (so `quantity=AbstractGap`
# catches the whole family).
_match(::Nothing, @nospecialize(T)) = true
_match(want::Symbol, @nospecialize(T)) = occursin(_norm(want), _norm(_label(T)))
_match(want::AbstractString, @nospecialize(T)) = occursin(_norm(want), _norm(_label(T)))
_match(want::Type, @nospecialize(T)) = T <: want

# A reference facet matches if any of a row's bibkeys contains it (case/underscore-insensitive).
_ref_match(want, refs) = any(r -> occursin(_norm(want), _norm(r)), refs)

# Models tied into an executable cross-check (duality / limit / LSM symmetry). Model-level for now.
function _cross_checked_models()
    s = Set{Type}()
    for d in DUALITIES
        push!(s, d.source)
        push!(s, d.target)
    end
    for l in LIMIT_EDGES
        push!(s, l.source)
        push!(s, l.target)
    end
    for p in SYMMETRY_PROFILES
        push!(s, p.model)
    end
    return s
end

"""
    search(; model, quantity, bc, regime, status, reliability, method, reference, cross_checked) -> QueryResult

Does the atlas have data matching the given facets? Any facet left `nothing` is unconstrained; all
given facets are AND-combined.

- `model`/`quantity`/`bc` — a Type (exact, or a family like `AbstractGap`) or a Symbol/String
  (fuzzy: case- and underscore-insensitive substring on the name).
- `regime` — one of `keys(REGIMES)` (`:ground_state`, `:finite_temperature`, `:dynamics`,
  `:universality`).
- `status` (`:exact`/`:bound`/`:approx`/`:universal`) and `reliability` (`:high`/`:medium`/`:low`)
  — exact Symbol match.
- `method` — a Symbol/String, fuzzy (so `:bethe` finds `:bethe_ansatz`).
- `reference` — a bibkey substring, fuzzy (so `"Pfeuty"` finds `"Pfeuty1970"`).
- `cross_checked` — `true`/`false`, the trust filter.

Returns a [`QueryResult`](@ref) — see [`search_jsonl`](@ref) for JSONL and [`available`](@ref) for
just the boolean.
"""
function search(;
    model=nothing,
    quantity=nothing,
    bc=nothing,
    regime=nothing,
    status=nothing,
    reliability=nothing,
    method=nothing,
    reference=nothing,
    cross_checked=nothing,
)
    regime === nothing ||
        haskey(REGIMES, regime) ||
        throw(ArgumentError("unknown regime $(repr(regime)); known: $(keys(REGIMES))"))
    xchecked = _cross_checked_models()
    hits = Hit[]
    for impl in REGISTRY
        _match(model, impl.model) || continue
        _match(quantity, impl.quantity) || continue
        _match(bc, impl.bc) || continue
        regime === nothing || REGIMES[regime](impl) || continue
        status === nothing || impl.status === status || continue
        reliability === nothing || impl.reliability === reliability || continue
        method === nothing || _match(method, impl.method) || continue
        reference === nothing || _ref_match(reference, impl.references) || continue
        xc = impl.model in xchecked
        cross_checked === nothing || xc === cross_checked || continue
        push!(
            hits,
            Hit(
                _label(impl.model),
                _label(impl.quantity),
                _label(impl.bc),
                impl.method,
                impl.status,
                impl.reliability,
                xc,
                copy(impl.references),
            ),
        )
    end
    sort!(hits; by=h -> (h.model, h.quantity, h.bc))
    return QueryResult(
        (;
            model,
            quantity,
            bc,
            regime,
            status,
            reliability,
            method,
            reference,
            cross_checked,
        ),
        !isempty(hits),
        hits,
    )
end

"""
    available(; model=nothing, quantity=nothing, bc=nothing, regime=nothing) -> Bool

The bare yes/no: does the atlas have anything matching the facets? Convenience over [`search`](@ref).
"""
available(; kwargs...) = search(; kwargs...).available

# ---- edge / relation search: a model's graph neighborhood ----

"""
    Relation

One edge in a model's graph neighborhood: `kind`
(`:dual`/`:limit_to`/`:limit_from`/`:reduces_to`/`:reduces_from`/`:realizes`/`:symmetry`), the
`from`/`to` endpoints (model names; a universality class for `:realizes`; `""` for `:symmetry`), a
`detail` summary, and `references`.
"""
struct Relation
    kind::Symbol
    from::String
    to::String
    detail::String
    references::Vector{String}
end

# The registered model Types matching a facet (a Type, or a fuzzy Symbol/String).
function _model_types(facet)
    seen = Type[]
    for impl in REGISTRY
        impl.model in seen && continue
        _match(facet, impl.model) || continue
        push!(seen, impl.model)
    end
    return seen
end

"""
    relations(model) -> Vector{Relation}

The graph neighborhood of `model` (a Type, or a fuzzy Symbol/String facet): every duality, limit,
reduction, realization, and symmetry edge it participates in, one `Relation` each. See
[`relations_jsonl`](@ref) for JSONL output.
"""
function relations(model)
    rels = Relation[]
    for M in _model_types(model)
        nm = _label(M)
        for d in DUALITIES
            (d.source === M || d.target === M) || continue
            partner = d.source === M ? d.target : d.source
            push!(
                rels,
                Relation(
                    :dual,
                    nm,
                    _label(partner),
                    "kind=$(d.kind); $(d.regime)",
                    copy(d.references),
                ),
            )
        end
        for l in LIMIT_EDGES
            l.source === M && push!(
                rels,
                Relation(
                    :limit_to,
                    nm,
                    _label(l.target),
                    "param=$(l.param); $(l.regime)",
                    copy(l.references),
                ),
            )
            l.target === M && push!(
                rels,
                Relation(
                    :limit_from,
                    _label(l.source),
                    nm,
                    "param=$(l.param); $(l.regime)",
                    copy(l.references),
                ),
            )
        end
        for r in REDUCES
            r.source === M && push!(
                rels,
                Relation(:reduces_to, nm, _label(r.target), r.regime, copy(r.references)),
            )
            r.target === M && push!(
                rels,
                Relation(:reduces_from, _label(r.source), nm, r.regime, copy(r.references)),
            )
        end
        for r in REALIZES
            r.model === M && push!(
                rels,
                Relation(:realizes, nm, string(r.class), r.regime, copy(r.references)),
            )
        end
        for p in SYMMETRY_PROFILES
            p.model === M && push!(
                rels,
                Relation(
                    :symmetry,
                    nm,
                    "",
                    "internal=$(p.internal), translation=$(p.translation)",
                    copy(p.references),
                ),
            )
        end
    end
    sort!(rels; by=r -> (string(r.kind), r.from, r.to))
    return rels
end

# ---- gap / absence search: a model's coverage holes ----

"""
    Gap

A coverage hole — something the atlas does NOT have. `kind` is `:regime` (a `REGIMES` capability the
`model` does not cover); `subject` names what is missing. Atlas-wide structural gaps are separate —
see `coherence_gaps()`.
"""
struct Gap
    kind::Symbol
    subject::String
    model::String
end

"""
    gaps(model) -> Vector{Gap}

What the atlas does NOT have for `model` (a Type, or a fuzzy Symbol/String facet): for each `regime`
in `REGIMES` the model does not cover, a `Gap(:regime, regime, model)`. The grounded negative of
[`search`](@ref) — no guessed "expected set", just the known capability facets the model lacks. See
[`gaps_jsonl`](@ref) for JSONL.
"""
function gaps(model)
    out = Gap[]
    for M in _model_types(model)
        nm = _label(M)
        for regime in keys(REGIMES)
            isempty(search(; model=M, regime=regime).hits) &&
                push!(out, Gap(:regime, string(regime), nm))
        end
    end
    sort!(out; by=g -> (g.model, g.subject))
    return out
end

# ---- JSONL serialization (hand-rolled; fields are clean identifiers / bibkeys) ----

function _json_escape(s::AbstractString)
    buf = IOBuffer()
    for c in s
        if c == '"'
            print(buf, "\\\"")
        elseif c == '\\'
            print(buf, "\\\\")
        elseif c == '\n'
            print(buf, "\\n")
        elseif c == '\t'
            print(buf, "\\t")
        elseif c == '\r'
            print(buf, "\\r")
        elseif c < ' '
            print(buf, "\\u", lpad(string(UInt16(c); base=16), 4, '0'))
        else
            print(buf, c)
        end
    end
    return String(take!(buf))
end
_jstr(s) = string('"', _json_escape(string(s)), '"')
_jarr(v) = string('[', join((_jstr(x) for x in v), ','), ']')

function _json_hit(h::Hit)
    return string(
        "{\"model\":",
        _jstr(h.model),
        ",\"quantity\":",
        _jstr(h.quantity),
        ",\"bc\":",
        _jstr(h.bc),
        ",\"method\":",
        _jstr(h.method),
        ",\"status\":",
        _jstr(h.status),
        ",\"reliability\":",
        _jstr(h.reliability),
        ",\"cross_checked\":",
        h.cross_checked,
        ",\"references\":",
        _jarr(h.references),
        "}",
    )
end

function _json_query(q::NamedTuple)
    parts = String[]
    for (k, v) in pairs(q)
        v === nothing && continue
        val = v isa Bool ? string(v) : _jstr(v isa Type ? _label(v) : v)
        push!(parts, string(_jstr(k), ":", val))
    end
    return string('{', join(parts, ','), '}')
end

"""
    search_jsonl([io=stdout]; model, quantity, bc, regime) -> nothing

Stream a [`search`](@ref) as JSONL. The first line is a summary —
`{"available":…,"query":…,"count":N,"verified":K}` — so a consumer can branch on it without reading
the rest; each following line is one hit object (model, quantity, bc, status, reliability,
cross_checked, references). `verified` counts hits that are cross-checked or exact/universal.
"""
function search_jsonl(io::IO=stdout; kwargs...)
    r = search(; kwargs...)
    nverified = count(
        h -> h.cross_checked || h.status === :exact || h.status === :universal, r.hits
    )
    print(
        io,
        "{\"available\":",
        r.available,
        ",\"query\":",
        _json_query(r.query),
        ",\"count\":",
        length(r.hits),
        ",\"verified\":",
        nverified,
        "}\n",
    )
    for h in r.hits
        print(io, _json_hit(h), "\n")
    end
    return nothing
end

function _json_relation(r::Relation)
    return string(
        "{\"kind\":",
        _jstr(r.kind),
        ",\"from\":",
        _jstr(r.from),
        ",\"to\":",
        _jstr(r.to),
        ",\"detail\":",
        _jstr(r.detail),
        ",\"references\":",
        _jarr(r.references),
        "}",
    )
end

"""
    relations_jsonl([io=stdout], model) -> nothing

Stream [`relations`](@ref) as JSONL: a `{available, query, count}` summary line then one Relation
object per line.
"""
function relations_jsonl(io::IO, model)
    rels = relations(model)
    print(
        io,
        "{\"available\":",
        !isempty(rels),
        ",\"query\":",
        _json_query((; model)),
        ",\"count\":",
        length(rels),
        "}\n",
    )
    for r in rels
        print(io, _json_relation(r), "\n")
    end
    return nothing
end
relations_jsonl(model) = relations_jsonl(stdout, model)

function _json_gap(g::Gap)
    return string(
        "{\"kind\":",
        _jstr(g.kind),
        ",\"subject\":",
        _jstr(g.subject),
        ",\"model\":",
        _jstr(g.model),
        "}",
    )
end

"""
    gaps_jsonl([io=stdout], model) -> nothing

Stream [`gaps`](@ref) as JSONL: a `{has_gaps, query, count}` summary line then one Gap object per
line.
"""
function gaps_jsonl(io::IO, model)
    gs = gaps(model)
    print(
        io,
        "{\"has_gaps\":",
        !isempty(gs),
        ",\"query\":",
        _json_query((; model)),
        ",\"count\":",
        length(gs),
        "}\n",
    )
    for g in gs
        print(io, _json_gap(g), "\n")
    end
    return nothing
end
gaps_jsonl(model) = gaps_jsonl(stdout, model)
