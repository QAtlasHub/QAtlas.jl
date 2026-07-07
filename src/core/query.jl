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

Coarse capability facets, kept as a convenience over the canonical query axes
(`family` × `thermal` × `dynamical` × `status`). Each predicate is now **derived
from those axes** rather than a hand-picked quantity list, so a whole quantity
class can no longer silently fall through: the thermal regimes are keyed on the
`thermal` axis (so `:ground_state` surfaces T=0-accessible correlations,
magnetization, … — not only gaps/entanglement), and each predicate is a superset
of its previous hand-list (no regression). `:dynamics` stays sparse (velocities are
`:transport`; genuine spectral/real-time data would be `:dynamic`); a search for it
honestly reports the gap. For a rigorous query, prefer the axes directly.
"""
const REGIMES = (
    ground_state=impl -> (
        impl.thermal === :zero ||
        impl.thermal === :both ||
        impl.quantity <: AbstractGap ||
        impl.quantity <: AbstractEntanglementMeasure ||
        impl.quantity === GroundStateEnergyDensity
    ),
    finite_temperature=impl -> (
        impl.thermal === :finite ||
        impl.thermal === :both ||
        impl.quantity <: AbstractThermalPotential
    ),
    universality=impl -> impl.status === :universal,
    dynamics=impl -> (
        impl.dynamical === :dynamic ||
        impl.dynamical === :transport ||
        impl.quantity <: AbstractVelocity
    ),
)

"""
    Hit

One available (model, quantity, bc) hub, with its `family` (the quantity's
super-family — `:correlation`, `:magnetization`, …; total, so it never silently
drops a class), its provenance, AND how to obtain the value. Provenance: `method`,
`status`, `reliability`, `references`, `cross_checked`
(`true` when the model is in an executable cross-check — a duality, limit, or
LSM-symmetry edge — the QAtlas trust signal, model-level here). Actionability:

- `params` — the fetch's declared keyword arguments (`i`, `j`, `beta`, `N`, `k`,
  couplings, …): the call signature, so a hit says not just *available* but *how to
  call it*. A fetch that declares `N`/`finite_N` is evaluable at any finite size
  (the OBC/PBC hubs are finite-N closed forms; `Infinite` is the thermodynamic
  limit) — the atlas has no separate `N` dimension, so this is where size shows up.
- `notes` / `valid_domain` — the closed form's remark and where it holds (from the
  `Implementation` row), so a consumer can match its numerics to the functional
  form / validity range. `valid_domain` is `""` when unset.
"""
struct Hit
    model::String
    quantity::String
    bc::String
    family::Symbol        # quantity super-family (:correlation, :magnetization, …); TOTAL
    method::Symbol        # how the value is obtained (:bdg, :bethe_ansatz, :analytic, …)
    status::Symbol        # :exact / :bound / :approx / :universal
    reliability::Symbol   # :high / :medium / :low
    thermal::Symbol       # orthogonal axis: :zero / :finite / :both / :unknown
    dynamical::Symbol     # orthogonal axis: :static / :transport / :dynamic / :unknown
    cross_checked::Bool
    references::Vector{String}
    params::Vector{String}  # fetch's declared kwargs — the call signature (how to obtain it)
    notes::String           # the Implementation row's note (functional form / caveats)
    valid_domain::String    # where the value holds (:approx rows); "" when unset
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

# An orthogonal-axis facet matches its stored value; `:both` (a hub spanning T=0 and T>0) matches
# either `thermal` query.
_axis_match(want, have) = want === have || have === :both

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

# The fetch's declared keyword arguments for a (model, quantity, bc) hub — the call
# signature a consumer needs (i/j sites, beta, N/finite_N, k, couplings…). Same
# introspection `axes.jl` uses; the `kwargs...` slurp is dropped. Best-effort: an
# empty vector if no method resolves.
function _fetch_params(M::Type, Q::Type, BC::Type)
    ps = String[]
    try
        for mth in methods(fetch, Tuple{M,Q,BC})
            for kw in Base.kwarg_decl(mth)
                kw === :kwargs && continue
                push!(ps, String(kw))
            end
        end
    catch
    end
    return sort(unique(ps))
end

"""
    search(; model, quantity, bc, regime, status, reliability, method, reference, cross_checked) -> QueryResult

Does the atlas have data matching the given facets? Any facet left `nothing` is unconstrained; all
given facets are AND-combined.

- `model`/`quantity`/`bc` — a Type (exact, or a family like `AbstractGap`) or a Symbol/String
  (fuzzy: case- and underscore-insensitive substring on the name).
- `family` — the quantity super-family (`:correlation`, `:structure_factor`, `:magnetization`,
  `:susceptibility`, `:thermodynamic`, `:gap`, `:entanglement`, `:velocity`, `:other`), exact
  Symbol match. TOTAL: every quantity has one, so no class silently drops (unlike `regime`).
- `regime` — one of `keys(REGIMES)` (`:ground_state`, `:finite_temperature`, `:dynamics`,
  `:universality`), a convenience derived from the axes; prefer `family`/`thermal`/`dynamical`.
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
    family=nothing,
    regime=nothing,
    status=nothing,
    reliability=nothing,
    method=nothing,
    reference=nothing,
    cross_checked=nothing,
    thermal=nothing,
    dynamical=nothing,
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
        family === nothing || quantity_family(impl.quantity) === family || continue
        regime === nothing || REGIMES[regime](impl) || continue
        status === nothing || impl.status === status || continue
        reliability === nothing || impl.reliability === reliability || continue
        method === nothing || _match(method, impl.method) || continue
        reference === nothing || _ref_match(reference, impl.references) || continue
        thermal === nothing || _axis_match(thermal, impl.thermal) || continue
        dynamical === nothing || _axis_match(dynamical, impl.dynamical) || continue
        xc = impl.model in xchecked
        cross_checked === nothing || xc === cross_checked || continue
        push!(
            hits,
            Hit(
                _label(impl.model),
                _label(impl.quantity),
                _label(impl.bc),
                quantity_family(impl.quantity),
                impl.method,
                impl.status,
                impl.reliability,
                impl.thermal,
                impl.dynamical,
                xc,
                copy(impl.references),
                _fetch_params(impl.model, impl.quantity, impl.bc),
                impl.notes,
                impl.valid_domain === nothing ? "" : impl.valid_domain,
            ),
        )
    end
    sort!(hits; by=h -> (h.model, h.quantity, h.bc))
    return QueryResult(
        (;
            model,
            quantity,
            bc,
            family,
            regime,
            status,
            reliability,
            method,
            reference,
            cross_checked,
            thermal,
            dynamical,
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

# ---- describe: the full per-model grounding record (the use-face payoff) ----

"""
    ModelRecord

Everything the atlas records about one model in a single grounding record: `summary` / `hamiltonian`
(from the `@about` card — `""` if uncarded, coverage is partial), the distinct `quantities`
available, the `relations` (graph neighborhood), and aggregated `references`. The structural content
(quantities + relations) lets an LLM disambiguate a model even when no prose card exists.
"""
struct ModelRecord
    model::String
    summary::String
    hamiltonian::String
    quantities::Vector{String}
    relations::Vector{Relation}
    references::Vector{String}
end

"""
    describe(model) -> Vector{ModelRecord}

The full record for each model matching `model` (a Type, or a fuzzy Symbol/String): summary +
Hamiltonian where carded, the distinct quantities available, the graph neighborhood, and the
aggregated references — the record an LLM grounds a physics statement on. See [`describe_jsonl`](@ref).
"""
function describe(model)
    out = ModelRecord[]
    for M in _model_types(model)
        hits = search(; model=M).hits
        card = about(M)
        quantities = sort(unique(h.quantity for h in hits))
        refs = String[]
        for h in hits
            append!(refs, h.references)
        end
        card === nothing || append!(refs, card.references)
        push!(
            out,
            ModelRecord(
                _label(M),
                card === nothing ? "" : card.summary,
                card === nothing ? "" : card.hamiltonian,
                quantities,
                relations(M),
                sort(unique(refs)),
            ),
        )
    end
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
        ",\"family\":",
        _jstr(h.family),
        ",\"method\":",
        _jstr(h.method),
        ",\"status\":",
        _jstr(h.status),
        ",\"reliability\":",
        _jstr(h.reliability),
        ",\"thermal\":",
        _jstr(h.thermal),
        ",\"dynamical\":",
        _jstr(h.dynamical),
        ",\"cross_checked\":",
        h.cross_checked,
        ",\"references\":",
        _jarr(h.references),
        ",\"params\":",
        _jarr(h.params),
        ",\"notes\":",
        _jstr(h.notes),
        ",\"valid_domain\":",
        _jstr(h.valid_domain),
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
cross_checked, references, and `params`/`notes`/`valid_domain` — the fetch call signature and
where the value holds). `verified` counts hits that are cross-checked or exact/universal.
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

function _json_record(r::ModelRecord)
    rels = string('[', join((_json_relation(x) for x in r.relations), ','), ']')
    return string(
        "{\"model\":",
        _jstr(r.model),
        ",\"summary\":",
        _jstr(r.summary),
        ",\"hamiltonian\":",
        _jstr(r.hamiltonian),
        ",\"quantities\":",
        _jarr(r.quantities),
        ",\"relations\":",
        rels,
        ",\"references\":",
        _jarr(r.references),
        "}",
    )
end

"""
    describe_jsonl([io=stdout], model) -> nothing

Stream [`describe`](@ref) as JSONL: a `{count, query}` header line then one rich model-record object
per line (summary, hamiltonian, quantities, relations, references).
"""
function describe_jsonl(io::IO, model)
    recs = describe(model)
    print(io, "{\"count\":", length(recs), ",\"query\":", _json_query((; model)), "}\n")
    for r in recs
        print(io, _json_record(r), "\n")
    end
    return nothing
end
describe_jsonl(model) = describe_jsonl(stdout, model)

# ---- realizing: the inverse of a model's :realizes edge — which models realize a class ----

"""
    realizing(class) -> Vector{Relation}

The models that realize a universality `class` (a Symbol/String, case-insensitive exact match) — the
inverse of the model→class `:realizes` edge `relations` returns. Answers "which physical models can I
study the Ising / Heisenberg / KPZ / … class with". See [`realizing_jsonl`](@ref).
"""
function realizing(class)
    want = lowercase(string(class))
    rels = Relation[]
    for r in REALIZES
        lowercase(string(r.class)) == want || continue
        push!(
            rels,
            Relation(
                :realizes, _label(r.model), string(r.class), r.regime, copy(r.references)
            ),
        )
    end
    sort!(rels; by=x -> x.from)
    return rels
end

"""
    realizing_jsonl([io=stdout], class) -> nothing

Stream [`realizing`](@ref) as JSONL: a `{available, query, count}` summary line then one `:realizes`
Relation per line.
"""
function realizing_jsonl(io::IO, class)
    rels = realizing(class)
    print(
        io,
        "{\"available\":",
        !isempty(rels),
        ",\"query\":",
        _json_query((; class)),
        ",\"count\":",
        length(rels),
        "}\n",
    )
    for r in rels
        print(io, _json_relation(r), "\n")
    end
    return nothing
end
realizing_jsonl(class) = realizing_jsonl(stdout, class)

# ---- query_schema: the query convention, self-describing (discover what you can ask) ----

"""
    Facet

One searchable facet of [`search`](@ref): its `name`, `kind`
(`:structural`/`:axis`/`:provenance`/`:derived`), the `match` rule
(`:type_or_fuzzy`/`:enum`/`:fuzzy`/`:bool`), and the enumerated `values` where the
domain is finite (`String[]` when open, e.g. free-text model/quantity names).
"""
struct Facet
    name::Symbol
    kind::Symbol
    match::Symbol
    values::Vector{String}
end

"""
    query_schema() -> Vector{Facet}

The self-describing catalog of [`search`](@ref)'s facets — every facet, its kind,
its match rule, and (where finite) its valid values, drawn live from the registry
(the registered `bc` / `family` / `method` / `reliability` sets) and the fixed
axis / status vocabularies. Lets a consumer (an LLM, a UI) discover *what it can
ask, and with which values*, without reading source — the anti-hallucination
principle applied to the query layer itself. See [`query_schema_jsonl`](@ref).
"""
function query_schema()
    reg_vals(f) = sort(unique(String(f(i)) for i in REGISTRY))
    return Facet[
        Facet(:model, :structural, :type_or_fuzzy, String[]),
        Facet(:quantity, :structural, :type_or_fuzzy, String[]),
        Facet(
            :bc, :structural, :type_or_fuzzy, sort(unique(_label(i.bc) for i in REGISTRY))
        ),
        Facet(:family, :structural, :enum, reg_vals(i -> quantity_family(i.quantity))),
        Facet(:thermal, :axis, :enum, ["zero", "finite", "both", "unknown"]),
        Facet(:dynamical, :axis, :enum, ["static", "transport", "dynamic", "unknown"]),
        Facet(:regime, :derived, :enum, sort(String.(collect(keys(REGIMES))))),
        Facet(:status, :provenance, :enum, sort(String.(collect(STATUS_VALUES)))),
        Facet(:reliability, :provenance, :enum, reg_vals(i -> i.reliability)),
        Facet(:method, :provenance, :fuzzy, reg_vals(i -> i.method)),
        Facet(:reference, :provenance, :fuzzy, String[]),
        Facet(:cross_checked, :provenance, :bool, ["true", "false"]),
    ]
end

function _json_facet(f::Facet)
    return string(
        "{\"name\":",
        _jstr(f.name),
        ",\"kind\":",
        _jstr(f.kind),
        ",\"match\":",
        _jstr(f.match),
        ",\"values\":",
        _jarr(f.values),
        "}",
    )
end

"""
    query_schema_jsonl([io=stdout]) -> nothing

Stream [`query_schema`](@ref) as JSONL: a `{"facets":N}` header then one Facet
object per line (`name`, `kind`, `match`, `values`).
"""
function query_schema_jsonl(io::IO=stdout)
    fs = query_schema()
    print(io, "{\"facets\":", length(fs), "}\n")
    for f in fs
        print(io, _json_facet(f), "\n")
    end
    return nothing
end
