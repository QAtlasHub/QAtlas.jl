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

One available (model, quantity, bc) hub, with its provenance. `cross_checked` is `true` when the
model participates in an executable cross-check (a duality, limit, or LSM-symmetry edge) — the
QAtlas-specific trust signal (model-level in this version).
"""
struct Hit
    model::String
    quantity::String
    bc::String
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
    search(; model=nothing, quantity=nothing, bc=nothing, regime=nothing) -> QueryResult

Does the atlas have data matching the given facets? Any facet left `nothing` is unconstrained.
`model`/`quantity`/`bc` accept a Type (exact / family) or a Symbol/String (fuzzy substring on the
name); `regime` is one of `keys(REGIMES)` (`:ground_state`, `:finite_temperature`, `:dynamics`,
`:universality`). Returns a [`QueryResult`](@ref) — see [`search_jsonl`](@ref) for JSONL output and
[`available`](@ref) for just the boolean.
"""
function search(; model=nothing, quantity=nothing, bc=nothing, regime=nothing)
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
        push!(
            hits,
            Hit(
                _label(impl.model),
                _label(impl.quantity),
                _label(impl.bc),
                impl.status,
                impl.reliability,
                impl.model in xchecked,
                copy(impl.references),
            ),
        )
    end
    sort!(hits; by=h -> (h.model, h.quantity, h.bc))
    return QueryResult((; model, quantity, bc, regime), !isempty(hits), hits)
end

"""
    available(; model=nothing, quantity=nothing, bc=nothing, regime=nothing) -> Bool

The bare yes/no: does the atlas have anything matching the facets? Convenience over [`search`](@ref).
"""
available(; kwargs...) = search(; kwargs...).available

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
        push!(parts, string(_jstr(k), ":", _jstr(v isa Type ? _label(v) : v)))
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
