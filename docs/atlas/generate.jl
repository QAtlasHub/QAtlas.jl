# docs/atlas/generate.jl
# docs-as-view (repo-wide). Pure VIEW over {all *_registry.jl claims} +
# {whole test/ INVENTORY}. No src change, no test execution, fetch/register
# framework untouched.
# R1 (ratified, Discussion #379): five-level assurance taxonomy
#   {universality-corroborated, corroborated-at-p, coherent, cited-only,
#    uncorroborated-but-feasible} with an ED-feasible / ED-infeasible
#   denominator split. Only `uncorroborated-but-feasible` is actionable;
#   `cited-only` is the honest frontier, not a penalty.
# R5 (ratified): Model/Quantity/BC@regime is the locked schema. Faceted
#   indices under docs/src/atlas/by/ (model / quantity / bc / level /
#   mechanism / regime); full-text search is Documenter's built-in bar.
# Step 4a/4b kept: one per-hub page (docs/src/atlas/hubs/<slug>.md) with
# the exact reconstructed verify(...) call; index drills down by model.
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "test", "harness", "atlas", "AtlasInventory.jl"))
include(joinpath(ROOT, "test", "harness", "atlas", "AtlasRegistry.jl"))
using .AtlasInventory, .AtlasRegistry

cards = AtlasInventory.scan_dir(joinpath(ROOT, "test"))
isempty(AtlasInventory.PARSE_FAILS) || error(
    "AtlasInventory PARSE_FAILS — refusing to write a partial INVENTORY or ",
    "emit a misleadingly-incomplete view: ",
    AtlasInventory.PARSE_FAILS,
)
AtlasInventory.write_inventory(joinpath(ROOT, "test", "INVENTORY.jsonl"), cards)

regfiles = String[]
for (root, _, fs) in walkdir(joinpath(ROOT, "src"))
    for f in fs
        endswith(f, "_registry.jl") && push!(regfiles, joinpath(root, f))
    end
end
claims = AtlasRegistry.Claim[]
regfail = String[]
for rf in regfiles
    try
        append!(claims, AtlasRegistry.scan_registry(rf))
    catch err
        rel = replace(rf, ROOT * "/" => "")
        @warn "registry parse failed (atlas will abort below)" file = rel exception = err
        push!(regfail, rel)
    end
end
isempty(regfail) || error(
    "Registry parse failures — refusing to emit a misleadingly-empty atlas ",
    "(a malformed *_registry.jl silently drops its whole claim set): ",
    regfail,
)

bycard = Dict{String,Vector{AtlasInventory.Card}}()
for c in cards
    push!(get!(bycard, c.hub, AtlasInventory.Card[]), c)
end
cardsof(h) = get(bycard, h, AtlasInventory.Card[])
mechsof(h) = Set(c.mechanism for c in cardsof(h))

claimed = sort(unique(c.hub for c in claims))
modelof(h) = first(split(h, "/"))
quantof(h) = (p=split(h, "/"); length(p) >= 2 ? p[2] : "?")
bcof(h) = (p=split(h, "/"); length(p) >= 3 ? p[3] : "?")
slugof(h) = replace(h, r"[^A-Za-z0-9]" => "_")
_quant_slugof(q::AbstractString) = replace(q, r"[^A-Za-z0-9]" => "_")
function _md_escape_dollar(s::AbstractString)
    s = replace(s, r"\[([^\]]+)\]\(@ref[^\)]*\)" => s"`\1`")
    return replace(s, r"(?<!\\)\$" => "\\\$")
end

# ── R1 taxonomy (single source of truth lives in AtlasInventory) ─────
# generate.jl is now a thin TYPED consumer of AtlasInventory's
# assurance_level / level_display: a mistyped level is a compile/lookup
# error, not a silently mis-bucketed hub. levelcode is memoised (was
# recomputed 15N+ times via the filters and per-model counts).
ed_infeasible(h) = modelof(h) in AtlasInventory.ED_INFEASIBLE_MODELS
const _LEVEL_CACHE = Dict{String,AtlasInventory.AssuranceLevel}()
function levelcode(h)
    get!(
        () -> AtlasInventory.assurance_level(mechsof(h), ed_infeasible(h)), _LEVEL_CACHE, h
    )
end
levelof(h) = AtlasInventory.level_display(levelcode(h))
levname(h) = levelof(h)[1]
badgeof(h) = levelof(h)[2]

models = sort(unique(modelof(h) for h in claimed))
let stale = setdiff(AtlasInventory.ED_INFEASIBLE_MODELS, Set(models))
    isempty(stale) || @warn string(
        "ED_INFEASIBLE_MODELS has entries absent from claimed models ",
        "(stale / renamed in src?) — risk denominator may be skewed: ",
        sort(collect(stale)),
    )
end
clby = Dict{String,Vector{AtlasRegistry.Claim}}()
for c in claims
    push!(get!(clby, c.hub, AtlasRegistry.Claim[]), c)
end

L_UNIV = filter(h -> levelcode(h) == AtlasInventory.UNIVERSALITY_CORROBORATED, claimed)
L_EDP = filter(h -> levelcode(h) == AtlasInventory.CORROBORATED_AT_P, claimed)
L_COH = filter(h -> levelcode(h) == AtlasInventory.COHERENT, claimed)
L_CITED = filter(h -> levelcode(h) == AtlasInventory.CITED_ONLY, claimed)
L_RISK = filter(h -> levelcode(h) == AtlasInventory.UNCORROBORATED_BUT_FEASIBLE, claimed)

ed_feasible_claimed = filter(h -> !ed_infeasible(h), claimed)
nfeas = length(ed_feasible_claimed)
n_struct = length(L_UNIV) + length(L_EDP)            # external independent
n_inrepo = n_struct + length(L_COH)                  # any executed card
rate_struct = nfeas == 0 ? 0.0 : round(100 * n_struct / nfeas; digits=1)
rate_inrepo = nfeas == 0 ? 0.0 : round(100 * n_inrepo / nfeas; digits=1)

# ── R5 facets ────────────────────────────────────────────────────────
function facet_link(h)
    string(badgeof(h), " [`", h, "`](../hubs/", slugof(h), ".md) — ", levname(h))
end
function group_by(keyfn)
    g = Dict{String,Vector{String}}()
    for h in claimed
        for k in keyfn(h)
            push!(get!(g, string(k), String[]), h)
        end
    end
    return g
end
G_model = group_by(h -> [modelof(h)])
G_quant = group_by(h -> [quantof(h)])
G_bc = group_by(h -> [bcof(h)])
G_level = group_by(h -> [levname(h)])
G_mech = group_by(h -> (M=collect(mechsof(h)); isempty(M) ? ["(no card)"] : M))
G_regime = group_by(
    h -> (r=unique(c.regime for c in cardsof(h)); isempty(r) ? ["(no card)"] : r)
)

const BANNER = string(
    "!!! note \"Provisional v2 view — RES not wired\"\n",
    "    Generated by `docs/atlas/generate.jl` — a pure VIEW over the ",
    "`*_registry.jl` claims + the static `test/INVENTORY.jsonl` AST ",
    "scan. **No test is executed and no `src` is run**; ",
    "`test/INVENTORY.jsonl` is regenerated in-place (idempotently) from ",
    "that static scan; `fetch`/`@register` untouched. Assurance labels ",
    "are PROVISIONAL: residuals / confidence are not shown yet (RES not ",
    "wired). Badges reflect the **committed test AST**, not the latest ",
    "CI run — a hub can read green while its `@test` is red between ",
    "regenerations. `@sweep` = a graceful regime-resolution gap, not ",
    "card omission.",
)

const LEGEND = string(
    "!!! info \"Assurance taxonomy (R1, ratified — Discussion #379)\"\n",
    "    Five honest levels, highest achieved tier wins:\n\n",
    "    - 🟣 **universality-corroborated** — agreement checked against ",
    "the universality class (exponents / central charge). Strongest, ",
    "regime-spanning evidence.\n",
    "    - 🟢 **corroborated-at-p** — an *independent* computation ",
    "(finite-size ED extrapolation or a second closed form) reproduces ",
    "the value at concrete parameter point(s) p.\n",
    "    - 🔵 **coherent** — an independent in-repo card exists and the ",
    "value satisfies an internal invariant (sum rule / limiting case / ",
    "delegation / retype, or an unrecognised/missing route), but no ",
    "external value re-derives it.\n",
    "    - ⚪ **cited-only** — backed only by a literature citation, or ",
    "the model is ED-infeasible so a citation is the ceiling. The ",
    "honest frontier — *neutral, not a penalty*.\n",
    "    - 🟠 **uncorroborated-but-feasible** — `src` claims the hub, ",
    "dense ED *is* feasible, yet no card checks it. **The only ",
    "actionable risk.**\n\n",
    "    Denominator split: the corroboration rate is taken over ",
    "ED-*feasible* claimed hubs only. ED-infeasible models ",
    "(`",
    join(sort(collect(AtlasInventory.ED_INFEASIBLE_MODELS)), "`, `"),
    "`) ",
    "are excluded from the risk denominator — their ceiling is the ",
    "published / DMRG value.",
)

# HOISTED: TIER-1 EXT HELPERS

# Quantity definition extraction: scan src/**/*.jl for docstrings preceding
# `struct <Name>[{params}] <: AbstractQuantity`, build a Dict{base_name => first paragraph}.
const _QUANTITY_DEFS = let
    defs = Dict{String,String}()
    for (root, _, fs) in walkdir(joinpath(ROOT, "src"))
        for f in fs
            endswith(f, ".jl") || continue
            txt = read(joinpath(root, f), String)
            for m in eachmatch(
                r"\"\"\"((?:(?!\"\"\").)+)\"\"\"\s*\n\s*struct\s+([A-Z][A-Za-z0-9_]*)(?:\{[^}]*\})?\s*<:\s*AbstractQuantity"s,
                txt,
            )
                docblock = strip(m.captures[1])
                name = m.captures[2]
                haskey(defs, name) && continue
                lines = split(docblock, "\n")
                body_start = 1
                for (i, ln) in enumerate(lines)
                    if i > 1 && isempty(strip(ln))
                        body_start = i + 1
                        break
                    end
                end
                if body_start > length(lines)
                    defs[name] = ""
                    continue
                end
                para = String[]
                for ln in lines[body_start:end]
                    isempty(strip(ln)) && !isempty(para) && break
                    isempty(strip(ln)) && continue
                    push!(para, strip(ln))
                end
                defs[name] = join(para, " ")
            end
        end
    end
    defs
end

function _quantity_base_name(q::AbstractString)
    return replace(q, r"\{.*\}" => "")
end

# Hub-string normalization for orphan-card matching.
function _normalize_hub(h::AbstractString)
    parts = split(h, "/")
    length(parts) == 3 || return String(h)
    m, q, b = parts[1], parts[2], parts[3]
    m = replace(m, r"^QAtlas\." => "")
    q = replace(q, r"^QAtlas\." => "")
    b = replace(b, r"^QAtlas\." => "")
    q = replace(q, r"\{[^}]*\}" => "")
    return string(m, "/", q, "/", b)
end

# ── TIER-1 EXT HELPERS (universality + calc + refs) ─────────────────────────

# CONVENTION block extraction: parse the standard header comment
#   # CONVENTION
#   #   Hamiltonian: ...
#   #   Observable: ...
#   #   Reference: ...
# from src/models/<class>/<Model>/<Model>.jl.  The block is enforced by
# the CI lint (`test/lint/`) on new model files; older files (e.g.
# IsingSquare) may not have it — return an empty list in that case.
const _STRUCT_PATH_CACHE = Dict{String,String}()
const _STRUCT_CACHE_BUILT = Ref(false)

function _build_struct_path_cache()
    _STRUCT_CACHE_BUILT[] && return nothing
    for (root, _, fs) in walkdir(joinpath(ROOT, "src", "models"))
        for f in fs
            endswith(f, ".jl") || continue
            endswith(f, "_registry.jl") && continue
            p = joinpath(root, f)
            txt = read(p, String)
            for m in eachmatch(r"^struct\s+([A-Z][A-Za-z0-9_]*)(?:\{[^}]*\})?\s*<:"m, txt)
                name = m.captures[1]
                haskey(_STRUCT_PATH_CACHE, name) || (_STRUCT_PATH_CACHE[name] = p)
            end
        end
    end
    _STRUCT_CACHE_BUILT[] = true
    return nothing
end

function _convention_path(model_name::AbstractString)
    for class in ("classical", "quantum")
        p = joinpath(ROOT, "src", "models", class, model_name, model_name * ".jl")
        isfile(p) && return p
    end
    _build_struct_path_cache()
    return get(_STRUCT_PATH_CACHE, model_name, "")
end

function _parse_convention(model_name::AbstractString)
    path = _convention_path(model_name)
    isempty(path) && return Pair{String,String}[]
    lines = readlines(path)
    out = Pair{String,String}[]
    in_block = false
    for ln in lines
        s = rstrip(ln)
        if !in_block
            if startswith(s, "# CONVENTION")
                in_block = true
            end
            continue
        end
        if isempty(s) || !startswith(s, "#")
            break
        end
        m = match(r"^#\s{2,}([^:]+?):\s*(.*)$", s)
        m === nothing && break
        push!(out, strip(m.captures[1]) => strip(m.captures[2]))
    end
    return out
end

# Universality detection: scan the model's CriticalExponents claim notes/refs
# for known QAtlas universality-class tokens (src/universalities/*.jl).  Used
# by ModelList; substrate-derived (zero @register annotation added).
const _UNIV_TOKENS = [
    "MeanField",
    "Ising",
    "Potts",
    "XY",
    "Heisenberg",
    "KPZ",
    "Percolation",
    "E8",
    "CardyEntanglement",
    "MinimalModel",
    "WZW",
    "Poisson",
    "RMT",
]

function _universality_of(model_name::AbstractString)
    for c in claims
        c.model == model_name || continue
        c.quantity in ("CriticalExponents", "CentralCharge") || continue
        haystack = string(c.notes, " | ", c.refs)
        for tok in _UNIV_TOKENS
            occursin(tok, haystack) && return tok
        end
    end
    return ""
end

# Calc cross-reference: substring match of model name (case-insensitive,
# canonicalized: strip trailing dimensionality digits / 'D') in calc filename.
const _CALC_DIR = joinpath(ROOT, "docs", "src", "calc")
const _CALC_FILES =
    isdir(_CALC_DIR) ? sort(filter(f -> endswith(f, ".md"), readdir(_CALC_DIR))) : String[]

function _model_calc_aliases(model_name::AbstractString)
    s = lowercase(model_name)
    aliases = String[s]
    m = match(r"^(.+?)(1d|2d|3d)$", s)
    m === nothing || push!(aliases, m.captures[1])
    m2 = match(r"^s\d+(.+?)(1d|2d|3d)?$", s)
    m2 === nothing || push!(aliases, m2.captures[1])
    return unique(aliases)
end

function _calc_files_for_model(model_name::AbstractString)
    aliases = _model_calc_aliases(model_name)
    out = String[]
    for f in _CALC_FILES
        flo = lowercase(f)
        any(a -> length(a) >= 3 && occursin(a, flo), aliases) || continue
        push!(out, f)
    end
    return out
end

# Hub-level (model + quantity) calc match: stricter than the per-model
# matcher.  Returns only those calc/*.md whose filename mentions BOTH
# the model alias AND the quantity name (lower-cased substring).
function _calc_files_for_hub(model_name::AbstractString, quant_name::AbstractString)
    aliases = _model_calc_aliases(model_name)
    qlo = lowercase(quant_name)
    strict = String[]
    for f in _CALC_FILES
        flo = lowercase(f)
        any(a -> length(a) >= 3 && occursin(a, flo), aliases) || continue
        length(qlo) >= 3 && occursin(qlo, flo) || continue
        push!(strict, f)
    end
    return strict
end

function _split_refs(s::AbstractString)
    isempty(s) && return String[]
    return [strip(r) for r in split(s, "|") if !isempty(strip(r))]
end

# Per the Discussion #379 framework: stateless view-generators over the
# fixed substrate (registry + INVENTORY + assurance taxonomy).  No src
# change, no @register annotation added — every derived column comes from
# what's already in the substrate.
const _ALL_BCS_CANON = ["OBC", "PBC", "Infinite"]

function _bc_order(bcs_present)
    canon = filter(b -> b in bcs_present, _ALL_BCS_CANON)
    extra = sort(filter(b -> !(b in _ALL_BCS_CANON), collect(bcs_present)))
    return vcat(canon, extra)
end

# ── per-model index page (Quantity × BC matrix) ──────────────────────
modelsdir = joinpath(ROOT, "docs", "src", "atlas", "models")
mkpath(modelsdir)

function render_model_index(model_name::AbstractString)
    hs = sort(filter(h -> modelof(h) == model_name, claimed))
    qts = sort(unique(quantof(h) for h in hs))
    bcs = _bc_order(unique(bcof(h) for h in hs))
    cell = Dict{Tuple{String,String},String}()
    for h in hs
        cell[(quantof(h), bcof(h))] = h
    end

    io = IOBuffer()
    P(s...) = println(io, string(s...))
    P("# `", model_name, "` — model index")
    P("")
    P(BANNER)
    P("")
    P(
        "All `(Quantity, BC)` hubs `src` claims for **`",
        model_name,
        "`**.  Cells link to the per-hub card; `—` = not yet implemented ",
        "at that BC.  The shape of the matrix is the *gap visualisation*: ",
        "empty cells are where physics could be added next.",
    )
    P("")
    levcounts = Dict{AtlasInventory.AssuranceLevel,Int}()
    for h in hs
        levcounts[levelcode(h)] = get(levcounts, levelcode(h), 0) + 1
    end
    conv = _parse_convention(model_name)
    P("## Convention")
    P("")
    if isempty(conv)
        P(
            "_No `CONVENTION` header found in `src/models/<class>/",
            model_name,
            "/",
            model_name,
            ".jl` (model file may predate the lint; see ",
            "`docs/src/conventions.md` for the project-wide ",
            "convention policy)._",
        )
    else
        P("| Field | Value |")
        P("|---|---|")
        for (k, v) in conv
            P("| ", k, " | ", _md_escape_dollar(v), " |")
        end
    end
    P("")
    P("## Coverage")
    P("")
    P("| Level | Count |")
    P("|---|---|")
    P(
        "| 🟣 universality-corroborated | ",
        get(levcounts, AtlasInventory.UNIVERSALITY_CORROBORATED, 0),
        " |",
    )
    P("| 🟢 corroborated-at-p | ", get(levcounts, AtlasInventory.CORROBORATED_AT_P, 0), " |")
    P("| 🔵 coherent | ", get(levcounts, AtlasInventory.COHERENT, 0), " |")
    P("| ⚪ cited-only | ", get(levcounts, AtlasInventory.CITED_ONLY, 0), " |")
    P(
        "| 🟠 uncorroborated-but-feasible | ",
        get(levcounts, AtlasInventory.UNCORROBORATED_BUT_FEASIBLE, 0),
        " |",
    )
    P("| **total claimed hubs** | **", length(hs), "** |")
    P("")
    if model_name in AtlasInventory.ED_INFEASIBLE_MODELS
        P("!!! note \"ED-infeasible model\"")
        P(
            "    This model is in `ED_INFEASIBLE_MODELS` (true 2D / frontier).  Its `cited-only` hubs are the published ceiling, **not** an actionable gap.",
        )
        P("")
    end
    methods_seen = sort(
        unique(string(first(clby[h]).method) for h in hs if haskey(clby, h))
    )
    refs_seen = String[]
    for h in hs
        haskey(clby, h) || continue
        r = first(clby[h]).refs
        isempty(r) || push!(refs_seen, r)
    end
    refs_unique = sort(unique(refs_seen))
    P(
        "**Methods** (from `@register`, derived): ",
        isempty(methods_seen) ? "—" : join(("`" * m * "`" for m in methods_seen), ", "),
    )
    P("")
    if !isempty(refs_unique)
        P("**References** (aggregated):")
        for r in refs_unique
            P("- ", _md_escape_dollar(r))
        end
        P("")
    end

    P("## Quantity × BC matrix")
    P("")
    if isempty(qts)
        P("_No hubs registered._")
    else
        P("| Quantity | ", join(("`" * b * "`" for b in bcs), " | "), " |")
        P("|---|", repeat("---|", length(bcs)))
        for q in qts
            cells = String[]
            push!(cells, string("[`", q, "`](../quantities/", _quant_slugof(q), ".md)"))
            for bc in bcs
                h = get(cell, (q, bc), "")
                if isempty(h)
                    push!(cells, "—")
                else
                    push!(cells, string(badgeof(h), " [hub](../hubs/", slugof(h), ".md)"))
                end
            end
            P("| ", join(cells, " | "), " |")
        end
    end
    P("")
    # Derivation notes from docs/src/calc/ (substring match on model name).
    cf = _calc_files_for_model(model_name)
    if !isempty(cf)
        P("## Derivation notes")
        P("")
        P("Matched by filename substring (no annotation; substrate-derived):")
        P("")
        for f in cf
            P("- [`", f, "`](../../calc/", f, ")")
        end
        P("")
    end
    P("[← Atlas index](../index.md) · [Model list →](../ModelList.md)")
    return String(take!(io))
end

for m in models
    write(joinpath(modelsdir, m * ".md"), render_model_index(m))
end

# ── per-quantity index page (Model × BC matrix) ──────────────────────
quantsdir = joinpath(ROOT, "docs", "src", "atlas", "quantities")
mkpath(quantsdir)

function render_quantity_index(quant_name::AbstractString)
    hs = sort(filter(h -> quantof(h) == quant_name, claimed))
    mdls = sort(unique(modelof(h) for h in hs))
    bcs = _bc_order(unique(bcof(h) for h in hs))
    cell = Dict{Tuple{String,String},String}()
    for h in hs
        cell[(modelof(h), bcof(h))] = h
    end

    io = IOBuffer()
    P(s...) = println(io, string(s...))
    P("# `", quant_name, "` — quantity index")
    P("")
    P(BANNER)
    P("")
    P(
        "All `(Model, BC)` hubs `src` claims for the **`",
        quant_name,
        "`** observable.  Empty cells = this model doesn't yet have a `",
        quant_name,
        "` registered at that BC — i.e. where this quantity could be added ",
        "to other models.",
    )
    P("")
    base_name = _quantity_base_name(quant_name)
    qdef = get(_QUANTITY_DEFS, base_name, "")
    if !isempty(qdef)
        P("## Definition")
        P("")
        P(_md_escape_dollar(qdef))
        P("")
        if base_name != quant_name
            P(
                "_(extracted from `src/core/quantities.jl` docstring for the base `",
                base_name,
                "`; this page covers the `",
                quant_name,
                "` variant.)_",
            )
            P("")
        else
            P("_(extracted from `src/core/quantities.jl` docstring.)_")
            P("")
        end
    end
    P("## Coverage")
    P("")
    P("- **Models with this quantity registered**: ", length(mdls))
    P("- **Total hubs (Model, BC pairs)**: ", length(hs))
    q_methods = sort(unique(string(first(clby[h]).method) for h in hs if haskey(clby, h)))
    P(
        "- **Methods** (derived from `@register`): ",
        isempty(q_methods) ? "—" : join(("`" * x * "`" for x in q_methods), ", "),
    )
    q_univs = sort(unique(filter(!isempty, [_universality_of(m) for m in mdls])))
    P(
        "- **Universality classes** (where applicable): ",
        isempty(q_univs) ? "—" : join(("`" * u * "`" for u in q_univs), ", "),
    )
    q_citemap = Dict{String,Int}()
    for h in hs
        haskey(clby, h) || continue
        for r in _split_refs(first(clby[h]).refs)
            q_citemap[r] = get(q_citemap, r, 0) + 1
        end
    end
    if !isempty(q_citemap)
        top_refs = sort(collect(keys(q_citemap)); by=k -> (-q_citemap[k], k))
        top_refs_disp = first(top_refs, min(5, length(top_refs)))
        P("")
        P("**Top references** (by hub count):")
        for r in top_refs_disp
            P(
                "- ",
                _md_escape_dollar(r),
                " — ",
                q_citemap[r],
                " hub",
                (q_citemap[r] == 1 ? "" : "s"),
            )
        end
    end
    P("")
    P("## Model × BC matrix")
    P("")
    if isempty(mdls)
        P("_No hubs registered._")
    else
        P("| Model | ", join(("`" * b * "`" for b in bcs), " | "), " |")
        P("|---|", repeat("---|", length(bcs)))
        for mname in mdls
            cells = String[]
            push!(cells, string("[`", mname, "`](../models/", mname, ".md)"))
            for bc in bcs
                h = get(cell, (mname, bc), "")
                if isempty(h)
                    push!(cells, "—")
                else
                    push!(cells, string(badgeof(h), " [hub](../hubs/", slugof(h), ".md)"))
                end
            end
            P("| ", join(cells, " | "), " |")
        end
    end
    P("")
    P("[← Atlas index](../index.md) · [Model list →](../ModelList.md)")
    return String(take!(io))
end

quants_all = sort(unique(quantof(h) for h in claimed))
for q in quants_all
    write(joinpath(quantsdir, _quant_slugof(q) * ".md"), render_quantity_index(q))
end

# ── ModelList.md — top searchable catalog ────────────────────────────
function render_model_list()
    io = IOBuffer()
    P(s...) = println(io, string(s...))
    P("# Model list — searchable catalog")
    P("")
    P(BANNER)
    P("")
    P(
        "Top-level catalog of all **",
        length(models),
        " models** with claimed hubs.  One row per model; the columns are ",
        "*derived* from the existing substrate (registry + INVENTORY + ",
        "ED-feasibility set + R1 assurance taxonomy) — no extra ",
        "annotation per model.  Use the browser's full-text search ",
        "(Ctrl+F) or Documenter's search bar to filter.  Click a model ",
        "name to drill into its `Quantity × BC` matrix.",
    )
    P("")
    P("| Model | Universality | #K | Methods | 🟣 | 🟢 | 🔵 | ⚪ | 🟠 | ED | Regimes (top 3) |")
    P("|---|---|---|---|---|---|---|---|---|---|---|")
    for m in models
        hs = sort(filter(h -> modelof(h) == m, claimed))
        methods_seen = sort(
            unique(string(first(clby[h]).method) for h in hs if haskey(clby, h))
        )
        regs = Set{String}()
        for h in hs
            for c in cardsof(h)
                push!(regs, c.regime)
            end
        end
        regs_sorted = sort(collect(regs))
        regs_disp = first(regs_sorted, min(3, length(regs_sorted)))
        P(
            "| [`",
            m,
            "`](models/",
            m,
            ".md) | ",
            (u=_universality_of(m); isempty(u) ? "—" : "`" * u * "`"),
            " | ",
            length(hs),
            " | ",
            isempty(methods_seen) ? "—" : join(("`" * x * "`" for x in methods_seen), ", "),
            " | ",
            count(h -> levelcode(h) == AtlasInventory.UNIVERSALITY_CORROBORATED, hs),
            " | ",
            count(h -> levelcode(h) == AtlasInventory.CORROBORATED_AT_P, hs),
            " | ",
            count(h -> levelcode(h) == AtlasInventory.COHERENT, hs),
            " | ",
            count(h -> levelcode(h) == AtlasInventory.CITED_ONLY, hs),
            " | ",
            count(h -> levelcode(h) == AtlasInventory.UNCORROBORATED_BUT_FEASIBLE, hs),
            " | ",
            (m in AtlasInventory.ED_INFEASIBLE_MODELS ? "infeasible" : "feasible"),
            " | ",
            isempty(regs_disp) ? "—" : join(("`" * x * "`" for x in regs_disp), ", "),
            " |",
        )
    end
    P("")
    P("## Quantity index")
    P("")
    P(
        "Each quantity has its own `Model × BC` matrix page (gap visualisation across models):",
    )
    P("")
    for q in quants_all
        nmodels = length(unique(modelof(h) for h in claimed if quantof(h) == q))
        P("- [`", q, "`](quantities/", _quant_slugof(q), ".md) — ", nmodels, " models")
    end
    P("")
    P("[← back to the Atlas index](index.md)")
    return String(take!(io))
end

write(joinpath(ROOT, "docs", "src", "atlas", "ModelList.md"), render_model_list())

println(
    "  + Tier-1 view-generators: ",
    length(models),
    " model pages + ",
    length(quants_all),
    " quantity pages + ModelList.md",
)

# ── per-hub pages ────────────────────────────────────────────────────
hubsdir = joinpath(ROOT, "docs", "src", "atlas", "hubs")
mkpath(hubsdir)
for h in claimed
    cl = first(clby[h])
    cs = cardsof(h)
    lev, bdg, adm = levelof(h)
    hio = IOBuffer()
    HP(s...) = println(hio, string(s...))
    HP("# ", bdg, " `", h, "`")
    HP("")
    HP(BANNER)
    HP("")
    HP("!!! ", adm, " \"Assurance level: ", lev, "\"")
    if lev == "uncorroborated-but-feasible"
        HP(
            "    `src` claims this hub and dense ED is feasible, but no ",
            "corroboration card exists. **Actionable**: add a ",
            "`route = :ed_finite_size` or `:second_closed_form` card.",
        )
    elseif lev == "cited-only"
        HP(
            "    Backed only by a literature citation",
            if ed_infeasible(h)
                " (model is ED-infeasible — this is the ceiling)."
            else
                " — no in-repo independent re-derivation yet."
            end,
        )
    elseif lev == "coherent"
        HP(
            "    An independent card exists and the value satisfies an ",
            "internal invariant; no external value re-derives it yet.",
        )
    else
        HP("    Independently corroborated. See the cards below.")
    end
    HP("")
    HP("## `src` claim")
    HP("")
    HP(
        "- method `",
        cl.method,
        "`, reliability `",
        cl.reliability,
        "`",
        isempty(cl.refs) ? "" : string(", refs: ", _md_escape_dollar(cl.refs)),
    )
    isempty(cl.notes) || HP("- ", _md_escape_dollar(cl.notes))
    HP("")
    HP("## Corroboration")
    HP("")
    if isempty(cs)
        HP(
            "_No corroboration card._ ",
            if ed_infeasible(h)
                "Model is ED-infeasible — frontier (cited-only), not a gap."
            else
                "Flagged by the R1 risk-linter (`src` claims this hub, ED is feasible, no independent card)."
            end,
        )
    else
        HP("| regime | mechanism | independence | refs | file |")
        HP("|---|---|---|---|---|")
        for c in cs
            b = c.independence == "structural" ? "🟢 structural" : "🟡 asserted"
            HP(
                "| `",
                c.regime,
                "` | `",
                c.mechanism,
                "` | ",
                b,
                " | ",
                _md_escape_dollar(c.refs),
                " | `",
                c.file,
                "` |",
            )
        end
    end
    if !isempty(cs)
        HP("")
        HP("## Test calls")
        HP("")
        HP(
            "_The exact `verify(...)` call the harness executed for this hub (reconstructed from the test AST):_",
        )
        HP("")
        for c in cs
            HP("```julia")
            HP(c.srctext)
            HP("```")
            HP("")
        end
    end
    HP("")
    HP("## Assurance (provisional)")
    HP("")
    HP("- level: **", lev, "** ", bdg)
    HP(
        "- cards: ",
        length(cs),
        " · model ED-",
        ed_infeasible(h) ? "infeasible (frontier)" : "feasible",
    )
    HP("- RES not wired — measured residuals / confidence are not shown yet.")
    HP("")
    cf_hub = _calc_files_for_hub(modelof(h), quantof(h))
    if !isempty(cf_hub)
        HP("")
        HP("## Derivation note")
        HP("")
        HP("Matched by filename substring (model + quantity); substrate-derived:")
        HP("")
        for f in cf_hub
            HP("- [`", f, "`](../../calc/", f, ")")
        end
    end
    HP("")
    HP(
        "[← Model: `",
        modelof(h),
        "`](../models/",
        modelof(h),
        ".md) · ",
        "[Quantity: `",
        quantof(h),
        "`](../quantities/",
        _quant_slugof(quantof(h)),
        ".md) · ",
        "[Atlas index](../index.md)",
    )
    write(joinpath(hubsdir, slugof(h) * ".md"), String(take!(hio)))
end

# ── faceted index pages (R5) ─────────────────────────────────────────
bydir = joinpath(ROOT, "docs", "src", "atlas", "by")
mkpath(bydir)
function write_facet(fname, title, groups, blurb)
    fio = IOBuffer()
    FP(s...) = println(fio, string(s...))
    FP("# ", title)
    FP("")
    FP(BANNER)
    FP("")
    FP(blurb)
    FP("")
    for k in sort(collect(keys(groups)))
        hs = sort(groups[k])
        FP("## `", k, "` (", length(hs), ")")
        FP("")
        for h in hs
            FP("- ", facet_link(h))
        end
        FP("")
    end
    FP("[← back to the Atlas index](../index.md)")
    write(joinpath(bydir, fname), String(take!(fio)))
end
write_facet(
    "model.md", "Atlas — by model", G_model, "Every `src`-claimed hub grouped by model."
)
write_facet(
    "quantity.md",
    "Atlas — by quantity",
    G_quant,
    "Grouped by the observable (the `Quantity` axis of the locked Model/Quantity/BC schema).",
)
write_facet(
    "bc.md",
    "Atlas — by boundary condition",
    G_bc,
    "Grouped by boundary condition (`Infinite` / `OBC` / `PBC` …).",
)
write_facet(
    "level.md",
    "Atlas — by assurance level",
    G_level,
    "Grouped by the R1 assurance level. `uncorroborated-but-feasible` is the only actionable bucket.",
)
write_facet(
    "mechanism.md",
    "Atlas — by corroboration mechanism",
    G_mech,
    "Grouped by the `route` the verify card used. A hub appears under each mechanism it has a card for.",
)
write_facet(
    "regime.md",
    "Atlas — by regime",
    G_regime,
    "Grouped by the named physical regime resolved from the test call (`@sweep` = loop-variable, not yet a named point).",
)
byidx = IOBuffer()
BI(s...) = println(byidx, string(s...))
BI("# Atlas — faceted search")
BI("")
BI(BANNER)
BI("")
BI(
    "Full-text search is the bar at the top of every page (Documenter ",
    "built-in — indexes every hub and facet page). Faceted indices over ",
    "the locked **Model / Quantity / BC @ regime** schema:",
)
BI("")
BI("- [By model](model.md) — ", length(G_model), " models")
BI("- [By quantity](quantity.md) — ", length(G_quant), " observables")
BI("- [By boundary condition](bc.md) — ", length(G_bc), " BCs")
BI("- [By assurance level](level.md) — R1 taxonomy")
BI("- [By corroboration mechanism](mechanism.md) — verify `route`")
BI("- [By regime](regime.md) — resolved physical regimes")
BI("")
BI("[← back to the Atlas index](../index.md)")
write(joinpath(bydir, "index.md"), String(take!(byidx)))

# ── atlas index ──────────────────────────────────────────────────────
io = IOBuffer()
P(s...) = println(io, string(s...))
P("# QAtlas — Verified Exact-Solution Atlas")
P("")
P(BANNER)
P("")
P(LEGEND)
P("")
P("## Coverage (all models)")
P("")
P("| | count |")
P("|---|---|")
P("| Hubs `src` claims (registry) | ", length(claimed), " |")
P("| ED-feasible claimed (risk denominator) | ", nfeas, " |")
P("| ED-infeasible claimed (frontier, excluded) | ", length(claimed) - nfeas, " |")
P("| 🟣 universality-corroborated | ", length(L_UNIV), " |")
P("| 🟢 corroborated-at-p | ", length(L_EDP), " |")
P("| 🔵 coherent | ", length(L_COH), " |")
P("| ⚪ cited-only (frontier — neutral) | ", length(L_CITED), " |")
P("| 🟠 uncorroborated-but-feasible (**actionable risk**) | ", length(L_RISK), " |")
P("| Inventory cards scanned (whole test/) | ", length(cards), " |")
P(
    "| Registry files parsed | ",
    length(regfiles) - length(regfail),
    " / ",
    length(regfiles),
    " |",
)
P("| Models | ", length(models), " |")
P("")
P(
    "**Externally-corroborated rate** (🟣+🟢 over ED-feasible claimed): **",
    rate_struct,
    "%** · **in-repo-verified rate** (incl. 🔵 coherent): **",
    rate_inrepo,
    "%**",
)
P("")
P("## Browse by facet")
P("")
P(
    "[**Faceted search →**](by/index.md) · ",
    "[by model](by/model.md) · [by quantity](by/quantity.md) · ",
    "[by BC](by/bc.md) · [by level](by/level.md) · ",
    "[by mechanism](by/mechanism.md) · [by regime](by/regime.md). ",
    "Full-text search is the top bar (Documenter built-in).",
)

P("")
P("")
_audit_counts = let
    no_conv = 0
    no_conv_no_file = 0
    for m in models
        cpath = _convention_path(m)
        if isempty(cpath)
            no_conv_no_file += 1
        elseif isempty(_parse_convention(m))
            no_conv += 1
        end
    end
    no_def = 0
    for q in quants_all
        base = _quantity_base_name(q)
        if !(haskey(_QUANTITY_DEFS, base) && !isempty(_QUANTITY_DEFS[base]))
            no_def += 1
        end
    end
    orphan_calc = 0
    for f in _CALC_FILES
        matched = false
        for m in models
            if f in _calc_files_for_model(m)
                matched = true
                break
            end
        end
        matched || (orphan_calc += 1)
    end
    zero_hubs = count(m -> isempty(filter(h -> modelof(h) == m, claimed)), models)
    norm_claims = Set(_normalize_hub(h) for h in claimed)
    orphan_cards = length(
        unique(filter(h -> !(_normalize_hub(h) in norm_claims), [c.hub for c in cards]))
    )
    (; conv=no_conv + no_conv_no_file, def=no_def, orphan_calc, zero_hubs, orphan_cards)
end

P("## Doc-health audit")
P("")
P("Actionable gap surface — see **[Audit](Audit.md)** for the itemised list.")
P("")
P("| Section | Count |")
P("|---|---|")
P("| 1. Models without CONVENTION header | ", _audit_counts.conv, " |")
P("| 2. Quantities without extracted Definition | ", _audit_counts.def, " |")
P("| 3. Orphan calc notes (matched to no model) | ", _audit_counts.orphan_calc, " |")
P("| 4. Models registered but with 0 hubs | ", _audit_counts.zero_hubs, " |")
P("| 5. INVENTORY card hubs with no `@register` claim | ", _audit_counts.orphan_cards, " |")
P("")
P("## Reference & derivation indices")
P("")
P(
    "Three more substrate-derived indices: ",
    "**[Bibliography](Bibliography.md)** — every citation with hub backlinks; ",
    "**[Methods](Methods.md)** — every `@register` `method=:X` value with hub backlinks; ",
    "**[Derivation-note index](CalcIndex.md)** — each `docs/src/calc/*.md` mapped to its model(s).",
)
P("")
P("## Model & Quantity matrices (Zettelkasten layer)")
P("")
P(
    "Each model has a per-model index showing its hubs as a ",
    "`Quantity × BC` matrix; each quantity has the inverse view ",
    "(`Model × BC`).  Empty cells = gap visualisation (physics not ",
    "yet implemented).  Use the **[Model list](ModelList.md)** for a ",
    "searchable top-catalog.",
)
P("")
P("## 🟠 R1 risk-linter — actionable only")
P("")
P(
    "`src` claims the hub, the model is ED-**feasible**, yet zero ",
    "corroboration cards exist. `cited-only` (frontier) and ED-infeasible ",
    "hubs are **not** listed here — they are the honest ceiling, not a gap.",
)
P("")
if isempty(L_RISK)
    P("!!! tip \"No actionable risk\"")
    P("    Every ED-feasible claimed hub has at least one corroboration card.")
else
    P("!!! warning \"", length(L_RISK), " actionable hub(s)\"")
    for h in sort(L_RISK)
        P("    - [`", h, "`](hubs/", slugof(h), ".md)")
    end
end
P("")
P("## Per-model breakdown")
P("")
P("| model | claimed | 🟣 | 🟢 | 🔵 | ⚪ | 🟠 | ED |")
P("|---|---|---|---|---|---|---|---|")
for m in models
    hs = filter(h -> modelof(h) == m, claimed)
    P(
        "| `",
        m,
        "` | ",
        length(hs),
        " | ",
        count(h -> levelcode(h) == AtlasInventory.UNIVERSALITY_CORROBORATED, hs),
        " | ",
        count(h -> levelcode(h) == AtlasInventory.CORROBORATED_AT_P, hs),
        " | ",
        count(h -> levelcode(h) == AtlasInventory.COHERENT, hs),
        " | ",
        count(h -> levelcode(h) == AtlasInventory.CITED_ONLY, hs),
        " | ",
        count(h -> levelcode(h) == AtlasInventory.UNCORROBORATED_BUT_FEASIBLE, hs),
        " | ",
        (m in AtlasInventory.ED_INFEASIBLE_MODELS ? "infeasible" : "feasible"),
        " |",
    )
end
P("")
P("## Hubs (", length(claimed), ") — select to drill down")
P("")
for m in models
    hs = sort(filter(h -> modelof(h) == m, claimed))
    P("### `", m, "` (", length(hs), ")")
    P("")
    for h in hs
        P("- ", badgeof(h), " [`", h, "`](hubs/", slugof(h), ".md) — ", levname(h))
    end
    P("")
end

out = joinpath(ROOT, "docs", "src", "atlas", "index.md")
mkpath(dirname(out))
write(out, String(take!(io)))
println(
    "wrote ",
    out,
    " + ",
    length(claimed),
    " per-hub pages + 7 facet pages  models=",
    length(models),
    " hubs=",
    length(claimed),
    " cards=",
    length(cards),
    " regfail=",
    length(regfail),
    "  R1[univ=",
    length(L_UNIV),
    " edp=",
    length(L_EDP),
    " coh=",
    length(L_COH),
    " cited=",
    length(L_CITED),
    " risk=",
    length(L_RISK),
    "] feas=",
    nfeas,
    " rate_struct=",
    rate_struct,
    " rate_inrepo=",
    rate_inrepo,
)

# ─── Auto-inject "Verified hubs" section into hand-written model pages ──────
# Closes the model-page <-> atlas-hub coherence gap. For each mapped docs
# page, writes a fresh `## Verified hubs` table between START / END
# markers; if markers are absent, the section is appended at end of file.

const MODEL_DOC_MAP = Dict{String,String}(
    "TFIM" => "docs/src/models/quantum/tfim.md",
    "Heisenberg1D" => "docs/src/models/quantum/heisenberg.md",
    "S1Heisenberg1D" => "docs/src/models/quantum/heisenberg.md",
    "HeisenbergXYZ" => "docs/src/models/quantum/heisenberg.md",
    "DMIHeisenberg1D" => "docs/src/models/quantum/heisenberg.md",
    "J1J2Heisenberg1D" => "docs/src/models/quantum/heisenberg.md",
    "MajumdarGhosh" => "docs/src/models/quantum/majumdar_ghosh.md",
    "Hubbard1D" => "docs/src/models/quantum/hubbard1d.md",
    "ExtendedHubbard1D" => "docs/src/models/quantum/hubbard1d.md",
    "KitaevHoneycomb" => "docs/src/models/quantum/kitaev-honeycomb.md",
    "KitaevHeisenberg" => "docs/src/models/quantum/kitaev-honeycomb.md",
    "Kitaev1D" => "docs/src/models/quantum/kitaev1d.md",
    "ToricCode" => "docs/src/models/quantum/toric-code.md",
    "XXZ1D" => "docs/src/models/quantum/xxz.md",
    "S1XXZ1D" => "docs/src/models/quantum/xxz.md",
    "IsingSquare" => "docs/src/models/classical/ising-square.md",
    "IsingTriangular" => "docs/src/models/classical/ising-triangular.md",
    "SixVertex" => "docs/src/models/classical/six-vertex.md",
)

const _PAGE_TO_MODELS = let
    d = Dict{String,Vector{String}}()
    for (m, page) in MODEL_DOC_MAP
        push!(get!(d, page, String[]), m)
    end
    for (_, ms) in d
        sort!(ms)
    end
    d
end

_posix(p::AbstractString) = replace(p, '\\' => '/')

const _ATLAS_HUBS_BEGIN = "<!-- ATLAS:HUBS:START -- auto-generated by docs/atlas/generate.jl. Do not edit by hand; edits between these markers are overwritten on next regen. -->"
const _ATLAS_HUBS_END = "<!-- ATLAS:HUBS:END -->"

function render_hub_section(page_rel, model_names)
    page_dir = dirname(joinpath(ROOT, page_rel))
    io = IOBuffer()
    P = (xs...) -> println(io, xs...)
    relevant = sort(filter(h -> modelof(h) in model_names, claimed))
    atlas_rel = _posix(relpath(joinpath(ROOT, "docs/src/atlas/index.md"), page_dir))
    P(_ATLAS_HUBS_BEGIN)
    P()
    P("## Verified hubs")
    P()
    n = length(model_names)
    subj = n == 1 ? "this model registers" : "these $(n) models register"
    nh = length(relevant)
    hub_word = nh == 1 ? "hub" : "hubs"
    P(
        "In the [Verified Atlas](",
        atlas_rel,
        "), ",
        subj,
        " ",
        nh,
        " ",
        hub_word,
        " (quantity / BC pair). The badge column shows the R1 assurance level; click a hub link to see the exact `verify(...)` calls, references, and corroboration mechanism.",
    )
    P()
    if nh == 0
        P("_No hubs registered yet._")
        P()
    else
        if n == 1
            P("| Quantity | BC | Assurance | Cards |")
            P("|---|---|---|---|")
            for h in relevant
                hub_page = joinpath(ROOT, "docs/src/atlas/hubs/$(slugof(h)).md")
                rp = _posix(relpath(hub_page, page_dir))
                P(
                    "| [`",
                    quantof(h),
                    "`](",
                    rp,
                    ") | `",
                    bcof(h),
                    "` | ",
                    badgeof(h),
                    " ",
                    levname(h),
                    " | ",
                    length(cardsof(h)),
                    " |",
                )
            end
        else
            P("| Model | Quantity | BC | Assurance | Cards |")
            P("|---|---|---|---|---|")
            for h in relevant
                hub_page = joinpath(ROOT, "docs/src/atlas/hubs/$(slugof(h)).md")
                rp = _posix(relpath(hub_page, page_dir))
                P(
                    "| `",
                    modelof(h),
                    "` | [`",
                    quantof(h),
                    "`](",
                    rp,
                    ") | `",
                    bcof(h),
                    "` | ",
                    badgeof(h),
                    " ",
                    levname(h),
                    " | ",
                    length(cardsof(h)),
                    " |",
                )
            end
        end
        P()
    end
    P(_ATLAS_HUBS_END)
    return String(take!(io))
end

function inject_hub_section!(page_rel, model_names)
    p = joinpath(ROOT, page_rel)
    if !isfile(p)
        @warn "ATLAS hub-inject skip: page not found" page = page_rel
        return nothing
    end
    content = read(p, String)
    section = render_hub_section(page_rel, model_names)
    b = findfirst(_ATLAS_HUBS_BEGIN, content)
    e = findfirst(_ATLAS_HUBS_END, content)
    # Explicit marker-pair validation: refuse to write when only one of
    # START/END is present or when they are misordered. Without this an
    # orphaned START at the top of the page would be paired with the
    # appended-section END at the bottom on the NEXT regen, silently
    # nuking everything between them.
    if (b === nothing) != (e === nothing)
        @warn "ATLAS hub-inject skip: page has only one of START/END marker; refusing to mutate to avoid data loss" page =
            page_rel
        return nothing
    end
    if b !== nothing && e !== nothing && first(e) <= last(b)
        @warn "ATLAS hub-inject skip: START/END markers in wrong order; refusing to mutate" page =
            page_rel
        return nothing
    end
    if b !== nothing && e !== nothing
        new = content[1:(first(b) - 1)] * section * content[(last(e) + 1):end]
    else
        new = rstrip(content, '\n') * "\n\n---\n\n" * section * "\n"
    end
    # Normalize trailing newlines so re-runs are byte-stable (without this
    # the append path leaves a trailing "\n\n" that grows by 1 newline on
    # each subsequent replacement-path regen).
    new = rstrip(new, '\n') * "\n"
    if new != content
        write(p, new)
        println("  injected hubs into ", page_rel)
    end
end

for page in sort(collect(keys(_PAGE_TO_MODELS)))
    inject_hub_section!(page, _PAGE_TO_MODELS[page])
end

# ── TIER-1 VIEW GENERATORS (ModelList + per-model + per-quantity) ──────────

# ── TIER-1 EXT EMITTERS (Bibliography + CalcIndex) ──────────────────────────
function render_bibliography()
    citemap = Dict{String,Vector{String}}()
    for c in claims
        for r in _split_refs(c.refs)
            push!(get!(citemap, r, String[]), c.hub)
        end
    end

    io = IOBuffer()
    P(s...) = println(io, string(s...))
    P("# Bibliography — citations across the atlas")
    P("")
    P(BANNER)
    P("")
    P(
        "Every distinct citation string appearing in any `@register(..., ",
        "references=[...], ...)` across the atlas (",
        length(citemap),
        " unique citations).  For each one we list which hubs cite it — a ",
        "*load-bearing-ness* view: high-count entries are central, ",
        "low-count entries are local.",
    )
    P("")
    P(
        "Citations are kept as-is from `@register` (no normalization, no ",
        "DOI lookup).  Free-form strings are an honest substrate; a ",
        "structured citation database is deferred to a follow-up.",
    )
    P("")
    keys_sorted = sort(collect(keys(citemap)); by=k -> (-length(citemap[k]), k))
    for k in keys_sorted
        hs = sort(unique(citemap[k]))
        P(
            "## ",
            _md_escape_dollar(k),
            " — ",
            length(hs),
            " hub",
            (length(hs) == 1 ? "" : "s"),
        )
        P("")
        for h in hs
            P("- ", badgeof(h), " [`", h, "`](hubs/", slugof(h), ".md)")
        end
        P("")
    end
    P("[← back to the Atlas index](index.md)")
    return String(take!(io))
end

write(joinpath(ROOT, "docs", "src", "atlas", "Bibliography.md"), render_bibliography())

function render_calc_index()
    io = IOBuffer()
    P(s...) = println(io, string(s...))
    P("# Derivation-note index (`docs/src/calc/`)")
    P("")
    P(BANNER)
    P("")
    P(
        "All ",
        length(_CALC_FILES),
        " step-by-step derivation notes under `docs/src/calc/`, with the ",
        "model(s) each matches by filename substring.  Substrate-derived ",
        "(no annotation in calc files or `@register`).",
    )
    P("")
    P("| Derivation note | Models matched |")
    P("|---|---|")
    for f in _CALC_FILES
        matched = String[]
        for m in models
            cf = _calc_files_for_model(m)
            f in cf || continue
            push!(matched, m)
        end
        P(
            "| [`",
            f,
            "`](../calc/",
            f,
            ") | ",
            if isempty(matched)
                "—"
            else
                join(("[`" * m * "`](models/" * m * ".md)" for m in matched), ", ")
            end,
            " |",
        )
    end
    P("")
    P("[← back to the Atlas index](index.md)")
    return String(take!(io))
end

write(joinpath(ROOT, "docs", "src", "atlas", "CalcIndex.md"), render_calc_index())

println(
    "  + Tier-1 ext: Bibliography.md + CalcIndex.md (", length(_CALC_FILES), " calc notes)"
)

# ── TIER-1 EXT AUDIT (Audit.md) ─────────────────────────────────────────────
function render_audit()
    io = IOBuffer()
    P(s...) = println(io, string(s...))
    P("# Audit — doc-health gap surface")
    P("")
    P(BANNER)
    P("")
    P(
        "Substrate-derived audit of actionable gaps.  Each section is a ",
        "concrete to-do list: an item here is either a hand-fixable ",
        "doc/code issue or a tracked physics task.  This page is the one ",
        "place to look for \"what's not great yet\".",
    )
    P("")

    P("## 1. Models without `CONVENTION` header")
    P("")
    P(
        "The CI lint enforces `# CONVENTION` headers on new model files, but ",
        "older files predate the lint.  These models show an absence note on ",
        "their per-model page; backfilling adds a one-block comment.",
    )
    P("")
    no_conv = String[]
    for m in models
        isempty(_convention_path(m)) && continue
        isempty(_parse_convention(m)) || continue
        push!(no_conv, m)
    end
    no_conv_no_file = String[m for m in models if isempty(_convention_path(m))]
    if isempty(no_conv) && isempty(no_conv_no_file)
        P("!!! tip \"All models have a CONVENTION header.\"")
    else
        if !isempty(no_conv)
            P(
                "**Has source file but missing/unparseable `CONVENTION` block** (",
                length(no_conv),
                "):",
            )
            P("")
            for m in no_conv
                P("- [`", m, "`](models/", m, ".md)")
            end
            P("")
        end
        if !isempty(no_conv_no_file)
            P(
                "**Source file not found at `src/models/<class>/<Model>/<Model>.jl`** (",
                length(no_conv_no_file),
                ") — model may live elsewhere or be defined inline:",
            )
            P("")
            for m in no_conv_no_file
                P("- [`", m, "`](models/", m, ".md)")
            end
            P("")
        end
    end

    P("## 2. Quantities without auto-extracted `Definition`")
    P("")
    P(
        "Quantities whose `struct X[{params}] <: AbstractQuantity` docstring ",
        "wasn't matched by the regex extractor (likely defined as bare ",
        "`struct X end` without `<: AbstractQuantity`, or with alternate ",
        "formatting).  Adding the supertype + docstring makes them appear ",
        "on the per-quantity page automatically.",
    )
    P("")
    no_def = String[]
    for q in quants_all
        base = _quantity_base_name(q)
        haskey(_QUANTITY_DEFS, base) && !isempty(_QUANTITY_DEFS[base]) && continue
        push!(no_def, q)
    end
    if isempty(no_def)
        P("!!! tip \"All quantities have an extracted Definition.\"")
    else
        P("**", length(no_def), " quantities**:")
        P("")
        for q in no_def
            P("- [`", q, "`](quantities/", _quant_slugof(q), ".md)")
        end
        P("")
    end

    P("## 3. Orphan calc notes (matched to no model)")
    P("")
    P(
        "`docs/src/calc/*.md` whose filename doesn't substring-match any ",
        "registered model.  Likely true derivation notes that describe a ",
        "method (e.g. `calabrese-cardy-obc-vs-pbc.md`, ",
        "`ad-thermodynamics-from-z.md`) rather than a model, but worth ",
        "scanning to confirm.",
    )
    P("")
    orphan_calc = String[]
    for f in _CALC_FILES
        matched = false
        for m in models
            if f in _calc_files_for_model(m)
                matched = true
                break
            end
        end
        matched || push!(orphan_calc, f)
    end
    if isempty(orphan_calc)
        P("!!! tip \"All calc notes match at least one model.\"")
    else
        P("**", length(orphan_calc), " orphan calc note(s)**:")
        P("")
        for f in orphan_calc
            P("- [`", f, "`](../calc/", f, ")")
        end
        P("")
    end

    P("## 4. Models registered but with 0 hubs")
    P("")
    zero_hubs = String[m for m in models if isempty(filter(h -> modelof(h) == m, claimed))]
    if isempty(zero_hubs)
        P("!!! tip \"Every registered model has at least one hub.\"")
    else
        P("**", length(zero_hubs), " models**:")
        for m in zero_hubs
            P("- `", m, "`")
        end
        P("")
    end

    P("## 5. INVENTORY cards with no matching registry claim")
    P("")
    P(
        "Verify cards exist for `(M, Q, BC)` triples that no `@register` ",
        "claims.  These hubs are visible only through the cards, not the ",
        "registry view; the absence of an `@register` claim means the ",
        "primary `src` claim isn't there.",
    )
    P("")
    norm_claims = Set(_normalize_hub(h) for h in claimed)
    orphan_cards = sort(
        unique(filter(h -> !(_normalize_hub(h) in norm_claims), [c.hub for c in cards]))
    )
    if isempty(orphan_cards)
        P("!!! tip \"All INVENTORY card hubs have a matching `@register` claim.\"")
    else
        P("**", length(orphan_cards), " orphan card hub(s)**:")
        P("")
        for h in orphan_cards
            P("- `", h, "`")
        end
        P("")
    end

    P("[← back to the Atlas index](index.md)")
    return String(take!(io))
end

write(joinpath(ROOT, "docs", "src", "atlas", "Audit.md"), render_audit())
println("  + Tier-1 ext: Audit.md")

# ── TIER-1 EXT METHODS (Methods.md) ──────────────────────────────────────
function render_methods()
    methodmap = Dict{String,Vector{String}}()
    for c in claims
        m = string(c.method)
        isempty(m) && continue
        push!(get!(methodmap, m, String[]), c.hub)
    end

    io = IOBuffer()
    P(s...) = println(io, string(s...))
    P("# Methods — solution mechanisms across the atlas")
    P("")
    P(BANNER)
    P("")
    P(
        "Every distinct `method=:X` value used in any `@register(...)` ",
        "across the atlas (",
        length(methodmap),
        " unique methods).  For each one we list which hubs use it - a ",
        "*solution-technique* view: high-count entries are central ",
        "machinery (analytic, dense_ed, ...), low-count entries are ",
        "model-specific (e.g. `bethe_ansatz` for the integrable models).",
    )
    P("")
    P(
        "Symmetric with **[Bibliography](Bibliography.md)** ",
        "(which groups by citation): both decompose hubs along orthogonal ",
        "axes of the `@register` metadata.",
    )
    P("")
    keys_sorted = sort(collect(keys(methodmap)); by=k -> (-length(methodmap[k]), k))
    for k in keys_sorted
        hs = sort(unique(methodmap[k]))
        P("## `", k, "` - ", length(hs), " hub", (length(hs) == 1 ? "" : "s"))
        P("")
        for h in hs
            P("- ", badgeof(h), " [`", h, "`](hubs/", slugof(h), ".md)")
        end
        P("")
    end
    P("[← back to the Atlas index](index.md)")
    return String(take!(io))
end

write(joinpath(ROOT, "docs", "src", "atlas", "Methods.md"), render_methods())
println("  + Tier-1 ext: Methods.md")
