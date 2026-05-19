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
    catch
        push!(regfail, replace(rf, ROOT * "/" => ""))
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

# ── R1 taxonomy ──────────────────────────────────────────────────────
# Models where dense exact diagonalisation at a physically meaningful
# size is computationally infeasible: the published / DMRG literature
# value is the ceiling, so the absence of an in-repo ED card is the
# honest frontier (cited-only), NOT an actionable gap.
const ED_INFEASIBLE_MODELS = Set([
    "KagomeHeisenbergAFM",
    "ToricCode",
    "XCube",
    "SYK",
    "ChernSimons3D",
    "FibonacciAnyons",
    "PpIp2DSC",
    "AKLT2D",
    "KitaevHoneycomb",
])
ed_infeasible(h) = modelof(h) in ED_INFEASIBLE_MODELS

const MECH_UNIV = Set(["universality_consistency"])
const MECH_EDP = Set(["ed_finite_size", "second_closed_form"])
const MECH_COH = Set([
    "delegation_invariant", "limiting_case", "sum_rule", "retype_formula", "unknown"
])
const MECH_CITED = Set(["literature_value"])

# Highest achieved tier wins. Returns (level, badge, admonition).
function levelof(h)
    M = mechsof(h)
    if !isempty(intersect(M, MECH_UNIV))
        return ("universality-corroborated", "🟣", "tip")
    elseif !isempty(intersect(M, MECH_EDP))
        return ("corroborated-at-p", "🟢", "tip")
    elseif !isempty(intersect(M, MECH_COH))
        return ("coherent", "🔵", "note")
    elseif !isempty(intersect(M, MECH_CITED))
        return ("cited-only", "⚪", "note")
    elseif ed_infeasible(h)
        return ("cited-only", "⚪", "note")
    else
        return ("uncorroborated-but-feasible", "🟠", "warning")
    end
end
const _LEVEL_CACHE = Dict{String,Tuple{String,String,String}}()
levelof_cached(h) = get!(() -> levelof(h), _LEVEL_CACHE, h)
levname(h) = levelof_cached(h)[1]
badgeof(h) = levelof_cached(h)[2]

models = sort(unique(modelof(h) for h in claimed))
let stale = setdiff(ED_INFEASIBLE_MODELS, Set(models))
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

L_UNIV = filter(h -> levname(h) == "universality-corroborated", claimed)
L_EDP = filter(h -> levname(h) == "corroborated-at-p", claimed)
L_COH = filter(h -> levname(h) == "coherent", claimed)
L_CITED = filter(h -> levname(h) == "cited-only", claimed)
L_RISK = filter(h -> levname(h) == "uncorroborated-but-feasible", claimed)

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
    "registry claims + static `test/INVENTORY.jsonl`. No source/test ",
    "executed; `fetch`/`@register` untouched. Assurance labels are ",
    "PROVISIONAL: residuals / confidence are not shown yet (RES not ",
    "wired); `@sweep` = a graceful regime-resolution gap, not card ",
    "omission.",
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
    "delegation), but no external value re-derives it.\n",
    "    - ⚪ **cited-only** — backed only by a literature citation, or ",
    "the model is ED-infeasible so a citation is the ceiling. The ",
    "honest frontier — *neutral, not a penalty*.\n",
    "    - 🟠 **uncorroborated-but-feasible** — `src` claims the hub, ",
    "dense ED *is* feasible, yet no card checks it. **The only ",
    "actionable risk.**\n\n",
    "    Denominator split: the corroboration rate is taken over ",
    "ED-*feasible* claimed hubs only. ED-infeasible models ",
    "(`",
    join(sort(collect(ED_INFEASIBLE_MODELS)), "`, `"),
    "`) ",
    "are excluded from the risk denominator — their ceiling is the ",
    "published / DMRG value.",
)

# ── per-hub pages ────────────────────────────────────────────────────
hubsdir = joinpath(ROOT, "docs", "src", "atlas", "hubs")
mkpath(hubsdir)
for h in claimed
    cl = first(clby[h])
    cs = cardsof(h)
    lev, bdg, adm = levelof_cached(h)
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
        isempty(cl.refs) ? "" : string(", refs: ", cl.refs),
    )
    isempty(cl.notes) || HP("- ", cl.notes)
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
                c.refs,
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
    HP("[← back to the Atlas index](../index.md)")
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
        count(h -> levname(h) == "universality-corroborated", hs),
        " | ",
        count(h -> levname(h) == "corroborated-at-p", hs),
        " | ",
        count(h -> levname(h) == "coherent", hs),
        " | ",
        count(h -> levname(h) == "cited-only", hs),
        " | ",
        count(h -> levname(h) == "uncorroborated-but-feasible", hs),
        " | ",
        (m in ED_INFEASIBLE_MODELS ? "infeasible" : "feasible"),
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
