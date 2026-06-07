# docs/atlas/graph_export.jl — render the QAtlas knowledge graph as an
# Obsidian-style force-directed graph view inside the Documenter site, using the
# Quartz graph engine (d3-force + PixiJS) ported to a standalone, data-driven
# bundle: docs/src/assets/qatlas-graph.js (built from
# obsidian-remote-ssh/docs-site/qatlas-graph-entry.ts via esbuild).
#
# Encoding (deliberately minimal — two essential axes):
#   * node color  = node type (model / universality class / bound domain /
#                    quantity)
#   * edge color  = verified status of the result a node provides: exact &
#                    universal share one color, bound, approx; `realizes`
#                    (model → class) edges use a neutral color
#   * edge style  = solid (a test validates it) / dashed (no dedicated test)
#
# Two relations only: a model *realizes* a class (model ↔ universality), and a
# namespace *provides* a quantity (quantity ↔ verified status).  Coherence gaps
# are NOT drawn — an isolated node (an undeveloped class, a model realizing no
# class) *is* the gap, visible by its lack of edges.
#
# Run: julia --project=. docs/atlas/graph_export.jl  →  docs/src/atlas/graph.md

using QAtlas

const REG = QAtlas.REGISTRY
const REAL = QAtlas.REALIZES

_short(T) = replace(string(T), "QAtlas." => "")
_js(s) = replace(string(s), "\\" => "\\\\", "\"" => "\\\"")
_clean(s) = occursin(r"^[A-Za-z0-9_]+$", s)

# ── collect nodes + edges (only nodes that appear in an edge) ────────────────
nodes = Dict{String,NamedTuple}()      # id => (text, group, url)
edges = Vector{NamedTuple}()           # (from, to, kind, status, verified)

addnode!(id, text, group, url) = get!(nodes, id, (text=text, group=group, url=url))

_is_uni(T) = T <: QAtlas.Universality
_is_bnd(T) = T <: QAtlas.Bound
_class(T) = T.parameters[1]
_concrete(T) = !_is_uni(T) && !_is_bnd(T)

_model_url(name) = _clean(name) ? "../models/$(name)/" : ""
_quantity_url(name) = _clean(name) ? "../quantities/$(name)/" : ""

# realizes : model → class  (membership is asserted, so always solid)
for r in REAL
    m = _short(r.model)
    addnode!("M:" * m, m, "model", _model_url(m))
    addnode!("C:" * string(r.class), string(r.class), "class", "")
    push!(
        edges,
        (
            from="M:" * m,
            to="C:" * string(r.class),
            kind="realizes",
            status="exact",
            verified=true,
        ),
    )
end

# provides : namespace → quantity, one per registry row.  The namespace is the
# universality class (`:universal`), the bound domain (`:bound`), or the
# concrete model (`:exact`/`:bound`/`:approx`); verified iff a test backs it.
for e in REG
    q = _short(e.quantity)
    addnode!("Q:" * q, q, "quantity", _quantity_url(q))
    ver = e.tested_in !== nothing
    if _is_uni(e.model)
        c = _class(e.model)
        addnode!("C:" * string(c), string(c), "class", "")
        push!(
            edges,
            (
                from="C:" * string(c),
                to="Q:" * q,
                kind="provides",
                status="universal",
                verified=ver,
            ),
        )
    elseif _is_bnd(e.model)
        d = _class(e.model)
        addnode!("B:" * string(d), string(d), "bound", "")
        push!(
            edges,
            (
                from="B:" * string(d),
                to="Q:" * q,
                kind="provides",
                status=string(e.status),
                verified=ver,
            ),
        )
    else
        m = _short(e.model)
        addnode!("M:" * m, m, "model", _model_url(m))
        push!(
            edges,
            (
                from="M:" * m,
                to="Q:" * q,
                kind="provides",
                status=string(e.status),
                verified=ver,
            ),
        )
    end
end

# merge edges by (from, to, kind): a (namespace, quantity) appears at several
# bc/scheme rows but is one edge — keep the strongest claim's color and mark it
# solid if *any* of its definitions is tested.
const _STATUS_RANK = Dict("exact" => 4, "universal" => 4, "bound" => 3, "approx" => 2)
_rank(s) = get(_STATUS_RANK, s, 1)

merged = Dict{Tuple{String,String,String},NamedTuple}()
for e in edges
    k = (e.from, e.to, e.kind)
    if haskey(merged, k)
        p = merged[k]
        best = _rank(e.status) > _rank(p.status) ? e.status : p.status
        merged[k] = (
            from=e.from,
            to=e.to,
            kind=e.kind,
            status=best,
            verified=p.verified || e.verified,
        )
    else
        merged[k] = e
    end
end
edges = sort!(collect(values(merged)); by=e -> (e.from, e.to, e.kind))

# ── emit QAtlasGraph JSON ─────────────────────────────────────────────────────
function node_json(id, n)
    "{id:\"$(_js(id))\",text:\"$(_js(n.text))\",group:\"$(n.group)\",url:\"$(_js(n.url))\"}"
end
nodes_js = join(
    [node_json(id, nodes[id]) for id in sort!(collect(keys(nodes)))], ",\n      "
)
function edge_json(e)
    "{source:\"$(_js(e.from))\",target:\"$(_js(e.to))\"," *
    "kind:\"$(e.kind)\",status:\"$(e.status)\",verified:$(e.verified)}"
end
edges_js = join([edge_json(e) for e in edges], ",\n      ")

n_models = count(id -> startswith(id, "M:"), keys(nodes))
n_classes = count(id -> startswith(id, "C:"), keys(nodes))
n_bounds = count(id -> startswith(id, "B:"), keys(nodes))
n_quant = count(id -> startswith(id, "Q:"), keys(nodes))

page = """
# Knowledge graph

The QAtlas vault as a force-directed network of two relations: a **model**
*belongs to* a **universality class** (model ↔ universality), and a namespace —
model, class, or **bound** domain — *provides* a **quantity** (quantity ↔
verified status).

Node color marks the node type.  Edge color marks the kind of claim the
provided result makes — **exact / universal** (one color), **bound**, or
**approx** — and edge style marks verification: a **solid** edge has a
dedicated test, a **dashed** edge does not.  Coherence *gaps* are not drawn:
an isolated node (an undeveloped class, or a model belonging to no class) *is*
the gap, visible by its lack of edges.

Type in the search box to highlight nodes by name.  Drag nodes, scroll to zoom,
hover to highlight neighbours and reveal labels, click a model or quantity to
open its page.

*Rendered by the Quartz graph engine (d3-force + PixiJS), ported to a standalone
bundle (`assets/qatlas-graph.js`).  Generated by `docs/atlas/graph_export.jl`
from the live registry — $(length(nodes)) nodes ($(n_models) models,
$(n_classes) classes, $(n_bounds) bound domains, $(n_quant) quantities),
$(length(edges)) edges.*

```@raw html
<div id="qatlas-graph" style="width:100%;height:720px;border:1px solid var(--light-color,#ddd);border-radius:6px;"></div>
<script>
(function(){
  var DATA = {
    nodes: [
      $(nodes_js)
    ],
    links: [
      $(edges_js)
    ]
  };
  var CFG = {
    drag: true, zoom: true, repelForce: 0.85, centerForce: 0.2,
    linkDistance: 105, linkStrength: 0.07, fontSize: 0.62, opacityScale: 1.1,
    focusOnHover: true, realizesLabel: "belongs to",
    labelColor: "#eaeaea", labelStroke: "#15171a"
  };
  function go(){
    var el = document.getElementById("qatlas-graph");
    if (!el) return;
    if (window.QAtlasGraph) { window.QAtlasGraph.render(el, DATA, CFG); return; }
    var root = location.pathname.split("atlas/graph")[0];
    var s = document.createElement("script");
    s.src = root + "assets/qatlas-graph.js";
    s.onload = function(){ if (window.QAtlasGraph) window.QAtlasGraph.render(el, DATA, CFG); };
    s.onerror = function(){ el.innerHTML = "<p style='padding:1em'>graph bundle failed to load from " + s.src + "</p>"; };
    document.head.appendChild(s);
  }
  if (document.readyState === "complete") { go(); }
  else { window.addEventListener("load", go); }
})();
</script>
```
"""

open(joinpath(@__DIR__, "..", "src", "atlas", "graph.md"), "w") do io
    write(io, page)
end
println("graph.md written: ", length(nodes), " nodes, ", length(edges), " edges")
