# core/links.jl — the knowledge-graph query layer over REGISTRY + REALIZES.
#
# QAtlas as a vault: models, quantities, universality classes, bounds and
# references are nodes; `implements` / `predicts` / `bounds` / `realizes` /
# `delegates` / `cites` are edges.  Every relation is queryable in BOTH
# directions (backlinks), so any node can pull the rest of the network — the
# substrate for graph-derived coherence checks (`coherence.jl`) and the atlas
# graph view.
#
# Pure queries over the existing edge stores (`REGISTRY`, `REALIZES`); no new
# state.  `_as_type` is shared with registry.jl / realizes.jl.

# ── node-type helpers ────────────────────────────────────────────────────────
# `_is_universality` lives in core/universality.jl (register! needs it earlier).
_is_bound(::Type{T}) where {T} = T <: Bound
_class_of(::Type{T}) where {T} = T <: Universality ? T.parameters[1]::Symbol : nothing
_domain_of(::Type{T}) where {T} = T <: Bound ? T.parameters[1]::Symbol : nothing
_kgshort(T) = replace(string(T), "QAtlas." => "")

# ── predicts : universality class → quantities  (and inverse) ─────────────────
"""
    predicts(class::Symbol) -> Vector{NamedTuple}

The quantities universality `class` predicts — its `status=:universal` rows:
`(quantity, bc, scheme)`.  Inverse of [`predicted_by`](@ref).
"""
function predicts(class::Symbol)
    return [
        (quantity=e.quantity, bc=e.bc, scheme=e.scheme) for
        e in REGISTRY if e.model === Universality{class} && e.status === :universal
    ]
end

"""
    predicted_by(quantity) -> Vector{Symbol}

The universality classes whose `:universal` rows predict `quantity`.
"""
function predicted_by(quantity)
    q_T = _as_type(quantity)
    cls = Symbol[]
    for e in REGISTRY
        if e.quantity === q_T && _is_universality(e.model)
            c = _class_of(e.model)
            c === nothing || push!(cls, c)
        end
    end
    return unique(cls)
end

# ── bounds : quantity ↔ bound ────────────────────────────────────────────────
"""
    bounds_on(quantity) -> Vector{NamedTuple}

Every registered bound constraining `quantity`: `(domain, direction, scheme,
references)` — the bound nodes reachable from a quantity node.
"""
function bounds_on(quantity)
    q_T = _as_type(quantity)
    return [
        (
            domain=_domain_of(e.model),
            direction=e.direction,
            scheme=e.scheme,
            references=e.references,
        ) for e in REGISTRY if e.quantity === q_T && e.status === :bound
    ]
end

# ── cites : implementation ↔ reference ───────────────────────────────────────
"""
    cited_by(ref::AbstractString) -> Vector{NamedTuple}

Every implementation row citing `ref`: `(model, quantity, bc, scheme, status)`.
The backlink inverse of the per-row `references` field.
"""
function cited_by(ref::AbstractString)
    return [
        (model=e.model, quantity=e.quantity, bc=e.bc, scheme=e.scheme, status=e.status) for
        e in REGISTRY if ref in e.references
    ]
end

# ── delegates : (model, quantity) → class / model ─────────────────────────────
"""
    delegations(model) -> Vector{NamedTuple}

The quantities `model` computes by delegating (`method=:delegation` rows),
paired with the declared targets of the delegation: the universality
`classes` the model realizes (model→class) and the `targets` it reduces to
(model→model): `(quantity, bc, classes, targets)`.
"""
function delegations(model)
    m_T = _as_type(model)
    classes = Symbol[r.class for r in REALIZES if r.model === m_T]
    targets = Type[r.target for r in REDUCES if r.source === m_T]
    return [
        (quantity=e.quantity, bc=e.bc, classes=classes, targets=targets) for
        e in REGISTRY if e.model === m_T && _is_delegation(e.method)
    ]
end

# ── implements : quantity → models ───────────────────────────────────────────
"""
    implementations_of(quantity) -> Vector{NamedTuple}

Every model with a registered row for `quantity`: `(model, bc, scheme, status)`.
"""
function implementations_of(quantity)
    q_T = _as_type(quantity)
    return [
        (model=e.model, bc=e.bc, scheme=e.scheme, status=e.status) for
        e in REGISTRY if e.quantity === q_T
    ]
end
