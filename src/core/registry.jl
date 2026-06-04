# core/registry.jl — declarative implementation registry.
#
# Each (model, quantity, bc) triple that has a directly-implemented
# `fetch` method also gets a one-liner `@register` declaration next to
# it, capturing the human metadata that `methods(fetch)` cannot:
#   * `method`      — algorithm tag (`:bdg`, `:dense_ed`, `:analytic`,
#                      `:transfer_matrix`, `:bethe_ansatz`, `:tba`, `:pfaffian`,
#                      `:not_implemented`, …)
#   * `status`      — claim kind, orthogonal to `reliability`: `:exact`
#                      (analytic closed form), `:bound` (one-sided
#                      inequality / saturating universal constant), or
#                      `:approx` (domain-limited approximation, e.g. a
#                      high-T expansion).  See [`STATUS_VALUES`](@ref).
#   * `reliability` — `:high` (closed-form + literature-tested),
#                      `:medium` (ED only / cross-check),
#                      `:low` (heuristic, not validated),
#                      `:not_implemented`.  Aligned with the
#                      (a)/(b)/(c) test categories in #118.
#   * `tested_in`   — relative path to the test file that validates
#                      this triple (or `nothing` if no dedicated test).
#   * `references`  — short literature pointers (author + year).
#   * `notes`       — caller-facing caveats (granularity, kwargs, …).
#
# `implementation_status()` returns the registry as
# `Vector{NamedTuple}`, which is `Tables.jl`-compatible without us
# taking a Tables dependency — downstream users wanting a `DataFrame`
# can call `DataFrame(implementation_status())` themselves.
#
# Conversion fallbacks (e.g. `Energy(:per_site)` at OBC routed through
# the `Energy(:total)` native + `÷ N`) are *not* registered separately:
# the registry reflects native implementations, and the routing is
# automatic by design.

"""
    Implementation

A single `(model, quantity, bc)` row of the QAtlas implementation
registry.  See [`@register`](@ref) for how rows are added and
[`implementation_status`](@ref) for how to query them.
"""
struct Implementation
    model::Type
    quantity::Type
    bc::Type
    method::Symbol
    status::Symbol
    direction::Union{Symbol,Nothing}
    reliability::Symbol
    tested_in::Union{String,Nothing}
    references::Vector{String}
    notes::String
end

"""
    STATUS_VALUES

The controlled vocabulary for the `status` axis of a registered
implementation — *what kind of mathematical claim* the row makes. This
is orthogonal to `reliability` (how confident the implementation is) and
to the test-corroboration level tracked by the atlas harness:

  * `:exact`  — analytic closed form; verified as an equality against the
                 literature value (the historical default; every legacy
                 row is `:exact`).
  * `:bound`  — a one-sided inequality. Either a *saturating* universal
                 constant (equality at the optimal state, e.g. a Tsirelson
                 bound) or a *variational* bound (an independently measured
                 quantity stays ≤/≥ the fetched value, with no saturation
                 guaranteed). The ≤/≥ direction lives on the verification
                 card, not here.
  * `:approx` — a domain-limited approximation (e.g. a high-temperature
                 expansion): correct on a stated region of validity with a
                 known leading error order.

`register!` rejects any `status` outside this tuple, so a typo fails at
package load time rather than silently mislabelling a claim.
"""
const STATUS_VALUES = (:exact, :bound, :approx)

"""
    BOUND_DIRECTIONS

The controlled vocabulary for the `direction` of a `status=:bound` row —
*which side* of the bounded quantity the fetched value constrains:

  * `:upper` — the fetched value is an upper bound; an independent witness
               stays `≤` it (verified with `verify_bound(...; relation=:leq)`).
  * `:lower` — the fetched value is a lower bound; a witness stays `≥` it
               (`relation=:geq`).

A bound is fully pinned by *what* it bounds (the registry `quantity`),
*which way* (`direction`), and *whose* bound it is (`references`, plus a
`source=` selector when several bounds share one quantity — e.g. the
classical / quantum / no-signalling CHSH bounds). Non-bound rows carry
`direction === nothing`; `register!` enforces both halves.
"""
const BOUND_DIRECTIONS = (:upper, :lower)

"""
    REGISTRY :: Vector{Implementation}

Module-level mutable vector populated at include-time by `@register`
calls scattered across `src/models/.../<Model>_registry.jl` files.
Public read API: [`implementation_status`](@ref).
"""
const REGISTRY = Implementation[]

"""
    register!(model_T, quantity_T, bc_T;
              method=:unknown, status=:exact, direction=nothing,
              reliability=:unknown, tested_in=nothing,
              references=String[], notes="")

Push a new [`Implementation`](@ref) row into [`REGISTRY`](@ref).
Usually called via the [`@register`](@ref) macro for ergonomics.
`status=:bound` requires a `direction` (`:upper`/`:lower`); any other
status must leave `direction === nothing`. See [`BOUND_DIRECTIONS`](@ref).
"""
function register!(
    model_T::Type,
    quantity_T::Type,
    bc_T::Type;
    method::Symbol=:unknown,
    status::Symbol=:exact,
    direction::Union{Symbol,Nothing}=nothing,
    reliability::Symbol=:unknown,
    tested_in::Union{String,Nothing}=nothing,
    references::AbstractVector{<:AbstractString}=String[],
    notes::AbstractString="",
)
    status in STATUS_VALUES || throw(
        ArgumentError("register!: status must be one of $(STATUS_VALUES); got :$(status)"),
    )
    # A bound must pin which way it constrains; conversely a non-bound must
    # not carry a direction (meaningless next to an equality/approximation).
    if status === :bound
        direction in BOUND_DIRECTIONS || throw(
            ArgumentError(
                "register!: status=:bound requires direction ∈ $(BOUND_DIRECTIONS) " *
                "(an upper/lower bound must say which way it constrains); got " *
                "$(repr(direction))",
            ),
        )
    else
        direction === nothing || throw(
            ArgumentError(
                "register!: direction is only meaningful for status=:bound; got " *
                "status=:$(status) with direction=$(repr(direction))",
            ),
        )
    end
    push!(
        REGISTRY,
        Implementation(
            model_T,
            quantity_T,
            bc_T,
            method,
            status,
            direction,
            reliability,
            tested_in,
            String[r for r in references],
            String(notes),
        ),
    )
    return nothing
end

"""
    @register Model Quantity BC method=… reliability=… tested_in=… references=… notes=…

Thin macro around [`register!`](@ref).  Lets each model file register
its native fetch methods declaratively, e.g.

```julia
@register TFIM Energy{:total} OBC method=:bdg reliability=:high \\
    tested_in="test/models/test_TFIM_thermal.jl" \\
    references=["Pfeuty 1970"]
```

The three positional arguments are spliced as types; the remaining
`key=value` pairs are forwarded as keyword arguments to
[`register!`](@ref).
"""
macro register(model_T, quantity_T, bc_T, kwargs...)
    kw_exprs = map(kwargs) do kw
        kw isa Expr && kw.head === :(=) || error("@register: expected key=value, got $kw")
        return Expr(:kw, kw.args[1], esc(kw.args[2]))
    end
    return Expr(
        :call,
        register!,
        Expr(:parameters, kw_exprs...),
        esc(model_T),
        esc(quantity_T),
        esc(bc_T),
    )
end

# ──────────────────────────────────────────────────────────────────────
# Smoke-check helper: distinguish "registered + fetch method exists" from
# "registered but only the catch-all error method matches".  The
# catch-all lives in core/type.jl; we capture it here at registry load
# time so future fetch additions can't accidentally re-shadow it.
# ──────────────────────────────────────────────────────────────────────

const _CATCH_ALL_FETCH_METHOD = which(
    fetch, Tuple{AbstractQAtlasModel,AbstractQuantity,BoundaryCondition}
)

"""
    has_native_fetch(impl::Implementation) -> Bool

`true` iff `which(fetch, (impl.model, impl.quantity, impl.bc))` resolves
to a method *more specific* than the catch-all in `core/type.jl`.
Conversion fallbacks (e.g. the generic `Energy{:per_site}` ↔
`Energy{:total}` router) count as "native" because they are a real
dispatchable implementation — they just live above the model layer.

Used by `test/core/test_registry.jl` to detect registry rows that
silently lost their backing fetch method.
"""
function has_native_fetch(impl::Implementation)
    m = which(fetch, Tuple{impl.model,impl.quantity,impl.bc})
    return m !== _CATCH_ALL_FETCH_METHOD
end

# ──────────────────────────────────────────────────────────────────────
# Query API
# ──────────────────────────────────────────────────────────────────────

function _to_nt(e::Implementation)
    (
        model=e.model,
        quantity=e.quantity,
        bc=e.bc,
        method=e.method,
        status=e.status,
        reliability=e.reliability,
        tested_in=e.tested_in,
        references=e.references,
        notes=e.notes,
    )
end

"""
    implementation_status() -> Vector{NamedTuple}
    implementation_status(model::AbstractQAtlasModel)
    implementation_status(::Type{<:AbstractQAtlasModel})
    implementation_status(quantity::AbstractQuantity)
    implementation_status(::Type{<:AbstractQuantity})
    implementation_status(queue::AbstractVector)

Return registry rows as `NamedTuple`s (`Tables.jl`-compatible without a
Tables dependency).

- No-arg: every registered triple.
- `model` / `quantity` (instance or type): rows whose corresponding type
  field matches exactly (no subtype walking — model parameters are part
  of the identity here).
- `queue`: a vector of `(model, quantity, bc)` triples (each component
  may be either an instance or a type).  Returns one row per queue
  entry that is registered, dropping entries that are not.

Use this to plan downstream work — e.g. before writing tests for a new
ThermalMPS workload, query the queue you intend to validate against.
"""
implementation_status() = [_to_nt(e) for e in REGISTRY]

function implementation_status(::Type{M}) where {M<:AbstractQAtlasModel}
    [_to_nt(e) for e in REGISTRY if e.model === M]
end
implementation_status(model::AbstractQAtlasModel) = implementation_status(typeof(model))

function implementation_status(::Type{Q}) where {Q<:AbstractQuantity}
    [_to_nt(e) for e in REGISTRY if e.quantity === Q]
end
implementation_status(quantity::AbstractQuantity) = implementation_status(typeof(quantity))

function implementation_status(queue::AbstractVector)
    out = NamedTuple[]
    for q in queue
        length(q) == 3 || error(
            "implementation_status(queue): each element must be a " *
            "(model, quantity, bc) triple; got length $(length(q))",
        )
        m_T = q[1] isa Type ? q[1] : typeof(q[1])
        q_T = q[2] isa Type ? q[2] : typeof(q[2])
        bc_T = q[3] isa Type ? q[3] : typeof(q[3])
        for e in REGISTRY
            if e.model === m_T && e.quantity === q_T && e.bc === bc_T
                push!(out, _to_nt(e))
                break
            end
        end
    end
    return out
end

# Coerce an instance-or-type argument to its `Type` for registry matching.
_as_type(x) = x isa Type ? x : typeof(x)

"""
    references_for(model, quantity, bc) -> Vector{String}
    references_for(model, quantity)     -> Vector{String}
    references_for(model)               -> Vector{String}

Return the literature references — `references.bib` bibkeys — that the
registered implementation(s) for the given arguments rest on. Use it as a
companion to [`fetch`](@ref): pass the same `model` (and optionally
`quantity` and boundary condition `bc`) you are calling to see *which
papers that closed-form value derives from and is checked against*.

These are exactly the keys rendered on the model's documentation page and
in the global [Reference List](@ref); resolve a key to its full entry
there (or in
[`docs/references.bib`](https://github.com/sotashimozono/QAtlas.jl/blob/main/docs/references.bib)).

Arguments may be instances or types, mirroring
[`implementation_status`](@ref). With fewer arguments the references are
aggregated over the unspecified axes (all boundary conditions for a
`(model, quantity)` pair; every registered quantity for a bare `model`).
The result is de-duplicated and sorted, and is empty when no matching
registry row carries references (including when the triple is not
registered at all).

# Examples
```julia
julia> references_for(TFIM(), Energy{:per_site}(), Infinite())
1-element Vector{String}:
 "Pfeuty1970"

julia> references_for(TFIM())            # every paper TFIM rests on
…
```
"""
function references_for(model, quantity, bc)
    m_T, q_T, bc_T = _as_type(model), _as_type(quantity), _as_type(bc)
    refs = String[]
    for e in REGISTRY
        if e.model === m_T && e.quantity === q_T && e.bc === bc_T
            append!(refs, e.references)
        end
    end
    return sort!(unique!(refs))
end

function references_for(model, quantity)
    m_T, q_T = _as_type(model), _as_type(quantity)
    refs = String[]
    for e in REGISTRY
        if e.model === m_T && e.quantity === q_T
            append!(refs, e.references)
        end
    end
    return sort!(unique!(refs))
end

function references_for(model)
    m_T = _as_type(model)
    refs = String[]
    for e in REGISTRY
        if e.model === m_T
            append!(refs, e.references)
        end
    end
    return sort!(unique!(refs))
end

# ──────────────────────────────────────────────────────────────────────
# Markdown rendering
# ──────────────────────────────────────────────────────────────────────

# Strip the leading `QAtlas.` (and any submodule prefix) from a `Type`
# so the rendered table reads `TFIM`, `Energy{:total}`, `OBC` rather
# than `QAtlas.TFIM` etc.
function _short_type(T::Type)
    s = string(T)
    return replace(s, r"^QAtlas\.|Main\.QAtlas\." => "")
end

"""
    implementation_status_markdown([io::IO=stdout], entries=implementation_status())

Render `entries` (any iterable of `NamedTuple` rows from
[`implementation_status`](@ref)) as a GitHub-flavoured Markdown table
to `io`.
"""
function implementation_status_markdown(io::IO=stdout, entries=implementation_status())
    println(
        io,
        "| Model | Quantity | BC | Method | Status | Reliability | Tested in | References |",
    )
    println(io, "|---|---|---|---|---|---|---|---|")
    for e in entries
        println(
            io,
            "| ",
            _short_type(e.model),
            " | ",
            _short_type(e.quantity),
            " | ",
            _short_type(e.bc),
            " | `",
            e.method,
            "`",
            " | `",
            e.status,
            "`",
            " | `",
            e.reliability,
            "`",
            " | ",
            something(e.tested_in, "—"),
            " | ",
            isempty(e.references) ? "—" : join(e.references, "; "),
            " |",
        )
    end
    return nothing
end
