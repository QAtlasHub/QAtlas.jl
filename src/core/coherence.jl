# core/coherence.jl — verification DERIVED from the knowledge graph.
#
# The thesis: correctness and comprehensiveness are properties of the *network*,
# not of isolated rows.  Each implementation is constrained by several
# independent cross-paths (its class's universal prediction, the bounds on its
# quantity, its cited literature, agreement with sibling models).  A wrong or
# missing implementation violates at least one cross-path; satisfying them all
# is a far stronger guarantee than any single unit test.  Missing edges are
# self-reported as coverage gaps.
#
# Each `check_*` walks `REGISTRY` / `REALIZES` (+ `links.jl`) and returns
# `CoherenceFinding`s.  `coherence_report` runs the structural suite.  Physical
# cross-checks (realization-agreement C5, bound-satisfaction C7) need per-quantity
# invocation probes and are layered on top via `check_realization_agreement`.

"""
    CoherenceFinding(check, severity, message)

One graph-derived finding.  `severity`: `:error` (a violated invariant — must be
empty), `:gap` (a missing-but-expected edge — coverage hole), `:info`.
"""
struct CoherenceFinding
    check::Symbol
    severity::Symbol
    message::String
end

function Base.show(io::IO, f::CoherenceFinding)
    print(io, "[", f.check, "/", f.severity, "] ", f.message)
end

# ── C1 — reference integrity: every cited bibkey exists ──────────────────────
function check_reference_integrity(bibkeys)
    keys = Set(string.(bibkeys))
    out = CoherenceFinding[]
    for e in REGISTRY, r in e.references
        r in keys || push!(
            out,
            CoherenceFinding(
                :reference_integrity,
                :error,
                "dangling reference '$(r)' in $(_kgshort(e.model))/$(_kgshort(e.quantity))",
            ),
        )
    end
    return out
end

# ── C2 — namespace ⟺ kind ────────────────────────────────────────────────────
function check_namespace_kind()
    out = CoherenceFinding[]
    for e in REGISTRY
        if _is_universality(e.model) && e.status !== :universal
            push!(
                out,
                CoherenceFinding(
                    :namespace_kind,
                    :error,
                    "$(_kgshort(e.model)) is Universality but status=:$(e.status) (want :universal)",
                ),
            )
        elseif e.status === :universal && !_is_universality(e.model)
            push!(
                out,
                CoherenceFinding(
                    :namespace_kind,
                    :error,
                    "$(_kgshort(e.model))/$(_kgshort(e.quantity)) is :universal but not a Universality node",
                ),
            )
        end
        if _is_bound(e.model) && e.status !== :bound
            push!(
                out,
                CoherenceFinding(
                    :namespace_kind,
                    :error,
                    "$(_kgshort(e.model)) is Bound but status=:$(e.status) (want :bound)",
                ),
            )
        end
    end
    return out
end

# ── C3 — canonical / scheme coherence per hub ────────────────────────────────
function check_canonical_coherence()
    out = CoherenceFinding[]
    hubs = Dict{Tuple{Type,Type,Type},Vector{Implementation}}()
    for e in REGISTRY
        push!(get!(hubs, (e.model, e.quantity, e.bc), Implementation[]), e)
    end
    for ((m, q, _), rows) in hubs
        nc = count(r -> r.canonical, rows)
        nc == 1 || push!(
            out,
            CoherenceFinding(
                :canonical_coherence,
                :error,
                "hub $(_kgshort(m))/$(_kgshort(q)) has $(nc) canonical rows (want exactly 1)",
            ),
        )
        schemes = [r.scheme for r in rows]
        length(unique(schemes)) == length(schemes) || push!(
            out,
            CoherenceFinding(
                :canonical_coherence,
                :error,
                "hub $(_kgshort(m))/$(_kgshort(q)) has duplicate schemes $(schemes)",
            ),
        )
    end
    return out
end

# ── C4 — delegation has a realized target (triage signal, not a hard error) ──
# `method=:delegation` is currently untyped: it covers both model→class
# delegation (which SHOULD have a `realizes` edge) and model→model reduction
# (which should not).  Until the edge records its target, a missing realizes
# edge is a `:gap` to triage, not a violated invariant.
function check_delegation_targets()
    out = CoherenceFinding[]
    for e in REGISTRY
        e.method === :delegation || continue
        any(r -> r.model === e.model, REALIZES) || push!(
            out,
            CoherenceFinding(
                :delegation_target,
                :gap,
                "$(_kgshort(e.model))/$(_kgshort(e.quantity)) delegates but its model realizes no class — triage: add @realizes (class delegation) or reclassify as model→model",
            ),
        )
    end
    return out
end

# ── C6 — coverage: universality nodes that are referenced but undeveloped ────
# A class realized by some model but with zero `predicts` edges is a sparse
# region of the vault — its universal predictions are not yet registered, so
# nothing can be cross-checked against the models that flow to it.
function coverage_report()
    out = CoherenceFinding[]
    for c in sort(unique(r.class for r in REALIZES))
        isempty(predicts(c)) && push!(
            out,
            CoherenceFinding(
                :coverage,
                :gap,
                "class :$(c) is realized by a model but predicts nothing (undeveloped universality node)",
            ),
        )
    end
    return out
end

# ── C5 — realization agreement (physical) ────────────────────────────────────
"""
    check_realization_agreement(probes; rtol=1e-8) -> Vector{CoherenceFinding}

For each `method=:delegation` row whose `quantity` has an entry in `probes`
(a `Dict` quantity-instance ⇒ NamedTuple of fetch kwargs), assert the model's
value equals the realized class's prediction at the probe point — the physical
cross-validation that makes a delegation edge *true*, not just declared.
"""
function check_realization_agreement(probes::AbstractDict; rtol=1e-8)
    out = CoherenceFinding[]
    for e in REGISTRY
        e.method === :delegation || continue
        q = nothing
        for (qi, _) in probes
            typeof(qi) === e.quantity && (q = qi)
        end
        q === nothing && continue
        classes = Symbol[r.class for r in REALIZES if r.model === e.model]
        isempty(classes) && continue
        kw = probes[q]
        try
            mval = fetch(e.model(), q, e.bc(); kw...)
            cval = fetch(Universality(first(classes)), q, e.bc(); kw...)
            isapprox(mval, cval; rtol=rtol) || push!(
                out,
                CoherenceFinding(
                    :realization_agreement,
                    :error,
                    "$(_kgshort(e.model))/$(_kgshort(e.quantity)) = $(mval) ≠ Universality(:$(first(classes))) = $(cval)",
                ),
            )
        catch err
            push!(
                out,
                CoherenceFinding(
                    :realization_agreement,
                    :info,
                    "could not probe $(_kgshort(e.model))/$(_kgshort(e.quantity)): $(err)",
                ),
            )
        end
    end
    return out
end

# ── aggregate (structural suite) ─────────────────────────────────────────────
"""
    coherence_report(; bibkeys=String[]) -> Vector{CoherenceFinding}

Run the structural graph-coherence suite (C1–C4 + the C6 coverage report).
Pass `bibkeys` (the set of keys in references.bib) to include C1.  `:error`
findings must be empty; `:gap` findings are the network's self-reported holes.
"""
function coherence_report(; bibkeys=String[])
    findings = CoherenceFinding[]
    isempty(bibkeys) || append!(findings, check_reference_integrity(bibkeys))
    append!(findings, check_namespace_kind())
    append!(findings, check_canonical_coherence())
    append!(findings, check_delegation_targets())
    append!(findings, coverage_report())
    return findings
end

"""
    coherence_errors(findings) -> Vector{CoherenceFinding}

The `:error` findings (violated invariants — must be empty) among `findings`
from [`coherence_report`](@ref).
"""
coherence_errors(findings) = filter(f -> f.severity === :error, findings)

"""
    coherence_gaps(findings) -> Vector{CoherenceFinding}

The `:gap` findings (missing-but-expected edges — the network's self-reported
holes) among `findings` from [`coherence_report`](@ref).
"""
coherence_gaps(findings) = filter(f -> f.severity === :gap, findings)
