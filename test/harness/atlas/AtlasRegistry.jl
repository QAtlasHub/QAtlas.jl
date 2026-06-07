# test/harness/atlas/AtlasRegistry.jl
#
# Static parse of src/models/.../$(MODEL)_registry.jl @register(...) calls
# -> the set of (model,quantity,bc) hubs src CLAIMS to solve. Used by the
# docs-as-view generator for the leverage risk-linter anti-join (M3):
#   src claims hub  AND  count(structural inventory cards) == 0  -> flag.
# Pure Base.Meta (no QAtlas load).
module AtlasRegistry

using Base.Meta: parseall

struct Claim
    hub::String
    model::String
    quantity::String
    bc::String
    method::String
    status::String
    reliability::String
    refs::String
    notes::String
end

function _sym(ex)
    if ex isa Symbol
        string(ex)
    else
        (
            if ex isa Expr && ex.head === :curly
                string(ex.args[1])
            elseif ex isa Expr && ex.head === :call
                string(ex.args[1])
            else
                string(ex)
            end
        )
    end
end

function _regkw(args)
    d = Dict{Symbol,Any}()
    for a in args
        if a isa Expr && a.head in (:kw, :(=)) && a.args[1] isa Symbol
            v = a.args[2]
            d[a.args[1]] = v isa QuoteNode ? v.value : v
        end
    end
    return d
end

include(joinpath(@__DIR__, "_atlas_common.jl"))   # _refs_text (shared)

function _walk!(out, ex)
    ex isa Expr || return nothing
    if ex.head === :macrocall && ex.args[1] === Symbol("@register")
        pos = filter(
            a -> !(a isa LineNumberNode) && !(a isa Expr && a.head in (:kw, :(=))), ex.args
        )[2:end]
        if length(pos) >= 3
            model, qty, bc = _sym(pos[1]), _sym(pos[2]), _sym(pos[3])
            kw = _regkw(ex.args)
            push!(
                out,
                Claim(
                    string(model, "/", qty, "/", bc),
                    model,
                    qty,
                    bc,
                    string(get(kw, :method, "")),
                    string(get(kw, :status, "")),
                    string(get(kw, :reliability, "")),
                    _refs_text(get(kw, :references, nothing)),
                    string(get(kw, :notes, "")),
                ),
            )
        end
        return nothing
    end
    for a in ex.args
        _walk!(out, a)
    end
end

function scan_registry(path::AbstractString)
    out = Claim[]
    _walk!(out, parseall(read(path, String); filename=path))
    return out
end

# ── @realizes (backend register): model ↔ universality-class membership ──
# `@register` wires the frontend (the hub/atlas the user sees); `@realizes`
# is the backend config recording which universality class a model flows to,
# and the regime where it does.  Statically parsed, same as @register.
struct Realization
    model::String
    class::String
    regime::String
    refs::String
end

function _walk_realizes!(out, ex)
    ex isa Expr || return nothing
    if ex.head === :macrocall && ex.args[1] === Symbol("@realizes")
        pos = filter(
            a -> !(a isa LineNumberNode) && !(a isa Expr && a.head in (:kw, :(=))), ex.args
        )[2:end]
        if length(pos) >= 2
            model = _sym(pos[1])
            class = pos[2] isa QuoteNode ? string(pos[2].value) : _sym(pos[2])
            kw = _regkw(ex.args)
            push!(
                out,
                Realization(
                    model,
                    class,
                    string(get(kw, :regime, "")),
                    _refs_text(get(kw, :references, nothing)),
                ),
            )
        end
        return nothing
    end
    for a in ex.args
        _walk_realizes!(out, a)
    end
end

function scan_realizes(path::AbstractString)
    out = Realization[]
    isfile(path) || return out
    _walk_realizes!(out, parseall(read(path, String); filename=path))
    return out
end

# ── @about (model description card): summary + Hamiltonian (LaTeX) ───────
# Authored with `raw"…"` so LaTeX backslashes survive; in the AST a `raw"…"`
# literal is an `@raw_str` macrocall, so unwrap it back to its String.
struct About
    model::String
    summary::String
    hamiltonian::String
    refs::String
end

function _strval(v)
    v isa String && return v
    if v isa Expr && v.head === :macrocall && v.args[1] === Symbol("@raw_str")
        return String(v.args[end])
    end
    if v isa Expr && v.head === :string   # interpolation/concat: keep literal parts
        return join((p isa String ? p : "" for p in v.args))
    end
    return string(v)
end

function _walk_about!(out, ex)
    ex isa Expr || return nothing
    if ex.head === :macrocall && ex.args[1] === Symbol("@about")
        pos = filter(
            a -> !(a isa LineNumberNode) && !(a isa Expr && a.head in (:kw, :(=))), ex.args
        )[2:end]
        if length(pos) >= 1
            model = _sym(pos[1])
            kw = _regkw(ex.args)
            push!(
                out,
                About(
                    model,
                    _strval(get(kw, :summary, "")),
                    _strval(get(kw, :hamiltonian, "")),
                    _refs_text(get(kw, :references, nothing)),
                ),
            )
        end
        return nothing
    end
    for a in ex.args
        _walk_about!(out, a)
    end
end

function scan_about(path::AbstractString)
    out = About[]
    isfile(path) || return out
    _walk_about!(out, parseall(read(path, String); filename=path))
    return out
end

end # module
