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

end # module
