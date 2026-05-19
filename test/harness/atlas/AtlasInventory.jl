# test/harness/atlas/AtlasInventory.jl
# v2 framework prototype. STATIC AST scan of existing verify(...) cards.
# Hardened: non-Symbol kwarg keys are skipped; per-verify extraction is
# wrapped so a non-parseable card is tallied (PARSE_FAILS) not a crash.
module AtlasInventory

using Base.Meta: parseall

const _BS = Char(92)
const _DQ = Char(34)
const _NL = Char(10)
const PARSE_FAILS = Tuple{String,String}[]

struct Regime
    token::String
    prose::String
    arity::String
    predicate::Function
end

const TFIM_VOCAB = Regime[
    Regime("@critical", "TFIM critical point h = J", "point",
           p -> haskey(p, :h) && haskey(p, :J) && isapprox(p.h, p.J; atol=1e-12)),
    Regime("@ordered", "ordered phase h < J", "point",
           p -> haskey(p, :h) && haskey(p, :J) && p.h < p.J),
    Regime("@disordered", "disordered phase h > J", "point",
           p -> haskey(p, :h) && haskey(p, :J) && p.h > p.J),
]

const SWEEP = Regime("@sweep", "parametric sweep (curve)", "curve", _ -> false)

_headsym(ex) = ex isa Expr && ex.head === :call ? ex.args[1] :
               ex isa Symbol ? ex : nothing

function _literal_params(modelex)
    pairs = Pair{Symbol,Float64}[]
    modelex isa Expr || return NamedTuple()
    for a in modelex.args
        if a isa Expr && a.head === :parameters
            for kw in a.args
                if kw isa Expr && kw.head === :kw && kw.args[1] isa Symbol &&
                   (kw.args[2] isa Number)
                    push!(pairs, kw.args[1] => Float64(kw.args[2]))
                end
            end
        elseif a isa Expr && a.head === :kw && a.args[1] isa Symbol &&
               (a.args[2] isa Number)
            push!(pairs, a.args[1] => Float64(a.args[2]))
        end
    end
    return (; pairs...)
end

function _kwargs(callex)
    d = Dict{Symbol,Any}()
    for a in callex.args
        if a isa Expr && a.head === :parameters
            for kw in a.args
                if kw isa Expr && kw.head === :kw && kw.args[1] isa Symbol
                    d[kw.args[1]] = kw.args[2]
                end
            end
        elseif a isa Expr && a.head === :kw && a.args[1] isa Symbol
            d[a.args[1]] = a.args[2]
        end
    end
    return d
end

function _independence(route::Symbol, refs)
    structural = route in (:second_closed_form, :ed_finite_size, :literature_value)
    disc = if route === :ed_finite_size
        "ed:dense-diagonalization"
    elseif route === :literature_value
        "lit:" * (refs === nothing ? "?" : string(refs))
    elseif route === :second_closed_form
        "cf:" * (refs === nothing ? "?" : string(refs))
    else
        "asserted:" * string(route)
    end
    return (structural ? "structural" : "asserted"), disc
end

function _jstr(s)
    t = replace(string(s), string(_BS) => string(_BS, _BS))
    t = replace(t, string(_DQ) => string(_BS, _DQ))
    t = replace(t, string(_NL) => " ")
    return string(_DQ, t, _DQ)
end

_kv(k, qv) = string(_DQ, k, _DQ, ":", qv)

struct Card
    hub::String
    regime::String
    arity::String
    file::String
    testset::String
    mechanism::String
    independence::String
    discriminant::String
    refs::String
end

function _refs_text(ex)
    ex === nothing && return ""
    if ex isa Expr && ex.head === :vect
        return join((x isa String ? x : string(x) for x in ex.args), " | ")
    end
    return string(ex)
end

function _handle_verify!(out, ex, file, testset)
    pos = filter(a -> !(a isa Expr && a.head in (:parameters, :kw)), ex.args)[2:end]
    length(pos) >= 3 || return
    model, qty, bc = pos[1], pos[2], pos[3]
    kw = _kwargs(ex)
    route = get(kw, :route, nothing)
    route isa QuoteNode && (route = route.value)
    route isa Symbol || (route = :unknown)
    refs = _refs_text(get(kw, :refs, nothing))
    params = _literal_params(model)
    reg = SWEEP
    for r in TFIM_VOCAB
        try
            r.predicate(params) && (reg = r; break)
        catch
        end
    end
    ind, disc = _independence(route, get(kw, :refs, nothing))
    hub = string(_headsym(model), "/", _headsym(qty), "/", _headsym(bc))
    push!(out, Card(hub, reg.token, reg.arity, file, testset,
                    string(route), ind, disc, refs))
    return
end

function _scan_expr!(out, ex, file, testset)
    ex isa Expr || return
    if ex.head === :macrocall && ex.args[1] === Symbol("@testset")
        nm = testset
        for a in ex.args
            a isa String && (nm = a)
        end
        for a in ex.args
            _scan_expr!(out, a, file, nm)
        end
        return
    end
    if ex.head === :call && _headsym(ex) === :verify
        try
            _handle_verify!(out, ex, file, testset)
        catch err
            push!(PARSE_FAILS, (file, string(testset, " :: ", typeof(err))))
        end
        return
    end
    for a in ex.args
        _scan_expr!(out, a, file, testset)
    end
end

function scan_file(path::AbstractString)
    out = Card[]
    src = read(path, String)
    top = parseall(src; filename=path)
    parts = split(path, "QAtlas.jl/")
    rel = length(parts) > 1 ? parts[end] : path
    _scan_expr!(out, top, rel, "")
    return out
end

function scan_dir(dir::AbstractString)
    out = Card[]
    for (root, _, files) in walkdir(dir)
        for f in files
            endswith(f, ".jl") || continue
            append!(out, scan_file(joinpath(root, f)))
        end
    end
    sort!(out; by=c -> (c.hub, c.regime, c.file, c.mechanism))
    return out
end

function to_jsonl(cards::Vector{Card})
    io = IOBuffer()
    for c in cards
        println(io, string("{",
            _kv("hub", _jstr(c.hub)), ",",
            _kv("regime", _jstr(c.regime)), ",",
            _kv("arity", _jstr(c.arity)), ",",
            _kv("plane", _jstr("why")), ",",
            _kv("mechanism", _jstr(c.mechanism)), ",",
            _kv("independence", _jstr(c.independence)), ",",
            _kv("discriminant", _jstr(c.discriminant)), ",",
            _kv("file", _jstr(c.file)), ",",
            _kv("testset", _jstr(c.testset)), ",",
            _kv("refs", _jstr(c.refs)), "}"))
    end
    return String(take!(io))
end

write_inventory(path, cards) = write(path, to_jsonl(cards))

end # module
