# test/harness/atlas/AtlasInventory.jl
# v2 framework. STATIC AST scan of existing verify(...) cards.
# Hardened: non-Symbol kwarg keys are skipped; per-verify extraction is
# wrapped so a non-parseable card is tallied (PARSE_FAILS) not a crash.
# Step 3: per-model regime vocabulary (MODEL_VOCAB) — literal-param point
# cards get a named physical regime; loop-variable sweeps stay @sweep
# (a curve, the honest graceful behaviour).
module AtlasInventory

using Base.Meta: parseall

const _BS = Char(92)
const _DQ = Char(34)
const _NL = Char(10)
const PARSE_FAILS = Tuple{String,String}[]
const _ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

struct Regime
    token::String
    prose::String
    arity::String
    predicate::Function
end

const SWEEP = Regime("@sweep", "parametric sweep (curve)", "curve", _ -> false)

# Defensive numeric accessor: missing key -> NaN so predicates on
# loop-variable cards (empty params) simply fail -> @sweep.
_pv(p, k) = haskey(p, k) ? float(getfield(p, k)) : NaN
_eq(a, b) = isfinite(a) && isapprox(a, b; atol=1e-12)

# Per-model controlled regime vocabulary. First matching predicate wins.
# A `_ -> true` final entry is a model-level constant regime (the model
# sits at one physical point regardless of scale params).
const MODEL_VOCAB = Dict{String,Vector{Regime}}(
    "TFIM" => Regime[
        Regime(
            "@critical",
            "TFIM quantum critical point h = J",
            "point",
            p -> _eq(_pv(p, :h), _pv(p, :J)),
        ),
        Regime(
            "@ordered",
            "ferromagnetic ordered phase h < J",
            "point",
            p -> _pv(p, :h) < _pv(p, :J),
        ),
        Regime(
            "@disordered",
            "paramagnetic disordered phase h > J",
            "point",
            p -> _pv(p, :h) > _pv(p, :J),
        ),
    ],
    "XXZ1D" => Regime[
        Regime(
            "@free_fermion",
            "XX free-fermion point Δ = 0",
            "point",
            p -> _eq(_pv(p, :Δ), 0.0),
        ),
        Regime(
            "@su2",
            "isotropic Heisenberg point Δ = 1",
            "point",
            p -> _eq(_pv(p, :Δ), 1.0),
        ),
        Regime("@fm", "ferromagnetic point Δ = -1", "point", p -> _eq(_pv(p, :Δ), -1.0)),
        Regime(
            "@gapless",
            "critical Luttinger liquid |Δ| < 1",
            "curve",
            p -> (d=_pv(p, :Δ); isfinite(d) && -1 < d < 1),
        ),
        Regime(
            "@gapped",
            "gapped regime |Δ| > 1",
            "curve",
            p -> (d=_pv(p, :Δ); isfinite(d) && abs(d) > 1),
        ),
    ],
    "Heisenberg1D" =>
        Regime[Regime("@su2", "SU(2) isotropic Heisenberg chain", "point", _ -> true)],
    "HeisenbergXYZ" => Regime[
        Regime(
            "@isotropic",
            "isotropic point Jx = Jy = Jz",
            "point",
            p -> _eq(_pv(p, :Jx), _pv(p, :Jy)) && _eq(_pv(p, :Jy), _pv(p, :Jz)),
        ),
        Regime(
            "@xx",
            "XX line Jz = 0, Jx = Jy",
            "point",
            p -> _eq(_pv(p, :Jz), 0.0) && _eq(_pv(p, :Jx), _pv(p, :Jy)),
        ),
        Regime(
            "@xxz",
            "uniaxial XXZ line Jx = Jy",
            "curve",
            p -> _eq(_pv(p, :Jx), _pv(p, :Jy)),
        ),
    ],
    "S1Heisenberg1D" =>
        Regime[Regime("@haldane", "spin-1 Haldane phase", "point", _ -> true)],
    "S1XXZ1D" => Regime[Regime("@haldane", "spin-1 Haldane (Δ=1)", "point", _ -> true)],
    "S1AnisotropicD1D" =>
        Regime[Regime("@haldane", "spin-1 Haldane (D=0)", "point", _ -> true)],
    "MajumdarGhosh" =>
        Regime[Regime("@dimer", "exact orthogonal-dimer point", "point", _ -> true)],
    "Cluster1D" =>
        Regime[Regime("@cluster", "cluster-state stabiliser point", "point", _ -> true)],
    "Compass1D" => Regime[
        Regime(
            "@isotropic",
            "isotropic compass J_x = J_y",
            "point",
            p -> _eq(_pv(p, :J_x), _pv(p, :J_y)),
        ),
        Regime("@anisotropic", "anisotropic compass J_x ≠ J_y", "curve", _ -> true),
    ],
    "Kitaev1D" => Regime[
        Regime(
            "@critical",
            "topological transition |μ| = 2|t|",
            "point",
            p -> _eq(abs(_pv(p, :μ)), 2 * abs(_pv(p, :t))),
        ),
        Regime(
            "@topological",
            "topological phase |μ| < 2|t|",
            "curve",
            p -> (
                m=_pv(p, :μ);
                t=_pv(p, :t);
                isfinite(m) && isfinite(t) && abs(m) < 2 * abs(t)
            ),
        ),
        Regime(
            "@trivial",
            "trivial phase |μ| > 2|t|",
            "curve",
            p -> (
                m=_pv(p, :μ);
                t=_pv(p, :t);
                isfinite(m) && isfinite(t) && abs(m) > 2 * abs(t)
            ),
        ),
    ],
    "SchwingerModel" => Regime[
        Regime("@massless", "massless Schwinger m = 0", "point", p -> _eq(_pv(p, :m), 0.0)),
        Regime(
            "@massive",
            "massive Schwinger m ≠ 0",
            "curve",
            p -> (m=_pv(p, :m); isfinite(m) && m != 0),
        ),
    ],
    "TightBinding1D" => Regime[
        Regime("@half_filling", "half filling μ = 0", "point", p -> _eq(_pv(p, :μ), 0.0)),
        Regime(
            "@band_insulator",
            "band insulator |μ| > 2t",
            "curve",
            p -> (
                m=_pv(p, :μ);
                t=_pv(p, :t);
                isfinite(m) && isfinite(t) && abs(m) > 2 * abs(t)
            ),
        ),
    ],
    "TightBindingV1D" => Regime[
        Regime(
            "@half_filling",
            "V=0 half filling μ = 0",
            "point",
            p -> _eq(_pv(p, :μ), 0.0),
        ),
        Regime(
            "@band_insulator",
            "V=0 band insulator |μ| > 2t",
            "curve",
            p -> (
                m=_pv(p, :μ);
                t=_pv(p, :t);
                isfinite(m) && isfinite(t) && abs(m) > 2 * abs(t)
            ),
        ),
    ],
    "XYh1D" => Regime[
        Regime("@xx", "XX limit h = 0", "point", p -> _eq(_pv(p, :h), 0.0)),
        Regime(
            "@polarized",
            "polarized |h| > 2J",
            "curve",
            p -> (h=_pv(p, :h); isfinite(h) && abs(h) > 2),
        ),
    ],
    "IsingSquare" =>
        Regime[Regime("@onsager", "2D Ising Onsager critical point", "point", _ -> true)],
    "IsingChain1D" =>
        Regime[Regime("@ising1d", "1D Ising (no finite-T order)", "point", _ -> true)],
    "IsingTriangular" => Regime[Regime(
        "@triangular", "triangular Ising (frustrated AFM)", "point", _ -> true
    )],
    "CurieWeissIsing" => Regime[Regime(
        "@mean_field", "mean-field complete-graph Ising", "point", _ -> true
    )],
)

_headsym(ex) =
    if ex isa Expr && ex.head === :call
        ex.args[1]
    elseif ex isa Symbol
        ex
    else
        nothing
    end

function _literal_params(modelex)
    pairs = Pair{Symbol,Float64}[]
    modelex isa Expr || return NamedTuple()
    for a in modelex.args
        if a isa Expr && a.head === :parameters
            for kw in a.args
                if kw isa Expr &&
                    kw.head === :kw &&
                    kw.args[1] isa Symbol &&
                    (kw.args[2] isa Number)
                    push!(pairs, kw.args[1] => Float64(kw.args[2]))
                end
            end
        elseif a isa Expr &&
            a.head === :kw &&
            a.args[1] isa Symbol &&
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

# NOTE: `refs` here is the raw parsed AST node (an Expr or `nothing`),
# so the lit:/cf: discriminant is the Julia-printed AST form. This
# INTENTIONALLY differs from verify.jl `_v2_independence`, which sees
# the runtime Vector{String} and pipe-joins it. The two artifacts
# (static INVENTORY vs runtime evidence JSONL) are not byte-joinable on
# `discriminant` — do not build a cross-join assuming they match.
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

# ── R1 assurance taxonomy — single source of truth, typed ────────────
# A mistyped level is now a compile/lookup error, not a silently
# mis-bucketed hub. Used by docs/atlas/generate.jl and unit-tested by
# test/harness/atlas/test_atlas_logic.jl.
@enum AssuranceLevel begin
    UNIVERSALITY_CORROBORATED
    CORROBORATED_AT_P
    COHERENT
    CITED_ONLY
    UNCORROBORATED_BUT_FEASIBLE
end

# Models where dense ED at a physically meaningful size is infeasible:
# the published / DMRG value is the ceiling, so a missing in-repo ED
# card is the honest frontier (cited-only), NOT an actionable gap.
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

const MECH_UNIV = Set(["universality_consistency"])
const MECH_EDP = Set(["ed_finite_size", "second_closed_form"])
const MECH_COH = Set([
    "delegation_invariant", "limiting_case", "sum_rule", "retype_formula", "unknown"
])
const MECH_CITED = Set(["literature_value"])

# Highest achieved tier wins. Pure: (card mechanisms, model infeasible?).
function assurance_level(mechs::AbstractSet, model_ed_infeasible::Bool)::AssuranceLevel
    isempty(intersect(mechs, MECH_UNIV)) || return UNIVERSALITY_CORROBORATED
    isempty(intersect(mechs, MECH_EDP)) || return CORROBORATED_AT_P
    isempty(intersect(mechs, MECH_COH)) || return COHERENT
    isempty(intersect(mechs, MECH_CITED)) || return CITED_ONLY
    model_ed_infeasible && return CITED_ONLY
    return UNCORROBORATED_BUT_FEASIBLE
end

const _LEVEL_DISPLAY = Dict{AssuranceLevel,NTuple{3,String}}(
    UNIVERSALITY_CORROBORATED => ("universality-corroborated", "🟣", "tip"),
    CORROBORATED_AT_P => ("corroborated-at-p", "🟢", "tip"),
    COHERENT => ("coherent", "🔵", "note"),
    CITED_ONLY => ("cited-only", "⚪", "note"),
    UNCORROBORATED_BUT_FEASIBLE => ("uncorroborated-but-feasible", "🟠", "warning"),
)
# (display name, badge emoji, Documenter admonition) for a level.
level_display(l::AssuranceLevel) = _LEVEL_DISPLAY[l]

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
    plane::String
    file::String
    testset::String
    mechanism::String
    independence::String
    discriminant::String
    refs::String
    srctext::String
    function Card(
        hub,
        regime,
        arity,
        plane,
        file,
        testset,
        mechanism,
        independence,
        discriminant,
        refs,
        srctext,
    )
        independence in ("structural", "asserted") ||
            throw(ArgumentError("Card: invalid independence $(independence)"))
        arity in ("point", "curve") || throw(ArgumentError("Card: invalid arity $(arity)"))
        plane == "why" || throw(ArgumentError("Card: invalid plane $(plane)"))
        return new(
            hub,
            regime,
            arity,
            plane,
            file,
            testset,
            mechanism,
            independence,
            discriminant,
            refs,
            srctext,
        )
    end
end

include(joinpath(@__DIR__, "_atlas_common.jl"))   # _refs_text (shared)

function _handle_verify!(out, ex, file, testset)
    pos = filter(a -> !(a isa Expr && a.head in (:parameters, :kw)), ex.args)[2:end]
    length(pos) >= 3 || return nothing
    model, qty, bc = pos[1], pos[2], pos[3]
    kw = _kwargs(ex)
    route = get(kw, :route, nothing)
    route isa QuoteNode && (route = route.value)
    route isa Symbol || (route = :unknown)
    refs = _refs_text(get(kw, :refs, nothing))
    params = _literal_params(model)
    reg = SWEEP
    for r in get(MODEL_VOCAB, string(_headsym(model)), Regime[])
        try
            r.predicate(params) && (reg=r; break)
        catch err
            @warn "regime predicate threw — card stays @sweep" model = string(
                _headsym(model)
            ) regime = r.token exception = err
        end
    end
    ind, disc = _independence(route, get(kw, :refs, nothing))
    hub = string(_headsym(model), "/", _headsym(qty), "/", _headsym(bc))
    srctext = try
        replace(string(ex), r"\s+" => " ")
    catch
        "verify(...)"
    end
    push!(
        out,
        Card(
            hub,
            reg.token,
            reg.arity,
            "why",
            file,
            testset,
            string(route),
            ind,
            disc,
            refs,
            srctext,
        ),
    )
    return nothing
end

function _scan_expr!(out, ex, file, testset)
    ex isa Expr || return nothing
    if ex.head === :macrocall && ex.args[1] === Symbol("@testset")
        nm = testset
        for a in ex.args
            a isa String && (nm = a)
        end
        for a in ex.args
            _scan_expr!(out, a, file, nm)
        end
        return nothing
    end
    if ex.head === :call && _headsym(ex) === :verify
        try
            _handle_verify!(out, ex, file, testset)
        catch err
            push!(PARSE_FAILS, (file, string(testset, " :: ", typeof(err))))
        end
        return nothing
    end
    for a in ex.args
        _scan_expr!(out, a, file, testset)
    end
end

function scan_file(path::AbstractString)
    out = Card[]
    rel = replace(relpath(path, _ROOT), _BS => '/')
    try
        src = read(path, String)
        top = parseall(src; filename=path)
        _scan_expr!(out, top, rel, "")
    catch err
        push!(PARSE_FAILS, (rel, string("scan_file :: ", typeof(err))))
        return Card[]
    end
    return out
end

function scan_dir(dir::AbstractString)
    empty!(PARSE_FAILS)
    out = Card[]
    for (root, _, files) in walkdir(dir)
        # test/util_verify/ holds unit tests OF the verify() harness
        # itself; its verify(...) calls target stub models and MUST
        # NOT enter the inventory (else stubs leak into hub counts).
        occursin(joinpath("test", "util_verify"), root) && continue
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
        println(
            io,
            string(
                "{",
                _kv("hub", _jstr(c.hub)),
                ",",
                _kv("regime", _jstr(c.regime)),
                ",",
                _kv("arity", _jstr(c.arity)),
                ",",
                _kv("plane", _jstr(c.plane)),
                ",",
                _kv("mechanism", _jstr(c.mechanism)),
                ",",
                _kv("independence", _jstr(c.independence)),
                ",",
                _kv("discriminant", _jstr(c.discriminant)),
                ",",
                _kv("file", _jstr(c.file)),
                ",",
                _kv("testset", _jstr(c.testset)),
                ",",
                _kv("refs", _jstr(c.refs)),
                ",",
                _kv("srctext", _jstr(c.srctext)),
                "}",
            ),
        )
    end
    return String(take!(io))
end

# Idempotent: only touch the file when content actually changed, so a
# docs build (a pure VIEW) never dirties a clean working tree.
function write_inventory(path, cards)
    s = to_jsonl(cards)
    (!isfile(path) || read(path, String) != s) && write(path, s)
    return nothing
end

end # module
