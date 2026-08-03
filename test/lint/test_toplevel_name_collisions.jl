# test/lint/test_toplevel_name_collisions.jl — no top-level name may be defined
# by two different `test_*.jl` files.
#
# Every test file is `include`d into the SAME `Main`, so a top-level helper is a
# global shared with whatever the shard planner places beside it.  #837 defined
# `_G(mf, m0) = …`; `test_hubbard1d_jks_eq53.jl` has `const _G = _U / 4`; Julia
# refuses to define a function over a non-function binding, and CI went red on
# #840 — a branch that touched neither file, because adding an unrelated test
# re-partitioned the shards and put the two together for the first time (#841).
#
# THE POINT IS THAT PASSING CI PROVED NOTHING.  #837 was green.  It proved the
# planner happened to separate the collision that day.  The failure then appeared
# on an innocent change and pointed at a file nobody had opened.  Short names
# (`_Q`, `_G`, `_Z`, `_T`) are exactly the ones several files reach for
# independently, so this is not a rare shape.
#
# Severity is not uniform, and the lint deliberately ignores that:
#
#   * two FUNCTIONS with the same name and arity — silent, last one wins;
#   * two FUNCTIONS with different arities — one generic function with two
#     unrelated meanings, so a future arity change makes one file call the
#     other's method (this was `_ed_thermo`, in the TFIM and AKLT suites);
#   * two `const`s with DIFFERENT values — Julia warns, but a warning does not
#     fail CI;
#   * two `const`s with the same value — harmless today, one edit from not being.
#
# Only the last is benign, and it is benign for a reason that expires.  A single
# rule anyone can obey beats four rules with a severity argument attached.

using QAtlas, Test

# ── the extractor ────────────────────────────────────────────────────────────
#
# `Meta.parseall` NEVER throws: a syntax error comes back as an `:error` /
# `:incomplete` node in the tree.  Without the check below an unparseable file
# would contribute zero names and be reported as clean — the failure mode this
# whole file exists to prevent, reproduced inside the lint itself.
function _lint_toplevel_names(path::AbstractString)
    ex = Meta.parseall(read(path, String); filename=String(path))
    parse_failed = Ref(false)
    check(e) =
        if e isa Expr
            (e.head === :error || e.head === :incomplete) && (parse_failed[] = true)
            foreach(check, e.args)
        end
    check(ex)
    parse_failed[] && return (Symbol[], true)

    names = Symbol[]
    # peel `f(x)`, `f(x)::T`, `f{T}(x)`, `f(x) where T` down to the bare name
    bare(x) =
        if x isa Symbol
            x
        elseif (x isa Expr && x.head in (:call, :(::), :curly, :where))
            bare(x.args[1])
        else
            nothing
        end
    function collect!(e)
        e isa Expr || return nothing
        h = e.head
        if h === :toplevel || h === :block
            foreach(collect!, e.args)          # a bare `begin … end` is still top level
        elseif h === :const || h === :global
            foreach(collect!, e.args)
        elseif h === :function || h === :macro || h === :(=)
            n = bare(e.args[1])
            n === nothing || push!(names, n)
        elseif h === :struct
            n = bare(e.args[2])
            n === nothing || push!(names, n)
        end
        # NOT descended into: `@testset`, function bodies, `let`.  Those introduce
        # a scope, so their definitions do not leak into `Main` and are not the
        # hazard.  Nor `include(...)` — the util/harness files it pulls in are the
        # DECLARED shared namespace, handled below.
        return nothing
    end
    collect!(ex)
    return (unique(names), false)
end

"""
    ALLOWED_SHARED :: Dict{Symbol,String}

Names permitted in more than one `test_*.jl`, each with the reason.

Kept deliberately small and deliberately annotated: an allow-list that grows
without justification is how the rot this lint prevents comes back through the
front door.
"""
const ALLOWED_SHARED = Dict(
    :ABQ =>
        "an identical module alias (`const ABQ = AbstractQAtlas`). A module " *
        "binding has no value to drift into, so the three copies cannot " *
        "disagree the way a constant or a method can.",
)

@testset "no top-level name is defined by two test files" begin
    root = joinpath(pkgdir(QAtlas), "test")
    tests = String[]
    shared = String[]
    for (d, _, fs) in walkdir(root), f in fs
        endswith(f, ".jl") || continue
        p = joinpath(d, f)
        rel = replace(relpath(p, root), '\\' => '/')
        if startswith(rel, "util/") || startswith(rel, "harness/") || startswith(rel, "ci/")
            push!(shared, p)
        elseif startswith(f, "test_")
            push!(tests, p)
        end
    end
    # the sweep really ran — an empty walk would satisfy every assertion below
    @test length(tests) > 200
    @test length(shared) > 5

    owner = Dict{Symbol,Vector{String}}()
    unparseable = String[]
    for p in tests
        (ns, failed) = _lint_toplevel_names(p)
        failed && (push!(unparseable, relpath(p, root)); continue)
        for n in ns
            push!(get!(owner, n, String[]), replace(relpath(p, root), '\\' => '/'))
        end
    end
    @test isempty(unparseable)
    @test length(owner) > 100                      # names were actually extracted

    dups = sort(
        [(k, v) for (k, v) in owner if length(v) > 1 && !haskey(ALLOWED_SHARED, k)];
        by=x -> string(x[1]),
    )
    if !isempty(dups)
        @info "top-level names defined in more than one test file" collisions = dups
    end
    @test isempty(dups)

    # …and a test file must not shadow the util/harness namespace either, which
    # is where `verify`, `_build_tfim_dense` and friends live: shadowing one is
    # the same failure with a longer fuse.
    sharednames = Set{Symbol}()
    for p in shared
        (ns, failed) = _lint_toplevel_names(p)
        failed || union!(sharednames, ns)
    end
    @test !isempty(sharednames)
    shadowed = sort([(k, v) for (k, v) in owner if k in sharednames]; by=x -> string(x[1]))
    if !isempty(shadowed)
        @info "test files shadowing a util/harness name" shadowed = shadowed
    end
    @test isempty(shadowed)

    # every allow-list entry must still be a real collision — a stale exemption
    # silently re-permits a name nobody is using any more
    for (k, why) in ALLOWED_SHARED
        @test haskey(owner, k)
        @test length(owner[k]) > 1
        @test !isempty(why)
    end
end

@testset "the extractor reports a parse failure instead of finding nothing" begin
    # `Meta.parseall` returning an `:error` node rather than throwing is the exact
    # trap this lint would otherwise fall into: no names, so no collisions, so a
    # clean report for a file that does not compile.
    mktemp() do path, io
        write(io, "function broken(\n")
        close(io)
        (ns, failed) = _lint_toplevel_names(path)
        @test failed
        @test isempty(ns)
    end
    mktemp() do path, io
        write(io, "_ok_helper(x) = x + 1\nconst _OK_CONST = 3\nstruct _OkType end\n")
        close(io)
        (ns, failed) = _lint_toplevel_names(path)
        @test !failed
        @test Set(ns) == Set([:_ok_helper, :_OK_CONST, :_OkType])
    end
    # a definition inside a @testset is scoped, so it must NOT be collected
    mktemp() do path, io
        write(io, "using Test\n@testset \"x\" begin\n  _inner(x) = x\nend\n")
        close(io)
        (ns, failed) = _lint_toplevel_names(path)
        @test !failed
        @test !(:_inner in ns)
    end
end
