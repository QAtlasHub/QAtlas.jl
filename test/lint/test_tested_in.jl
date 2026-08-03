# test/lint/test_tested_in.jl — every `tested_in` path names a file that exists.
#
# This is the `references` check (test_references_bib.jl) applied to the other
# pointer a registry row carries.  It was missing, and the field had rotted:
# 208 of 450 rows pointed at paths left behind by the test-tree reorganisation
# (`test/models/test_TFIM_thermal.jl` for a file now at
# `test/models/quantum/TFIM/test_TFIM_thermal.jl`), across 56 distinct values.
#
# WHY THIS IS WORTH A TEST rather than a one-off sweep.  `tested_in` is the field
# you follow when a row's number looks wrong and you want to see what pins it.
# Following it landed on nothing — and "no such file" reads exactly like "no
# test", so the natural next move is to write a SECOND test for something already
# covered.  It is also the obvious input to any coverage audit, which would have
# scored those 208 rows as untested while they were in fact fine.
#
# Nothing failed while it rotted, which is the whole point: a stale pointer emits
# no signal at all.

using QAtlas, Test

@testset "every tested_in path resolves" begin
    root = pkgdir(QAtlas)
    rows = [e for e in QAtlas.REGISTRY if e.tested_in !== nothing]
    @test !isempty(rows)

    bad = Tuple{Any,Any,String}[]
    for e in rows, p in e.tested_in
        isfile(joinpath(root, p)) || push!(bad, (e.model, e.quantity, p))
    end
    if !isempty(bad)
        @info "tested_in paths that do not resolve" count = length(bad) first_10 = bad[1:min(
            10, end
        )]
    end
    @test isempty(bad)

    # Shape, not just existence: a bare filename or an absolute path would
    # "resolve" on someone's machine and nowhere else.
    for e in rows, p in e.tested_in
        @test startswith(p, "test/")
        @test !isabspath(p)
        @test endswith(p, ".jl")
    end

    # The field is stored normalised, so no read site has to handle two shapes.
    @test all(e -> e.tested_in isa Vector{String}, rows)
    @test all(e -> !isempty(e.tested_in), rows)
    @test all(e -> allunique(e.tested_in), rows)

    # A `String` is still an accepted spelling at the declaration site — that is
    # what made the field's going plural a non-migration for the other 250 rows.
    @test QAtlas._normalise_tested_in("test/x.jl") == ["test/x.jl"]
    @test QAtlas._normalise_tested_in(["a", "b"]) == ["a", "b"]
    @test QAtlas._normalise_tested_in(nothing) === nothing
    @test QAtlas._normalise_tested_in(String[]) === nothing   # empty is absent, not []

    # Coverage floor. Most rows should name a test; a sharp drop means rows are
    # being added without one, which this file cannot otherwise see.
    frac = length(rows) / length(QAtlas.REGISTRY)
    @test frac > 0.8
end
