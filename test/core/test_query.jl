using QAtlas, Test

# Availability search (core/query.jl): "does the atlas have X?" → yes/no + JSONL, by
# model / quantity / bc / regime. These tests pin the contract the LLM "use" face depends on.

@testset "query: search + available agree; facets narrow" begin
    r = QAtlas.search(; model=:TFIM, regime=:finite_temperature)
    @test r isa QAtlas.QueryResult
    @test r.available
    @test !isempty(r.hits)
    @test QAtlas.available(; model=:TFIM, regime=:finite_temperature) == r.available
    @test all(h -> occursin("TFIM", h.model), r.hits)        # every hit really is TFIM
end

@testset "query: regime taxonomy filters correctly" begin
    uni = QAtlas.search(; regime=:universality)
    @test uni.available
    @test all(h -> h.status === :universal, uni.hits)        # :universality ⇒ only :universal rows
    @test QAtlas.search(; regime=:finite_temperature).available
    @test_throws ArgumentError QAtlas.search(; regime=:not_a_regime)
end

@testset "query: rich node facets — status / reliability / method / reference / cross_checked" begin
    ex = QAtlas.search(; model=:TFIM, status=:exact)
    @test ex.available
    @test all(h -> h.status === :exact, ex.hits)            # exact Symbol match
    @test all(h -> h.method isa Symbol, ex.hits)            # method now surfaces in the hit
    pf = QAtlas.search(; reference="Pfeuty")                 # bibkey substring (fuzzy)
    @test pf.available
    @test all(h -> any(r -> occursin("pfeuty", lowercase(r)), h.references), pf.hits)
    @test QAtlas.available(; method=:bethe)                 # fuzzy → finds :bethe_ansatz
    xc = QAtlas.search(; model=:TFIM, cross_checked=true)   # the trust filter
    @test xc.available && all(h -> h.cross_checked, xc.hits)
    @test !QAtlas.available(; model=:TFIM, status=:no_such_status)  # AND-combined, honest miss
    io = IOBuffer()
    QAtlas.search_jsonl(io; model=:TFIM, status=:exact)
    @test occursin("\"method\":", String(take!(io)))        # method in the JSONL hit
end

@testset "query: fuzzy facet matching (case/underscore-insensitive)" begin
    @test QAtlas.available(; model=:tfim)                    # lowercase Symbol
    @test QAtlas.available(; quantity=:specific_heat)        # underscore vs SpecificHeat
    @test QAtlas.available(; quantity="FreeEnergy")          # String facet
    @test !QAtlas.available(; model=:NoSuchModelXYZ)         # honest "not available"
end

@testset "query: a Type facet matches a whole family; a concrete type narrows" begin
    fam = QAtlas.search(; quantity=QAtlas.AbstractThermalPotential)
    fe = QAtlas.search(; quantity=QAtlas.FreeEnergy)
    @test fam.available && fe.available
    @test all(h -> h.quantity == "FreeEnergy", fe.hits)
    @test length(fe.hits) <= length(fam.hits)
end

@testset "query: search_jsonl shape (summary line + one object per hit)" begin
    r = QAtlas.search(; model=:TFIM, regime=:finite_temperature)
    io = IOBuffer()
    QAtlas.search_jsonl(io; model=:TFIM, regime=:finite_temperature)
    lines = split(strip(String(take!(io))), '\n')
    @test length(lines) == length(r.hits) + 1                # summary + one per hit
    @test startswith(lines[1], "{\"available\":true")        # branchable first line
    @test occursin("\"count\":$(length(r.hits))", lines[1])
    @test all(l -> startswith(l, "{\"model\":"), lines[2:end])
    @test all(l -> occursin("\"cross_checked\":", l), lines[2:end])
end

@testset "query: relations — a model's graph neighborhood (edge search)" begin
    rs = QAtlas.relations(:TFIM)
    @test !isempty(rs)
    @test :dual in (r.kind for r in rs)                      # TFIM is Kramers–Wannier self-dual
    @test all(r -> r.from == "TFIM" || r.to == "TFIM", rs)   # every relation touches TFIM
    @test all(r -> r.references isa Vector, rs)
    io = IOBuffer()
    QAtlas.relations_jsonl(io, :TFIM)
    lines = split(strip(String(take!(io))), '\n')
    @test length(lines) == length(rs) + 1                    # summary + one per relation
    @test startswith(lines[1], "{\"available\":true")
    @test all(l -> occursin("\"kind\":", l), lines[2:end])
    io2 = IOBuffer()
    QAtlas.relations_jsonl(io2, :NoSuchModelXYZ)             # no such model → false summary, no lines
    out2 = strip(String(take!(io2)))
    @test startswith(out2, "{\"available\":false")
    @test !occursin('\n', out2)
end

@testset "query: not-available search → single false summary, no hit lines" begin
    io = IOBuffer()
    QAtlas.search_jsonl(io; model=:NoSuchModelXYZ)
    lines = split(strip(String(take!(io))), '\n')
    @test length(lines) == 1
    @test startswith(lines[1], "{\"available\":false")
end
