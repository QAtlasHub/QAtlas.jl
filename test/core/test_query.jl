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

@testset "query: orthogonal axes — thermal × dynamical (captured at @register)" begin
    ft = QAtlas.search(; thermal=:finite)
    @test ft.available
    @test all(h -> h.thermal in (:finite, :both), ft.hits)  # :both spans, matches :finite query
    @test QAtlas.available(; thermal=:zero)
    tr = QAtlas.search(; dynamical=:transport)               # velocities are :transport, not :dynamic
    @test tr.available
    @test all(h -> h.dynamical === :transport, tr.hits)
    # THE FIX: real dynamics (:dynamic) is honestly absent → no more velocity overclaim
    @test !QAtlas.available(; dynamical=:dynamic)
    @test !QAtlas.available(; thermal=:finite, dynamical=:dynamic)  # the conjunction, honestly empty
    en = QAtlas.search(; quantity=QAtlas.Energy)             # Energy is :both (T=0 or thermal)
    @test en.available && any(h -> h.thermal === :both, en.hits)
    io = IOBuffer()
    QAtlas.search_jsonl(io; thermal=:finite)
    @test occursin("\"thermal\":", String(take!(io)))        # axes surface in the JSONL hit
end

@testset "query: hit carries the fetch call signature + notes/valid_domain" begin
    r = QAtlas.search(; model=:TFIM)
    @test r.available
    # every hit now carries the three actionability fields, well-typed
    @test all(h -> h.params isa Vector{String}, r.hits)
    @test all(h -> h.notes isa String && h.valid_domain isa String, r.hits)
    @test all(h -> h.params == sort(unique(h.params)), r.hits)   # sorted-unique, no `kwargs` slurp
    @test all(h -> !("kwargs" in h.params), r.hits)
    @test any(h -> !isempty(h.params), r.hits)                   # the call signature is exposed
    @test any(h -> !isempty(h.notes), r.hits)                    # notes surface
    # a finite-temperature quantity exposes its temperature kwarg → "how to call it"
    sh = QAtlas.search(; model=:TFIM, quantity=:SpecificHeat)
    @test sh.available
    @test any(
        h -> any(
            p ->
                occursin("beta", lowercase(p)) || p == "β" || lowercase(p) == "temperature",
            h.params,
        ),
        sh.hits,
    )
    # JSONL surfaces the new fields (additive keys)
    io = IOBuffer()
    QAtlas.search_jsonl(io; model=:TFIM, quantity=:SpecificHeat)
    s = String(take!(io))
    @test occursin("\"params\":", s)
    @test occursin("\"notes\":", s)
    @test occursin("\"valid_domain\":", s)
end

@testset "query: family facet is total; regime no longer drops correlations" begin
    # family is a TOTAL classification — every registered hub gets one
    all_hits = QAtlas.search().hits
    @test all(h -> h.family isa Symbol && h.family !== Symbol(""), all_hits)
    # the family facet finds correlations — the class `regime` used to silently drop
    corr = QAtlas.search(; family=:correlation)
    @test corr.available && all(h -> h.family === :correlation, corr.hits)
    @test all(h -> occursin("Correlation", h.quantity), corr.hits)
    @test QAtlas.available(; model=:Heisenberg, family=:correlation)
    # THE FIX: regime=:ground_state now surfaces T=0-accessible correlations
    gs = QAtlas.search(; model=:Heisenberg1D, regime=:ground_state)
    @test any(h -> occursin("Correlation", h.quantity), gs.hits)
    # NO REGRESSION: the previous ground_state members are still ground_state
    @test QAtlas.available(; quantity=:MassGap, regime=:ground_state)                   # <:AbstractGap
    @test QAtlas.available(; quantity=:GroundStateEnergyDensity, regime=:ground_state)  # thermal=:unknown, kept explicitly
    # family surfaces in the JSONL hit
    io = IOBuffer()
    QAtlas.search_jsonl(io; model=:Heisenberg1D, family=:correlation)
    @test occursin("\"family\":\"correlation\"", String(take!(io)))
end

@testset "query: query_schema — the query convention is self-describing" begin
    sch = QAtlas.query_schema()
    @test sch isa Vector{QAtlas.Facet}
    names = [f.name for f in sch]
    @test all(in(names), (:model, :quantity, :bc, :family, :thermal, :dynamical, :regime))
    # enumerated facets carry their valid values, drawn live from the registry
    fam = only(f for f in sch if f.name === :family)
    @test fam.kind === :structural && "correlation" in fam.values
    th = only(f for f in sch if f.name === :thermal)
    @test th.kind === :axis && Set(th.values) == Set(["zero", "finite", "both", "unknown"])
    bcf = only(f for f in sch if f.name === :bc)
    @test "OBC" in bcf.values && "Infinite" in bcf.values
    # JSONL shape: header + one facet per line
    io = IOBuffer()
    QAtlas.query_schema_jsonl(io)
    lines = split(strip(String(take!(io))), '\n')
    @test startswith(lines[1], "{\"facets\":")
    @test length(lines) == length(sch) + 1
    @test all(l -> occursin("\"name\":", l) && occursin("\"values\":", l), lines[2:end])
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

@testset "query: gaps — per-model coverage holes (absence search)" begin
    g = QAtlas.gaps(:TFIM)
    @test all(x -> x.kind === :regime && x.model == "TFIM", g)
    # grounded: gaps + covered regimes partition all REGIMES (no guessed "expected set")
    covered = count(r -> QAtlas.available(; model=:TFIM, regime=r), keys(QAtlas.REGIMES))
    @test length(g) + covered == length(QAtlas.REGIMES)
    # each reported gap is genuinely absent
    @test all(x -> !QAtlas.available(; model=:TFIM, regime=Symbol(x.subject)), g)
    io = IOBuffer()
    QAtlas.gaps_jsonl(io, :TFIM)
    lines = split(strip(String(take!(io))), '\n')
    @test length(lines) == length(g) + 1
    @test startswith(lines[1], "{\"has_gaps\":")
end

@testset "query: describe — the full per-model grounding record" begin
    rs = QAtlas.describe(:TFIM)
    @test length(rs) == 1
    r = rs[1]
    @test r.model == "TFIM"
    @test !isempty(r.quantities)                 # the observables — disambiguating structural content
    @test !isempty(r.relations)                  # the graph neighborhood
    @test r.summary isa String                   # "" if uncarded, honest
    # the dogfood trap: FibonacciAnyons is uncarded (no prose), but the structure still identifies it
    fib = QAtlas.describe(:fibonacci)
    @test length(fib) == 1 && fib[1].model == "FibonacciAnyons"
    @test fib[1].summary == ""                   # honestly uncarded (no @about card)
    # JSONL: a header line then one rich record per model
    io = IOBuffer()
    QAtlas.describe_jsonl(io, :TFIM)
    lines = split(strip(String(take!(io))), '\n')
    @test length(lines) == 2                      # header + 1 record
    @test startswith(lines[1], "{\"count\":1")
    @test occursin("\"quantities\":", lines[2]) && occursin("\"relations\":", lines[2])
end

@testset "query: realizing(class) — inverse edge query (class → models)" begin
    is = QAtlas.realizing(:Ising)
    @test !isempty(is)
    @test all(r -> r.kind === :realizes && r.to == "Ising", is)
    @test "TFIM" in (r.from for r in is)          # TFIM realizes the Ising class
    @test [r.from for r in QAtlas.realizing(:ising)] == [r.from for r in is]  # case-insensitive
    @test all(r -> r.to == "Ising", QAtlas.realizing(:ising))  # exact: does NOT catch IsingSDRG
    io = IOBuffer()
    QAtlas.realizing_jsonl(io, :Ising)
    lines = split(strip(String(take!(io))), '\n')
    @test length(lines) == length(is) + 1
    @test startswith(lines[1], "{\"available\":true")
    @test all(l -> occursin("\"kind\":\"realizes\"", l), lines[2:end])
end

@testset "query: not-available search → single false summary, no hit lines" begin
    io = IOBuffer()
    QAtlas.search_jsonl(io; model=:NoSuchModelXYZ)
    lines = split(strip(String(take!(io))), '\n')
    @test length(lines) == 1
    @test startswith(lines[1], "{\"available\":false")
end
