# test/core/test_references_bib.jl
#
# Key-consistency check: every bibkey-form reference cited in the
# `@register` REGISTRY must exist in references.bib.
#
# This is the first half of the two-stage reference check:
#   1. (here, Julia) every cited bibkey resolves to an entry in
#      references.bib — catches typos and dangling keys.
#   2. (CI, doiget verify) every references.bib entry's DOI / arXiv id
#      actually exists upstream — see .github/workflows/VerifyReferences.yml.
#
# Migration note: the registry is being migrated from free author-year
# strings (e.g. "Onsager 1944") to bibkeys (e.g. "AKLT1988"). Only
# bibkey-form references are checked here; not-yet-migrated free strings
# (which contain spaces) are skipped, so this test stays green during the
# migration and automatically widens its coverage as more models switch.

using QAtlas, Test

"""
    references_bib_path() -> String

Path to the canonical references.bib. Overridable via the
`QATLAS_REFERENCES_BIB` environment variable so CI (and `doiget verify`)
can point at the same file from any working directory.
"""
references_bib_path() = get(
    ENV,
    "QATLAS_REFERENCES_BIB",
    joinpath(pkgdir(QAtlas), "docs", "references.bib"),
)

"""
    bib_citation_keys(path) -> Set{String}

Parse the set of citation keys from a BibTeX file: the `KEY` in each
`@entrytype{KEY, ...}` line.
"""
function bib_citation_keys(path::AbstractString)
    keys = Set{String}()
    for line in eachline(path)
        m = match(r"^@\w+\{\s*([^,\s]+)\s*,", line)
        m === nothing || push!(keys, String(m.captures[1]))
    end
    return keys
end

"""
    is_bibkey(s) -> Bool

A reference is in bibkey form (vs. a not-yet-migrated free string) when it
is a single run of `[A-Za-z0-9_]` starting with a letter — i.e. no spaces.
"""
is_bibkey(s::AbstractString) = occursin(r"^[A-Za-z][A-Za-z0-9_]*$", s)

@testset "references.bib key consistency" begin
    path = references_bib_path()
    @test isfile(path)

    available = bib_citation_keys(path)
    @test !isempty(available)

    # Collect the bibkey-form references cited across the whole registry.
    cited = Set{String}()
    for e in QAtlas.REGISTRY, r in e.references
        is_bibkey(r) && push!(cited, String(r))
    end

    # Every cited bibkey must exist in references.bib.
    missing_keys = sort(collect(setdiff(cited, available)))
    if !isempty(missing_keys)
        @info "bibkeys cited in REGISTRY but absent from references.bib" missing_keys
    end
    @test missing_keys == String[]
end
