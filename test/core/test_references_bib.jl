# test/core/test_references_bib.jl
#
# Key-consistency check: every reference cited in the `@register` REGISTRY
# must resolve to an entry in references.bib.
#
# This is the first half of the two-stage reference check:
#   1. (here, Julia) every cited bibkey resolves to an entry in
#      references.bib — catches typos and dangling keys.
#   2. (CI, doiget verify) every references.bib entry's DOI / arXiv id
#      actually exists upstream — see .github/workflows/VerifyReferences.yml.
#
# The registry has been migrated from free author-year strings
# (e.g. "Onsager 1944") to bibkeys (e.g. "AKLT1988"). The migration is
# complete except for a small, explicit allowlist of references that have
# no citable DOI / arXiv record (textbooks, pre-DOI Soviet journals, or
# author strings too ambiguous to pin). Any OTHER free string is a test
# failure: a new reference must either be added to references.bib with a
# bibkey, or — if genuinely unciteable — added to `KNOWN_UNMIGRATED` below
# with a justification.

using QAtlas, Test

"""
    references_bib_path() -> String

Path to the canonical references.bib. Overridable via the
`QATLAS_REFERENCES_BIB` environment variable so CI (and `doiget verify`)
can point at the same file from any working directory.
"""
function references_bib_path()
    get(ENV, "QATLAS_REFERENCES_BIB", joinpath(pkgdir(QAtlas), "docs", "references.bib"))
end

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

A reference is in bibkey form when it is a single run of `[A-Za-z0-9_]`
starting with a letter — i.e. no spaces/punctuation.
"""
is_bibkey(s::AbstractString) = occursin(r"^[A-Za-z][A-Za-z0-9_]*$", s)

# References intentionally left as free strings: no reliable Crossref DOI /
# arXiv id exists to make them machine-verifiable. Each is a textbook, a
# pre-DOI Soviet-journal paper, or an author string too ambiguous to pin.
# Keep this list tight — prefer adding a bibkey to references.bib.
const KNOWN_UNMIGRATED = Set{String}([
    "Ashcroft-Mermin 1976",          # textbook, Solid State Physics (no DOI)
    "Baxter 1982",                   # textbook, Exactly Solved Models (no DOI)
    "Chen-Roncaglia 2008",           # ambiguous author string, no unique match
    "Fateev-Zamolodchikov 1985",     # Sov. Phys. JETP 62, 215 — pre-DOI
    "Kitaev 2015",                   # SYK KITP talks (no DOI)
    "Landau 1937",                   # ZhETF 7, 19 — pre-DOI original
    "Landau-Lifshitz §149",          # textbook, Statistical Physics (no DOI)
    "Stanley 1971",                  # textbook, Intro to Phase Transitions (no DOI)
    "Tzeng-Yang-Hsu 2017",           # ambiguous author string, no unique match
])

@testset "references.bib key consistency" begin
    path = references_bib_path()
    @test isfile(path)

    available = bib_citation_keys(path)
    @test !isempty(available)

    # Every bibkey-form reference cited across the registry must exist.
    cited = Set{String}()
    for e in QAtlas.REGISTRY, r in e.references
        is_bibkey(r) && push!(cited, String(r))
    end

    missing_keys = sort(collect(setdiff(cited, available)))
    if !isempty(missing_keys)
        @info "bibkeys cited in REGISTRY but absent from references.bib" missing_keys
    end
    @test missing_keys == String[]
end

@testset "references.bib completeness (no stray free strings)" begin
    available = bib_citation_keys(references_bib_path())

    # Every registry reference must be EITHER a bibkey present in
    # references.bib OR an explicitly-allowlisted unciteable free string.
    stray = Set{String}()
    for e in QAtlas.REGISTRY, r in e.references
        s = String(r)
        if is_bibkey(s)
            s in available || push!(stray, s)         # dangling bibkey
        else
            s in KNOWN_UNMIGRATED || push!(stray, s)   # un-allowlisted free string
        end
    end

    strays = sort(collect(stray))
    if !isempty(strays)
        @info "references not resolvable to references.bib (add a bibkey, or " *
            "allowlist in KNOWN_UNMIGRATED if genuinely unciteable)" strays
    end
    @test strays == String[]

    # Guard against the allowlist rotting: every KNOWN_UNMIGRATED entry must
    # still be cited somewhere (otherwise drop it).
    all_refs = Set(String(r) for e in QAtlas.REGISTRY for r in e.references)
    dead_allow = sort(collect(setdiff(KNOWN_UNMIGRATED, all_refs)))
    if !isempty(dead_allow)
        @info "KNOWN_UNMIGRATED entries no longer cited — remove them" dead_allow
    end
    @test dead_allow == String[]
end
