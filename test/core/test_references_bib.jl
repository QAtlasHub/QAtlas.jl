# test/core/test_references_bib.jl
#
# Key-consistency check: every reference cited in ANY declarative edge store
# (@register / @realizes / @reduces / @about and the constraint stores, via
# the EDGE_STORES registration of core/constraints.jl) must resolve to an
# entry in references.bib (#696).
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
    return get(
        ENV, "QATLAS_REFERENCES_BIB", joinpath(pkgdir(QAtlas), "docs", "references.bib")
    )
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

# Allowlist of references intentionally left as free strings. The bibkey
# migration is now COMPLETE — every `@register` reference resolves to a
# references.bib entry (DOI-bearing where one exists; otherwise @book /
# @article+ISSN / @misc+URL), so this list is empty. A new free-string
# reference therefore fails the completeness test below: add a bibkey to
# references.bib, or (only if genuinely unciteable) list it here with a
# justification.
const KNOWN_UNMIGRATED = Set{String}()

# Every reference cited anywhere in the knowledge graph, via the
# EDGE_STORES registration (one loop covers @register, @realizes, @reduces,
# @about, @symmetry, @identity, @dual, @limits_to — and any future store the
# moment it registers itself).
function all_cited_references()
    refs = String[]
    for spec in QAtlas.EDGE_STORES, row in spec.store
        append!(refs, String(r) for r in spec.references_of(row))
    end
    return refs
end

@testset "references.bib key consistency" begin
    path = references_bib_path()
    @test isfile(path)

    available = bib_citation_keys(path)
    @test !isempty(available)

    # The registration must actually cover the graph (a regression here means
    # a store stopped registering itself).
    @test Set(s.name for s in QAtlas.EDGE_STORES) ⊇
        Set([:registry, :realizes, :reduces, :about])

    # Every bibkey-form reference cited across every edge store must exist.
    cited = Set{String}()
    for r in all_cited_references()
        is_bibkey(r) && push!(cited, r)
    end

    missing_keys = sort(collect(setdiff(cited, available)))
    if !isempty(missing_keys)
        @info "bibkeys cited in edge stores but absent from references.bib" missing_keys
    end
    @test missing_keys == String[]
end

@testset "references.bib completeness (no stray free strings)" begin
    available = bib_citation_keys(references_bib_path())

    # Every edge-store reference must be EITHER a bibkey present in
    # references.bib OR an explicitly-allowlisted unciteable free string.
    stray = Set{String}()
    for s in all_cited_references()
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
    all_refs = Set(all_cited_references())
    dead_allow = sort(collect(setdiff(KNOWN_UNMIGRATED, all_refs)))
    if !isempty(dead_allow)
        @info "KNOWN_UNMIGRATED entries no longer cited — remove them" dead_allow
    end
    @test dead_allow == String[]
end
