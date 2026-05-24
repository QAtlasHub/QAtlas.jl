# test/harness/atlas/test_doc_structure.jl
#
# Structural consistency guard for the Tier-1 docs view (ModelList +
# per-model + per-quantity matrices).  Asserts that the regenerated atlas
# is 1:1 with the substrate:
#
#   1. Every model named in the registry has a docs/src/atlas/models/<Model>.md
#   2. Every quantity has a docs/src/atlas/quantities/<slug>.md
#   3. ModelList.md exists and contains one row per model
#   4. No orphan files (residue from a removed model / quantity)
#   5. Every per-hub page back-links to its model index
#
# Run after `julia --project=docs docs/atlas/generate.jl`.

using Test

const _ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(_ROOT, "test", "harness", "atlas", "AtlasInventory.jl"))
include(joinpath(_ROOT, "test", "harness", "atlas", "AtlasRegistry.jl"))
using .AtlasInventory, .AtlasRegistry

@testset "Tier-1 docs structure (ModelList + matrices)" begin
    claims = AtlasRegistry.Claim[]
    for (root, _, fs) in walkdir(joinpath(_ROOT, "src"))
        for f in fs
            endswith(f, "_registry.jl") || continue
            append!(claims, AtlasRegistry.scan_registry(joinpath(root, f)))
        end
    end
    @test !isempty(claims)

    claimed = sort(unique(c.hub for c in claims))
    modelof(h) = first(split(h, "/"))
    quantof(h) = (p=split(h, "/"); length(p) >= 2 ? p[2] : "?")
    _quant_slugof(q::AbstractString) = replace(q, r"[^A-Za-z0-9]" => "_")

    expected_models = sort(unique(modelof(h) for h in claimed))
    expected_quants = sort(unique(quantof(h) for h in claimed))

    models_dir = joinpath(_ROOT, "docs", "src", "atlas", "models")
    quants_dir = joinpath(_ROOT, "docs", "src", "atlas", "quantities")
    list_path = joinpath(_ROOT, "docs", "src", "atlas", "ModelList.md")

    @testset "models/<Model>.md exists for every claimed model" begin
        @test isdir(models_dir)
        for m in expected_models
            p = joinpath(models_dir, m * ".md")
            @test isfile(p)
        end
    end

    @testset "no orphan models/*.md" begin
        actual = sort(
            replace.(filter(f -> endswith(f, ".md"), readdir(models_dir)), ".md" => "")
        )
        @test actual == expected_models
    end

    @testset "quantities/<slug>.md exists for every claimed quantity" begin
        @test isdir(quants_dir)
        for q in expected_quants
            p = joinpath(quants_dir, _quant_slugof(q) * ".md")
            @test isfile(p)
        end
    end

    @testset "no orphan quantities/*.md" begin
        actual = sort(
            replace.(filter(f -> endswith(f, ".md"), readdir(quants_dir)), ".md" => "")
        )
        expected_slugs = sort(unique(_quant_slugof.(expected_quants)))
        @test actual == expected_slugs
    end

    @testset "ModelList.md exists with one row per model" begin
        @test isfile(list_path)
        body = read(list_path, String)
        for m in expected_models
            marker = string("[`", m, "`](models/", m, ".md)")
            @test occursin(marker, body)
        end
    end

    @testset "per-hub pages back-link to their model index" begin
        hubsdir = joinpath(_ROOT, "docs", "src", "atlas", "hubs")
        @test isdir(hubsdir)
        sample_misses = String[]
        for h in claimed
            slug = replace(h, r"[^A-Za-z0-9]" => "_")
            p = joinpath(hubsdir, slug * ".md")
            isfile(p) || (push!(sample_misses, "missing hub page: " * h); continue)
            body = read(p, String)
            marker = string("../models/", modelof(h), ".md")
            occursin(marker, body) || push!(sample_misses, "no model back-link: " * h)
        end
        if !isempty(sample_misses)
            @info "doc-structure misses (first 5)" miss = first(sample_misses, 5)
        end
        @test isempty(sample_misses)
    end
end
