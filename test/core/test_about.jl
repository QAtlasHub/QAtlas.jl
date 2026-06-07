# model description cards (@about / about) + the static scanner the docs
# generator uses to render them.

using QAtlas, Test
using QAtlas: TFIM, Heisenberg1D, about, about!, ABOUT, ModelCard

@testset "model description cards (@about / about)" begin
    @testset "about(model) returns the authored card" begin
        c = about(TFIM)
        @test c !== nothing
        @test occursin("transverse-field Ising", c.summary)
        @test occursin(raw"\sigma", c.hamiltonian)      # LaTeX backslashes preserved
        @test occursin(raw"$h = J$", c.summary)          # inline math survives (raw"")
        @test about(TFIM()) == c                          # instance or type
    end

    @testset "about(model) is nothing when no card was authored" begin
        @test about(MajumdarGhosh) !== nothing            # seeded
        @test about(SpinIce) === nothing                  # not seeded -> docstring fallback
    end

    @testset "about! validates and appends" begin
        n = length(ABOUT)
        about!(TFIM; summary="unit-test summary", hamiltonian=raw"H = 0")
        @test length(ABOUT) == n + 1
        @test_throws ArgumentError about!(TFIM; summary="   ")  # empty summary rejected
    end
end

# The docs generator (docs/atlas/generate.jl) does NOT load QAtlas — it reads
# the cards by a static AST scan. Exercise that exact production path.
module _AboutScanProbe
include(joinpath(@__DIR__, "..", "harness", "atlas", "AtlasRegistry.jl"))
end

@testset "AtlasRegistry.scan_about (static path used by the atlas)" begin
    reg = joinpath(@__DIR__, "..", "..", "src", "about_registry.jl")
    cards = _AboutScanProbe.AtlasRegistry.scan_about(reg)
    @test !isempty(cards)
    bym = Dict(c.model => c for c in cards)
    @test haskey(bym, "TFIM")
    @test occursin(raw"\sum", bym["TFIM"].hamiltonian)   # raw"" unwrapped, backslash kept
    @test occursin("Ising", bym["TFIM"].summary)
    @test _AboutScanProbe.AtlasRegistry.scan_about("/no/such/file.jl") == []
end
