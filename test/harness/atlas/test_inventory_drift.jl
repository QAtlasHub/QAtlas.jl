# test/harness/atlas/test_inventory_drift.jl
# M1 drift guard (repo-wide): committed test/INVENTORY.jsonl MUST equal a
# fresh whole-test static scan. Shard-independent, no test execution,
# fetch/register framework untouched.
using Test
include(joinpath(@__DIR__, "AtlasInventory.jl"))
using .AtlasInventory

@testset "INVENTORY.jsonl drift guard (repo-wide)" begin
    root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    cards = AtlasInventory.scan_dir(joinpath(root, "test"))
    fresh = AtlasInventory.to_jsonl(cards)
    committed = read(joinpath(root, "test", "INVENTORY.jsonl"), String)
    if fresh != committed
        @error "INVENTORY drift: regenerate via docs/atlas/generate.jl"
    end
    @test fresh == committed
    @test isempty(AtlasInventory.PARSE_FAILS)
end
