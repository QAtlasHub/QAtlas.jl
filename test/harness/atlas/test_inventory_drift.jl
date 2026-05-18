# test/harness/atlas/test_inventory_drift.jl
#
# M1 drift guard: the committed test/INVENTORY.jsonl MUST equal a fresh
# static scan. Shard-independent (no test execution), single job. Run via
# the next-scoped CI job, not the sharded WHY plane.
using Test
include(joinpath(@__DIR__, "AtlasInventory.jl"))
using .AtlasInventory

@testset "INVENTORY.jsonl drift guard (TFIM slice)" begin
    root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    cards = AtlasInventory.scan_dir(joinpath(root, "test", "models", "quantum", "TFIM"))
    fresh = AtlasInventory.to_jsonl(cards)
    committed = read(joinpath(root, "test", "INVENTORY.jsonl"), String)
    if fresh != committed
        @error "INVENTORY.jsonl drift: regenerate with docs/atlas/generate.jl or the scanner"
    end
    @test fresh == committed
end
