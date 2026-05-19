# test/harness/atlas/test_atlas_logic.jl
# Unit tests for the R1 assurance taxonomy (the policy logic that decides
# which hubs are flagged as actionable risk) and the static registry
# parser. Runs in the dedicated `inventory-drift` CI job — test/harness/
# is intentionally excluded from the sharded universe guard.
using Test
include(joinpath(@__DIR__, "AtlasInventory.jl"))
include(joinpath(@__DIR__, "AtlasRegistry.jl"))
using .AtlasInventory:
    assurance_level,
    level_display,
    AssuranceLevel,
    UNIVERSALITY_CORROBORATED,
    CORROBORATED_AT_P,
    COHERENT,
    CITED_ONLY,
    UNCORROBORATED_BUT_FEASIBLE,
    ED_INFEASIBLE_MODELS
using .AtlasRegistry: scan_registry

@testset "assurance_level — mechanism → level precedence" begin
    @test assurance_level(Set(["universality_consistency"]), false) ==
        UNIVERSALITY_CORROBORATED
    @test assurance_level(Set(["ed_finite_size"]), false) == CORROBORATED_AT_P
    @test assurance_level(Set(["second_closed_form"]), false) == CORROBORATED_AT_P
    @test assurance_level(Set(["sum_rule"]), false) == COHERENT
    @test assurance_level(Set(["limiting_case"]), false) == COHERENT
    @test assurance_level(Set(["delegation_invariant"]), false) == COHERENT
    @test assurance_level(Set(["retype_formula"]), false) == COHERENT
    @test assurance_level(Set(["unknown"]), false) == COHERENT
    @test assurance_level(Set(["literature_value"]), false) == CITED_ONLY
    # Highest tier wins when several mechanisms are present.
    @test assurance_level(Set(["ed_finite_size", "sum_rule"]), false) == CORROBORATED_AT_P
    @test assurance_level(Set(["universality_consistency", "ed_finite_size"]), false) ==
        UNIVERSALITY_CORROBORATED
    @test assurance_level(Set(["literature_value", "sum_rule"]), false) == COHERENT
    # Empty: feasible → actionable risk; infeasible → cited-only frontier.
    @test assurance_level(Set{String}(), false) == UNCORROBORATED_BUT_FEASIBLE
    @test assurance_level(Set{String}(), true) == CITED_ONLY
    # ed-infeasible never *upgrades* a hub that already has cards.
    @test assurance_level(Set(["ed_finite_size"]), true) == CORROBORATED_AT_P
    @test assurance_level(Set(["literature_value"]), true) == CITED_ONLY
end

@testset "level_display — exhaustive, stable strings" begin
    @test length(instances(AssuranceLevel)) == 5
    for l in instances(AssuranceLevel)
        name, badge, adm = level_display(l)
        @test !isempty(name)
        @test !isempty(badge)
        @test adm in ("tip", "note", "warning")
    end
    @test level_display(UNIVERSALITY_CORROBORATED) ==
        ("universality-corroborated", "🟣", "tip")
    @test level_display(CORROBORATED_AT_P) == ("corroborated-at-p", "🟢", "tip")
    @test level_display(COHERENT) == ("coherent", "🔵", "note")
    @test level_display(CITED_ONLY) == ("cited-only", "⚪", "note")
    @test level_display(UNCORROBORATED_BUT_FEASIBLE) ==
        ("uncorroborated-but-feasible", "🟠", "warning")
end

@testset "ED_INFEASIBLE_MODELS — expected frontier set" begin
    @test "KitaevHoneycomb" in ED_INFEASIBLE_MODELS
    @test "KagomeHeisenbergAFM" in ED_INFEASIBLE_MODELS
    @test !("TFIM" in ED_INFEASIBLE_MODELS)
    @test length(ED_INFEASIBLE_MODELS) == 9
end

@testset "scan_registry — synthetic @register (comma-style, real syntax)" begin
    mktempdir() do d
        f = joinpath(d, "Foo_registry.jl")
        write(
            f,
            """
            @register(
                Foo,
                Energy,
                Infinite,
                method=:bethe_ansatz,
                reliability=:high,
                references=["Foo 1999", "Bar 2001"],
                notes="synthetic",
            )
            @register(Foo, MassGap, OBC, method=:bdg, reliability=:medium)
            """,
        )
        claims = scan_registry(f)
        @test length(claims) == 2
        c1 = first(claims)
        @test c1.hub == "Foo/Energy/Infinite"
        @test c1.model == "Foo"
        @test c1.quantity == "Energy"
        @test c1.bc == "Infinite"
        @test c1.method == "bethe_ansatz"
        @test c1.reliability == "high"
        @test c1.refs == "Foo 1999 | Bar 2001"
        @test c1.notes == "synthetic"
        @test claims[2].hub == "Foo/MassGap/OBC"
        @test claims[2].refs == ""
    end
end

@testset "scan_dir resets PARSE_FAILS (in-process idempotency)" begin
    push!(AtlasInventory.PARSE_FAILS, ("dummy", "stale"))
    AtlasInventory.scan_dir(mktempdir())
    @test isempty(AtlasInventory.PARSE_FAILS)
end
