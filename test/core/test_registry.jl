using Test
using QAtlas
using QAtlas:
    TFIM,
    Energy,
    MassGap,
    CentralCharge,
    FreeEnergy,
    ThermalEntropy,
    SpecificHeat,
    MagnetizationX,
    MagnetizationXLocal,
    MagnetizationZLocal,
    EnergyLocal,
    SusceptibilityXX,
    SusceptibilityZZ,
    ZZStructureFactor,
    ZZCorrelation,
    XXCorrelation,
    VonNeumannEntropy,
    FidelitySusceptibility,
    OBC,
    PBC,
    Infinite,
    Implementation,
    implementation_status,
    implementation_status_markdown,
    has_native_fetch,
    REGISTRY

# Tests for the declarative implementation registry (`src/core/registry.jl`
# + `src/models/quantum/TFIM/TFIM_registry.jl`).  The registry exists so
# that downstream consumers can query "which (model, quantity, bc)
# triples does QAtlas implement, and how reliable are they?" without
# grepping `src/`.  Two safety properties live here:
#
#   1. Query API correctness — filters return the rows the caller asked
#      for and nothing else.
#   2. Drift detection — every registered row corresponds to a real
#      `fetch` method (more specific than the catch-all in
#      `src/core/type.jl`).  This catches the silent regression where
#      a future PR removes a fetch method but forgets the registry row.

@testset "REGISTRY is populated for TFIM at load time" begin
    @test !isempty(REGISTRY)
    @test all(e isa Implementation for e in REGISTRY)
    @test any(
        e.model === TFIM && e.quantity === Energy{:total} && e.bc === OBC for e in REGISTRY
    )
end

@testset "implementation_status() returns Tables.jl-shaped rows" begin
    rows = implementation_status()
    @test rows isa Vector
    @test !isempty(rows)
    sample = first(rows)
    # NamedTuple with the documented field set (status sits between the
    # algorithm tag `method` and the confidence tag `reliability`).
    expected_keys = (
        :model,
        :quantity,
        :bc,
        :method,
        :status,
        :reliability,
        :tested_in,
        :references,
        :notes,
    )
    @test keys(sample) === expected_keys
    @test sample.model isa Type
    @test sample.quantity isa Type
    @test sample.bc isa Type
    @test sample.method isa Symbol
    @test sample.status isa Symbol
    @test sample.reliability isa Symbol
    @test sample.tested_in isa Union{String,Nothing}
    @test sample.references isa Vector{String}
    @test sample.notes isa String
end

@testset "implementation_status filtering" begin
    # By model type
    tfim_rows = implementation_status(TFIM)
    @test !isempty(tfim_rows)
    @test all(r.model === TFIM for r in tfim_rows)

    # By model instance
    @test implementation_status(TFIM(; J=1.0, h=1.0)) == tfim_rows

    # By quantity type
    energy_total_rows = implementation_status(Energy{:total})
    @test !isempty(energy_total_rows)
    @test all(r.quantity === Energy{:total} for r in energy_total_rows)

    # By quantity instance
    @test implementation_status(Energy(:total)) == energy_total_rows

    # Unregistered quantity returns empty (WignerSurmise has no
    # fetch implementation yet, so the registry must not list it).
    @test isempty(implementation_status(WignerSurmise()))
end

@testset "implementation_status(queue) returns one row per registered triple" begin
    m = TFIM(; J=1.0, h=1.0)
    queue = [
        (m, Energy(:total), OBC(8)),       # registered
        (m, MassGap(), Infinite()),   # registered
        (m, WignerSurmise(), OBC(8)),   # NOT registered → dropped
    ]
    rows = implementation_status(queue)
    @test length(rows) == 2
    @test rows[1].quantity === Energy{:total}
    @test rows[2].quantity === MassGap

    # Type-only triples (no instances needed)
    type_queue = [(TFIM, Energy{:per_site}, Infinite), (TFIM, CentralCharge, Infinite)]
    rows_t = implementation_status(type_queue)
    @test length(rows_t) == 2
    @test rows_t[1].quantity === Energy{:per_site}
    @test rows_t[2].quantity === CentralCharge

    # Malformed queue element triggers a clear error
    @test_throws ErrorException implementation_status([(m, Energy(:total))])
end

@testset "Drift guard: every registered TFIM row has a non-catch-all fetch" begin
    # If this fails, a registry row is lying about its backing fetch
    # method — either delete the row or restore the implementation.
    for e in REGISTRY
        @test has_native_fetch(e)
    end
end

@testset "Markdown rendering produces a non-empty GFM table" begin
    io = IOBuffer()
    # Pick TFIM-only rows so this test stays stable as more models get
    # registered ahead of TFIM in include order.
    implementation_status_markdown(io, implementation_status(TFIM)[1:3])
    md = String(take!(io))
    # Header + alignment row + 3 data rows = 5 lines minimum
    lines = split(strip(md), '\n')
    @test length(lines) ≥ 5
    @test startswith(lines[1], "| Model | Quantity | BC")
    # The status column must surface in the rendered header.
    @test occursin("Status", lines[1])
    @test occursin("|---|", lines[2])
    # First TFIM row should be reachable via short-type rendering
    @test any(occursin("TFIM", line) for line in lines[3:end])
end

@testset "Reliability values are drawn from the documented vocabulary" begin
    allowed = (:high, :medium, :low, :not_implemented, :unknown)
    for e in REGISTRY
        @test e.reliability in allowed
    end
end

@testset "Status values are drawn from STATUS_VALUES" begin
    # The status axis (claim kind) is a controlled vocabulary owned by
    # src/core/registry.jl. register! rejects anything outside it, so this
    # is a belt-and-braces guard that every populated row honours the
    # single source of truth.
    for e in REGISTRY
        @test e.status in QAtlas.STATUS_VALUES
    end
    # Rows that omit `status=` default to :exact (e.g. the Energy rows);
    # the v0.24 worked examples declare :bound / :approx explicitly.
    @test all(
        e.status === :exact for e in REGISTRY if e.model === TFIM && e.quantity <: Energy
    )
    @test any(e.status === :bound for e in REGISTRY if e.model === TFIM)
end

@testset "register! rejects a status outside STATUS_VALUES" begin
    # Fail-fast at registration time, not silently at query time. The
    # rejection happens before the push!, so REGISTRY is left untouched.
    n_before = length(REGISTRY)
    @test_throws ArgumentError QAtlas.register!(TFIM, MassGap, Infinite; status=:nonsense)
    # A :bound row must pin a direction; non-bounds must not carry one.
    @test_throws ArgumentError QAtlas.register!(TFIM, MassGap, Infinite; status=:bound)
    @test_throws ArgumentError QAtlas.register!(
        TFIM, MassGap, Infinite; status=:bound, direction=:sideways
    )
    @test_throws ArgumentError QAtlas.register!(
        TFIM, MassGap, Infinite; status=:exact, direction=:upper
    )
    @test length(REGISTRY) == n_before
end

@testset "Method values are non-:unknown for the populated TFIM rows" begin
    # All TFIM rows in this PR have a real algorithm tag.  This test
    # pins the convention so future TFIM additions can't silently leave
    # `method=:unknown` behind.
    for e in REGISTRY
        e.model === TFIM || continue
        @test e.method !== :unknown
    end
end

@testset "references_for returns the registry bibkeys for a triple" begin
    # TFIM per-site energy at infinite size rests on Pfeuty 1970 (BdG).
    @test references_for(TFIM(), Energy{:per_site}(), Infinite()) == ["Pfeuty1970"]

    # Types and instances are interchangeable, mirroring implementation_status.
    @test references_for(TFIM, Energy{:per_site}, Infinite) ==
        references_for(TFIM(), Energy{:per_site}(), Infinite())

    # Unregistered / reference-free triples yield an empty vector, never throw.
    @test references_for(TFIM(), MeanRatio(), OBC(8)) == String[]
end

@testset "references_for aggregates over unspecified axes" begin
    # (model, quantity): union across boundary conditions, deduped + sorted.
    pair = references_for(TFIM(), Energy{:per_site}())
    @test "Pfeuty1970" in pair
    @test issorted(pair) && allunique(pair)

    # (model): every paper any TFIM row cites — a superset of any narrower query.
    all_tfim = references_for(TFIM())
    @test issorted(all_tfim) && allunique(all_tfim)
    @test "Pfeuty1970" in all_tfim
    @test issubset(Set(pair), Set(all_tfim))

    # Cross-check against the source of truth: identical to folding the
    # registry rows for TFIM by hand.
    expected = sort!(
        unique!(String[r for e in REGISTRY if e.model === TFIM for r in e.references])
    )
    @test all_tfim == expected
end
