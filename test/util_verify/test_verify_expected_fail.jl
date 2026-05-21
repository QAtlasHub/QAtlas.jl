# ─────────────────────────────────────────────────────────────────────────────
# test/util_verify/test_verify_expected_fail.jl
#
# Regression guard for the `expected_fail` branch in `test/util/verify.jl`
# for bug-surfacing cards. Three contracts:
#
#   (a) `expected_fail=true` + numerically-failing assertion →
#       `@test_broken` records a Broken result (passes the suite).
#   (b) `expected_fail=true` + numerically-passing assertion →
#       `@test_broken` triggers an "Unexpected Pass" Error (surfaces
#       the bug-fixed case so the card can be promoted back to @test).
#   (c) Under `QATLAS_EMIT=1`, the emitted JSONL card carries
#       `"expected_fail":true` so dashboards can distinguish "tracked
#       broken card" from "unintentional regression".
#
# Without this guard, a future rename of the kwarg or a dispatch break
# would silently no-op the broken status and revert all bug-surfacing
# cards to plain `status:"fail"`.
# ─────────────────────────────────────────────────────────────────────────────

using Test
using QAtlas: QAtlas, AbstractQAtlasModel, fetch, Energy, OBC

# Tiny stub model whose fetch returns a constant. Defined in this file
# only (no registry write); behaviour does not leak across test files.
struct _EFStub <: AbstractQAtlasModel
    value::Float64
end

QAtlas.fetch(m::_EFStub, ::Energy{:per_site}, ::OBC; kwargs...) = m.value

# Capture-only testset: records Pass/Fail/Broken/Error into .results and
# never propagates failures to the parent, so we can assert on the
# inner outcomes without aborting this file's outer testset.
mutable struct _CaptureSet <: Test.AbstractTestSet
    description::String
    results::Vector{Any}
end
_CaptureSet(desc::AbstractString; kw...) = _CaptureSet(String(desc), Any[])
Test.record(ts::_CaptureSet, t) = (push!(ts.results, t); t)
Test.finish(ts::_CaptureSet) = ts

@testset "verify expected_fail — @test_broken branch contracts" begin
    # (a) expected_fail=true + actual mismatch → Broken, not Error.
    let res = @testset _CaptureSet "ef_fail_broken" begin
            verify(
                _EFStub(1.0),
                Energy(:per_site),
                OBC(4);
                route=:ed_finite_size,
                independent=2.0,
                agree_within=1e-9,
                refs=["stub: expected_fail=true with mismatched independent"],
                expected_fail=true,
            )
        end
        @test count(r -> r isa Test.Broken, res.results) == 1
        @test count(r -> r isa Test.Pass, res.results) == 0
        @test count(r -> r isa Test.Error, res.results) == 0
    end

    # (b) expected_fail=true + actual match → "Unexpected Pass" Error.
    let res = @testset _CaptureSet "ef_unexpected_pass" begin
            verify(
                _EFStub(1.0),
                Energy(:per_site),
                OBC(4);
                route=:ed_finite_size,
                independent=1.0,
                agree_within=1e-9,
                refs=["stub: expected_fail=true with matching independent"],
                expected_fail=true,
            )
        end
        @test count(r -> r isa Test.Error, res.results) == 1
        @test count(r -> r isa Test.Broken, res.results) == 0
    end

    # Sanity: expected_fail=false (default) with a match is a plain Pass.
    let res = @testset _CaptureSet "ef_default_pass" begin
            verify(
                _EFStub(1.0),
                Energy(:per_site),
                OBC(4);
                route=:ed_finite_size,
                independent=1.0,
                agree_within=1e-9,
                refs=["stub: default expected_fail=false with matching independent"],
            )
        end
        @test count(r -> r isa Test.Pass, res.results) == 1
        @test count(r -> r isa Test.Broken, res.results) == 0
    end
end

# (c) Emit path — drive two cards under QATLAS_EMIT=1 into a scratch
#     directory and check the produced JSONL carries the new flag.
@testset "verify expected_fail — JSONL emit records the flag" begin
    mktempdir() do dir
        withenv(
            "QATLAS_EMIT" => "1", "QATLAS_CIOUT_DIR" => dir, "QATLAS_TEST_FILES" => ""
        ) do
            local res_broken = @testset _CaptureSet "emit_ef_true" begin
                verify(
                    _EFStub(1.0),
                    Energy(:per_site),
                    OBC(4);
                    route=:ed_finite_size,
                    independent=2.0,
                    agree_within=1e-9,
                    refs=["stub: emit-path expected_fail=true card"],
                    expected_fail=true,
                )
            end
            @test count(r -> r isa Test.Broken, res_broken.results) == 1

            verify(
                _EFStub(2.0),
                Energy(:per_site),
                OBC(4);
                route=:ed_finite_size,
                independent=2.0,
                agree_within=1e-9,
                refs=["stub: emit-path expected_fail=false card"],
            )
        end

        files = filter(f -> endswith(f, ".jsonl"), readdir(dir))
        @test !isempty(files)
        body = join((read(joinpath(dir, f), String) for f in files), "\n")

        @test occursin("\"expected_fail\":true", body)
        @test occursin("\"expected_fail\":false", body)

        @test occursin("\"status\":\"fail\"", body)
        @test occursin("\"status\":\"pass\"", body)
    end
end
