# ─────────────────────────────────────────────────────────────────────────────
# test/util_verify/test_verify_json.jl
#
# Regression guard for the verify() evidence-card JSON emitters
# (`test/util/verify.jl` `_json_num` / `_json_arr`). A Complex subject
# (correlators of Hermitian operators return ComplexF64 with ~0 imag)
# was emitted as the raw Julia token `0.408 + 0.0im` — invalid JSON —
# which broke the push:main "Record verification evidence" job (only
# the emit path, never exercised by the fast/no-emit PR CI). These
# tests pin: every emitter output is a valid JSON number/null token,
# never a Julia-repr token. `_json_num`/`_json_arr` are in Main scope
# via the runtests util-include block.
# ─────────────────────────────────────────────────────────────────────────────

using Test

# A strict JSON scalar token: decimal/exponent number, or `null`.
const _JSON_NUM_RX = r"^(-?(0|[1-9]\d*)(\.\d+)?([eE][-+]?\d+)?|-?\d+\.\d+([eE][-+]?\d+)?|null)$"

@testset "verify emit — _json_num is always valid JSON, never Julia repr" begin
    @test _json_num(0.408 + 0.0im) == "0.408"        # the exact main-CI bug
    @test _json_num(1.0) == "1.0"
    @test _json_num(3) == "3.0"
    @test _json_num(1 // 2) == "0.5"
    @test _json_num(-0.0) == "-0.0"
    @test _json_num(2.0 + 1.0im) == "null"            # genuinely complex
    @test _json_num(NaN) == "null"
    @test _json_num(Inf) == "null"
    @test _json_num(-Inf) == "null"

    for z in (0.4 + 0.0im, 1.0 + 0.0im, -2.5 + 0.0im, 1.0e-9 + 0.0im)
        tok = _json_num(z)
        @test !occursin("im", tok)                    # never a complex token
        @test !occursin(" ", tok)
        @test occursin(_JSON_NUM_RX, tok)
    end
    for x in (0.0, -1.25, 1.0e-12, 6.022e23, 12345, 7 // 8)
        @test occursin(_JSON_NUM_RX, _json_num(x))
    end
end

@testset "verify emit — _json_arr routes Complex through _json_num" begin
    @test _json_arr([1.0, 2.0 + 0.0im, 3]) == "[1.0,2.0,3.0]"
    @test _json_arr([NaN, Inf, -Inf]) == "[null,null,null]"
    @test _json_arr(Float64[]) == "[]"
    @test !occursin("im", _json_arr([0.7 + 0.0im, 1.3 + 0.0im]))
    @test _json_arr(["a", "b"]) == "[\"a\",\"b\"]"     # non-numbers still quoted
end
