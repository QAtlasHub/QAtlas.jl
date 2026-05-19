# ─────────────────────────────────────────────────────────────────────────────
# test/util_verify/test_extrapolate.jl
#
# Unit tests for the finite-size extrapolation helper
# (`test/util/extrapolate.jl`). `test/util_verify/` is the home for tests
# OF the `test/util/` helpers themselves (kept out of the model test
# files). `extrapolate_inf` is in scope via the runtests util-include
# block; this test depends on no QAtlas physics.
# ─────────────────────────────────────────────────────────────────────────────

using Test

@testset "extrapolate_inf — recovers known thermodynamic limits" begin
    @testset "O(1/N) tail (power=1)" begin
        Ns = [8, 16, 32, 64, 128]
        vals = [2.0 + 3.0 / N for N in Ns]
        r = extrapolate_inf(Ns, vals; power=1)
        @test isapprox(r.value, 2.0; atol=1e-9)
        @test abs(r.value - 2.0) <= 10 * r.uncertainty + 1e-9
    end

    @testset "O(1/N²) tail (power=2)" begin
        Ns = [6, 10, 14, 18, 22]
        vals = [-0.5 + 1.7 / N^2 for N in Ns]
        r = extrapolate_inf(Ns, vals; power=2)
        @test isapprox(r.value, -0.5; atol=1e-9)
    end

    @testset "mixed 1/N + 1/N² tail (power=1 still converges)" begin
        Ns = [8, 12, 16, 24, 32, 48]
        vals = [4.0 - 0.9 / N + 2.5 / N^2 for N in Ns]
        r = extrapolate_inf(Ns, vals; power=1)
        @test isapprox(r.value, 4.0; atol=1e-6)
        @test abs(r.value - 4.0) <= 20 * r.uncertainty + 1e-9
    end

    @testset "constant series → exact, ~zero uncertainty" begin
        Ns = [4, 6, 8, 10]
        vals = fill(float(pi), length(Ns))
        r = extrapolate_inf(Ns, vals; power=1)
        @test isapprox(r.value, pi; atol=1e-12)
        @test r.uncertainty <= 1e-10
    end

    @testset "exponential tail is over-bounded, not under (gapped ED)" begin
        # gapped observables converge ~exp(-N/ξ); a 1/N poly fit must
        # still land within its own reported uncertainty of the limit.
        Ns = [6, 8, 10, 12, 14]
        vals = [1.25 + 0.8 * exp(-N / 1.5) for N in Ns]
        r = extrapolate_inf(Ns, vals; power=1)
        @test abs(r.value - 1.25) <= 50 * r.uncertainty + 1e-6
    end

    @testset "two points → linear extrapolation" begin
        r = extrapolate_inf([10, 20], [3.0 + 5.0 / 10, 3.0 + 5.0 / 20]; power=1)
        @test isapprox(r.value, 3.0; atol=1e-9)
    end

    @testset "argument validation" begin
        @test_throws ErrorException extrapolate_inf([8, 16], [1.0]; power=1)
        @test_throws ErrorException extrapolate_inf([8], [1.0]; power=1)
        @test_throws ErrorException extrapolate_inf([0, 8], [1.0, 2.0]; power=1)
        @test_throws ErrorException extrapolate_inf([8, 16], [1.0, 2.0]; power=0)
    end
end
