# test/models/quantum/Hubbard1D/test_hubbard1d_jks_stage_c3_partial.jl
#
# Stage C.3 partial regression: complex analytic continuation of phi (#523).

using Test
using QAtlas
using QAtlas.Hubbard1DJKSNLIE: jks_log_phi, jks_log_phi_complex, jks_phi_complex

@testset "Hubbard1D — JKS Stage C.3 partial (#523)" begin
    @testset "Reduces to real-axis form for s on |s| > 1" begin
        for s in (1.5, 2.0, 5.0, -3.0), beta in (0.1, 1.0)
            real_value = jks_log_phi(s, beta)
            complex_value = jks_log_phi_complex(s + 0im, beta)
            @test isapprox(real(complex_value), real_value; atol=1e-12)
            @test abs(imag(complex_value)) < 1e-12
        end
    end

    @testset "Pure imaginary on the cut |s| < 1" begin
        for s in (0.0, 0.3, 0.5, 0.99), beta in (0.5, 1.0)
            v = jks_log_phi_complex(s + 0im, beta)
            @test abs(real(v)) < 1e-12
        end
    end

    @testset "Complex s off-axis is finite" begin
        for x in (2.0, 3.0, -2.5), alpha in (0.1, 0.3)
            v_real = jks_log_phi_complex(x + 0im, 1.0)
            v_shifted = jks_log_phi_complex(x + alpha * im, 1.0)
            @test isfinite(real(v_shifted))
            @test isfinite(imag(v_shifted))
            @test abs(v_shifted) < 2 * abs(v_real) + 5
        end
    end

    @testset "phi_complex = exp(log_phi_complex)" begin
        for s in (1.5 + 0im, 2.0 + 0.3im, 0.5 + 0im, -3.0 + 0.1im)
            @test jks_phi_complex(s, 1.0) ≈ exp(jks_log_phi_complex(s, 1.0))
        end
    end

    @testset "log_phi_complex validates beta > 0" begin
        @test_throws DomainError jks_log_phi_complex(2.0 + 0im, 0.0)
        @test_throws DomainError jks_log_phi_complex(2.0 + 0im, -1.0)
    end

    @testset "Even symmetry: ln phi(-s) = ln phi(s) on real axis" begin
        for s in (1.5, 3.0, 7.0)
            v_pos = jks_log_phi_complex(s + 0im, 1.0)
            v_neg = jks_log_phi_complex(-s + 0im, 1.0)
            @test isapprox(v_neg, v_pos; atol=1e-12)
        end
    end

    @testset "Conjugation symmetry: ln phi(conj(s)) = conj(ln phi(s))" begin
        for x in (1.5, 3.0), y in (0.2, 0.5)
            s = x + y * im
            v = jks_log_phi_complex(s, 1.0)
            v_conj = jks_log_phi_complex(conj(s), 1.0)
            @test isapprox(v_conj, conj(v); atol=1e-12)
        end
    end
end
