# =============================================================================
# TFIM dynamics — small-N ED correctness + sanity + spreading
#
# Split out of the legacy test_TFIM_dynamics.jl (PR vs next: refactor for
# CI shard balance — original file was ~17 min on s14 because Layers 2b+2c
# ran heavy BdG Pfaffian sweeps at N=200-240 inside the same file).
#
# Helpers _build_tfim_dense, _op_site, _SZ, _SX come from test/util/
# tfim_dense_ed.jl via runtests.jl ambient include.
# =============================================================================

@testset "TFIM dynamics — small-N ED + sanity + spreading" begin

    # ---- Layer 1: ED comparison, small N ------------------------------------
    @testset "ED comparison N=4, J=1, h=0.5" begin
        N, J, h = 4, 1.0, 0.5
        H = _build_tfim_dense(N, J, h)
        E, V = eigen(H)
        gs = V[:, 1]

        # Sanity: GS energy matches QAtlas BdG.
        @test E[1] ≈ QAtlas.fetch(TFIM(; J=J, h=h), Energy(), OBC(; N=N)) atol=1e-10

        # σz σz at t = 0
        for i in 1:N, j in i:N
            ED = real(gs' * _op_site(_SZ, i, N) * _op_site(_SZ, j, N) * gs)
            QAT = real(
                QAtlas.fetch(
                    TFIM(; J=J, h=h),
                    DynamicalCorrelation(:z, :z),
                    OBC(; N=N);
                    i=i,
                    j=j,
                    t=0.0,
                ),
            )
            @test QAT ≈ ED atol=1e-10
        end

        # σx σx at t = 0
        for i in 1:N, j in i:N
            ED = real(gs' * _op_site(_SX, i, N) * _op_site(_SX, j, N) * gs)
            QAT = real(
                QAtlas.fetch(
                    TFIM(; J=J, h=h),
                    DynamicalCorrelation(:x, :x),
                    OBC(; N=N);
                    i=i,
                    j=j,
                    t=0.0,
                ),
            )
            @test QAT ≈ ED atol=1e-10
        end

        # Unequal-time σz σz, several (i, j, t) — full complex check
        Udag = V'  # = V^†
        for (i, j, t) in [(1, 1, 0.3), (2, 3, 0.7), (1, 4, 1.5), (2, 4, 2.0)]
            Ut = V * Diagonal(exp.(-im * E * t)) * Udag
            ED =
                exp(im * E[1] * t) *
                (gs' * _op_site(_SZ, i, N) * Ut * _op_site(_SZ, j, N) * gs)
            QAT = QAtlas.fetch(
                TFIM(; J=J, h=h), DynamicalCorrelation(:z, :z), OBC(; N=N); i=i, j=j, t=t
            )
            @test real(QAT) ≈ real(ED) atol=1e-10
            @test imag(QAT) ≈ imag(ED) atol=1e-10
        end
    end

    # ---- Sanity: equal-position autocorrelator at t=0 is 1 ------------------
    @testset "(σᶻ_i)² = 1 sanity" begin
        N, J, h = 12, 1.0, 0.7
        for i in (1, 4, 8, 12)
            v = QAtlas.fetch(
                TFIM(; J=J, h=h), DynamicalCorrelation(:z, :z), OBC(; N=N); i=i, j=i, t=0.0
            )
            @test real(v) ≈ 1.0 atol=1e-10
            @test abs(imag(v)) < 1e-10
        end
    end

    # ---- spreading convenience function ------------------------------------
    @testset "sz_sz_spreading shape and consistency" begin
        N, J, h = 10, 1.0, 0.5
        center = 5
        times = [0.0, 0.5, 1.0]
        C = QAtlas.fetch(
            TFIM(; J=J, h=h),
            LightconeSpinCorrelation(:z, :z),
            OBC(; N=N);
            center=center,
            times=times,
        )
        @test size(C) == (3, N)

        # Each entry should match an individual call.
        for (it, t) in enumerate(times), ix in 1:N
            ref = QAtlas.fetch(
                TFIM(; J=J, h=h),
                DynamicalCorrelation(:z, :z),
                OBC(; N=N);
                i=ix,
                j=center,
                t=t,
            )
            @test C[it, ix] ≈ ref atol=1e-10
        end
    end
end
