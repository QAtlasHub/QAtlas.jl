# =============================================================================
# Tests for the TFIM static (equal-time) σˣσˣ correlator module
# (`src/models/quantum/TFIM/TFIM_xx_static.jl`).
#
# Layers:
#   1. T = 0 self-consistency: the new `SpinCorrelation{:x,:x}` equals the
#      real part of the existing `DynamicalCorrelation{(:x, :x)}` evaluated at
#      `t = 0`.
#   2. Boundary conditions of σˣ algebra:
#        ⟨(σˣ_i)²⟩ = 1  (since (σˣ)² = I)  ⇒ static at i = j returns 1.
#   3. Connected formula:  C^c_{ij} = C_{ij} − ⟨σˣ_i⟩ ⟨σˣ_j⟩
#      with `⟨σˣ_i⟩` from `MagnetizationXLocal` (which uses the same
#      Majorana-covariance internals).
#   4. Connected at i = j: 1 − ⟨σˣ_i⟩².
#   5. ED comparison at small N for thermal β ∈ {1.0, 2.5, Inf} —
#      uses the dense ED helpers from `test/util/tfim_dense_ed.jl`.
#   6. Infinite proxy: same call path as OBC(N_proxy), so exact equality
#      is required.
# =============================================================================

using QAtlas, Test, LinearAlgebra

@testset "TFIM XX static correlator" begin
    @testset "T=0 self-consistency: dynamic at t=0 = static" begin
        for h in (0.5, 1.0, 1.5), N in (8, 12)
            model = TFIM(; J=1.0, h=h)
            for i in 2:(N - 1), j in i:(N - 1)
                v_static = QAtlas.fetch(model, SpinCorrelation(:x, :x), OBC(N); i=i, j=j)
                v_dynamic_re = real(
                    QAtlas.fetch(
                        model, DynamicalCorrelation(:x, :x), OBC(N); i=i, j=j, t=0.0
                    ),
                )
                @test v_static ≈ v_dynamic_re atol=1e-10
            end
        end
    end

    @testset "i = j returns ⟨(σˣ)²⟩ = 1" begin
        for h in (0.5, 1.5), β in (Inf, 1.0)
            model = TFIM(; J=1.0, h=h)
            v = QAtlas.fetch(model, SpinCorrelation(:x, :x), OBC(8); beta=β, i=4, j=4)
            @test v ≈ 1.0 atol=1e-10
        end
    end

    @testset "Connected = static - ⟨σˣ_i⟩⟨σˣ_j⟩" begin
        h, N, β = 0.7, 10, 2.0
        model = TFIM(; J=1.0, h=h)
        mx_local = QAtlas.fetch(model, MagnetizationXLocal(), OBC(N); beta=β)
        for i in 3:(N - 2), j in (i + 1):(N - 2)
            v_st = QAtlas.fetch(model, SpinCorrelation(:x, :x), OBC(N); beta=β, i=i, j=j)
            v_cn = QAtlas.fetch(
                model, ConnectedSpinCorrelation(:x, :x), OBC(N); beta=β, i=i, j=j
            )
            @test v_cn ≈ v_st - mx_local[i] * mx_local[j] atol=1e-10
        end
    end

    @testset "Connected at i = j: 1 - ⟨σˣ⟩²" begin
        model = TFIM(; J=1.0, h=0.7)
        N, β = 8, 1.0
        mx_local = QAtlas.fetch(model, MagnetizationXLocal(), OBC(N); beta=β)
        for i in 2:(N - 1)
            v_cn = QAtlas.fetch(
                model, ConnectedSpinCorrelation(:x, :x), OBC(N); beta=β, i=i, j=i
            )
            @test v_cn ≈ 1.0 - mx_local[i]^2 atol=1e-10
        end
    end

    @testset "ED comparison at small N" begin
        # Uses helpers from test/util/tfim_dense_ed.jl
        # (_build_tfim_dense, _op_site, _SX) loaded once via runtests.jl.
        N, J, h = 4, 1.0, 0.7
        H = _build_tfim_dense(N, J, h)
        E, V = eigen(H)
        for β in (1.0, 2.5, Inf)
            ws = if isinf(β)
                [k == 1 ? 1.0 : 0.0 for k in 1:length(E)]
            else
                shifted = exp.(-β .* (E .- E[1]))
                shifted ./ sum(shifted)
            end
            ρ = V * Diagonal(ws) * V'
            for i in 1:N, j in i:N
                op = _op_site(_SX, i, N) * _op_site(_SX, j, N)
                ed_val = real(tr(ρ * op))
                qa_val = QAtlas.fetch(
                    TFIM(; J=J, h=h), SpinCorrelation(:x, :x), OBC(N); beta=β, i=i, j=j
                )
                @test qa_val ≈ ed_val atol=1e-10
            end
        end
    end

    @testset "Infinite proxy" begin
        # Same call path internally, so exact equality is the right test.
        h, β, i, j = 0.7, 2.0, 40, 50
        model = TFIM(; J=1.0, h=h)
        v_inf = QAtlas.fetch(
            model, SpinCorrelation(:x, :x), Infinite(); beta=β, i=i, j=j, N_proxy=80
        )
        v_obc = QAtlas.fetch(model, SpinCorrelation(:x, :x), OBC(80); beta=β, i=i, j=j)
        @test v_inf == v_obc
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "TFIM XX static — verification cards" begin
    let J = 1.0, h = 1.0, N = 6, i = 2, j = 4
        F = LinearAlgebra.eigen(_build_tfim_dense(N, J, h))
        ψ = F.vectors[:, 1]
        σx = ComplexF64[0 1; 1 0]
        xx_ed = real(ψ' * (_op_site(σx, i, N) * (_op_site(σx, j, N) * ψ)))
        verify(
            TFIM(; J=J, h=h),
            XXCorrelation(; mode=:static),
            OBC(N);
            route=:ed_finite_size,
            fetch_kw=(; i=i, j=j, beta=Inf),
            independent=xx_ed,
            agree_within=1e-8,
            refs=["Direct OBC dense-ED ⟨σx_i σx_j⟩ via _build_tfim_dense GS"],
        )
    end
end
