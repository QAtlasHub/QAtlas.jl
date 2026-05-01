using QAtlas, Test, LinearAlgebra

# Local σʸ Pauli matrix (the test util's `_SX, _SZ` consts in
# `tfim_dense_ed.jl` do not include `_SY`; we define it here so as not
# to perturb the shared util.).
const _SY = ComplexF64[0 -im; im 0]

@testset "TFIM YY OBC: static / connected / dynamic correlators" begin
    @testset "ED comparison: static ⟨σʸ_i σʸ_j⟩_β" begin
        for h in (0.5, 1.0, 1.5)
            N = 4
            for β in (Inf, 1.5, 0.5)
                H = _build_tfim_dense(N, 1.0, h)
                E, V = eigen(Hermitian(H))
                ws = if isinf(β)
                    [k == 1 ? 1.0 : 0.0 for k in 1:length(E)]
                else
                    w = exp.(-β .* (E .- E[1]))
                    w ./ sum(w)
                end
                ρ = V * (Diagonal(ComplexF64.(ws))) * V'
                for i in 1:N, j in i:N
                    op = _op_site(_SY, i, N) * _op_site(_SY, j, N)
                    ed_val = real(tr(ρ * op))
                    qa_val = QAtlas.fetch(
                        TFIM(; J=1.0, h=h),
                        YYCorrelation{:static}(),
                        OBC(N);
                        beta=β,
                        i=i,
                        j=j,
                    )
                    @test qa_val ≈ ed_val atol = 1e-10
                end
            end
        end
    end

    @testset "i = j returns ⟨(σʸ)²⟩ = 1" begin
        for h in (0.5, 1.0, 1.5), β in (Inf, 1.0)
            v = QAtlas.fetch(
                TFIM(; J=1.0, h=h), YYCorrelation{:static}(), OBC(8); beta=β, i=4, j=4
            )
            @test v ≈ 1.0 atol = 1e-12
        end
    end

    @testset "Connected = static (since ⟨σʸ⟩ = 0 by parity)" begin
        for h in (0.5, 1.0, 1.5), β in (1.0, Inf)
            for i in 2:7, j in i:7
                v_st = QAtlas.fetch(
                    TFIM(; J=1.0, h=h), YYCorrelation{:static}(), OBC(8); beta=β, i=i, j=j
                )
                v_cn = QAtlas.fetch(
                    TFIM(; J=1.0, h=h),
                    YYCorrelation{:connected}(),
                    OBC(8);
                    beta=β,
                    i=i,
                    j=j,
                )
                @test v_st ≈ v_cn atol = 1e-12
            end
        end
    end

    @testset "Dynamic at t = 0 reduces to static (real)" begin
        h, N, β = 0.7, 8, Inf
        model = TFIM(; J=1.0, h=h)
        for i in 2:(N - 1), j in i:(N - 1)
            c_dyn = QAtlas.fetch(
                model, YYCorrelation{:dynamic}(), OBC(N); i=i, j=j, t=0.0, beta=β
            )
            c_st = QAtlas.fetch(model, YYCorrelation{:static}(), OBC(N); i=i, j=j, beta=β)
            @test imag(c_dyn) ≈ 0 atol = 1e-12
            @test real(c_dyn) ≈ c_st atol = 1e-10
        end
    end

    @testset "Re even / Im odd in t (Hermitian autocorrelation)" begin
        h, N, β = 0.7, 8, 1.0
        model = TFIM(; J=1.0, h=h)
        for i in (3, 5), t in (0.4, 1.5)
            c_pos = QAtlas.fetch(
                model, YYCorrelation{:dynamic}(), OBC(N); i=i, j=i, t=t, beta=β
            )
            c_neg = QAtlas.fetch(
                model, YYCorrelation{:dynamic}(), OBC(N); i=i, j=i, t=(-t), beta=β
            )
            @test real(c_pos) ≈ real(c_neg) atol = 1e-10
            @test imag(c_pos) ≈ -imag(c_neg) atol = 1e-10
        end
    end
end

@testset "TFIM MagnetizationY OBC: identically zero" begin
    for h in (0.3, 0.7, 1.0, 1.5), β in (Inf, 0.5, 2.0)
        v = QAtlas.fetch(TFIM(; J=1.0, h=h), MagnetizationY(), OBC(6); beta=β)
        @test v == 0.0
    end
end

@testset "TFIM SusceptibilityYY OBC: ED variance comparison" begin
    for h in (0.5, 1.0, 1.5), β in (0.5, 1.5)
        N = 4
        H = _build_tfim_dense(N, 1.0, h)
        E, V = eigen(Hermitian(H))
        ws = exp.(-β .* (E .- E[1]));
        ws ./= sum(ws)
        ρ = V * (Diagonal(ComplexF64.(ws))) * V'
        M_y = sum(_op_site(_SY, k, N) for k in 1:N)
        M2 = real(tr(ρ * (M_y * M_y)))
        M1 = real(tr(ρ * M_y))
        χ_ed = β * (M2 - M1^2) / N
        χ_qa = QAtlas.fetch(TFIM(; J=1.0, h=h), SusceptibilityYY(), OBC(N); beta=β)
        @test χ_qa ≈ χ_ed atol = 1e-10
    end
end

@testset "SU(2)-broken TFIM: χ_yy ≠ χ_xx and ≠ χ_zz" begin
    # Sanity: TFIM is U(1)-broken (h breaks rotation), so the three
    # axis susceptibilities are distinct.  This catches an accidental
    # identity bug where one of (XX, YY, ZZ) wraps another.
    h, N, β = 0.7, 8, 1.0
    model = TFIM(; J=1.0, h=h)
    χ_xx = QAtlas.fetch(model, SusceptibilityXX(), OBC(N); beta=β)
    χ_yy = QAtlas.fetch(model, SusceptibilityYY(), OBC(N); beta=β)
    χ_zz = QAtlas.fetch(model, SusceptibilityZZ(), OBC(N); beta=β)
    @test !isapprox(χ_xx, χ_yy; atol=1e-3)
    @test !isapprox(χ_yy, χ_zz; atol=1e-3)
    @test χ_yy > 0
end
