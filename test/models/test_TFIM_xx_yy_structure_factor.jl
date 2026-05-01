using QAtlas, Test, LinearAlgebra

const _SY = ComplexF64[0 -im; im 0]

# Independent ED structure factor: build the dense H, partial-trace the
# thermal state, and compute `S_αα(q) = (1/N) Σ_{ij} e^{-iq(i-j)} ⟨σᵅ_i σᵅ_j⟩`.
function _ed_structure_factor(σ_local::AbstractMatrix, N::Int, J::Real, h::Real, β::Real, q::Real)
    H = _build_tfim_dense(N, J, h)
    E, V = eigen(Hermitian(H))
    ws = exp.(-β .* (E .- E[1])); ws ./= sum(ws)
    ρ = V * (Diagonal(ComplexF64.(ws))) * V'
    s = 0.0 + 0.0im
    for i in 1:N, j in 1:N
        op = _op_site(σ_local, i, N) * _op_site(σ_local, j, N)
        s += exp(-im * q * (i - j)) * tr(ρ * op)
    end
    return real(s) / N
end

@testset "TFIM XXStructureFactor OBC: ED comparison" begin
    for h in (0.5, 1.0, 1.5), β in (0.5, 1.5)
        N = 4
        for q in (0.0, π / 2, π)
            ed_val = _ed_structure_factor(_SX, N, 1.0, h, β, q)
            qa_val = QAtlas.fetch(
                TFIM(; J=1.0, h=h), XXStructureFactor(), OBC(N); beta=β, q=q
            )
            @test qa_val ≈ ed_val atol = 1e-10
        end
    end
end

@testset "TFIM YYStructureFactor OBC: ED comparison" begin
    for h in (0.5, 1.0, 1.5), β in (0.5, 1.5)
        N = 4
        for q in (0.0, π / 2, π)
            ed_val = _ed_structure_factor(_SY, N, 1.0, h, β, q)
            qa_val = QAtlas.fetch(
                TFIM(; J=1.0, h=h), YYStructureFactor(), OBC(N); beta=β, q=q
            )
            @test qa_val ≈ ed_val atol = 1e-10
        end
    end
end

@testset "Bounds: S_αα(q) ≥ 0  (positive definite at any q)" begin
    # S_αα is the expectation of the projected operator |M_α(q)|², so
    # it is non-negative by Cauchy-Schwarz.
    for h in (0.5, 1.0, 1.5), β in (Inf, 1.0)
        for q in (0.0, π / 3, π / 2, 2π / 3, π)
            S_xx = QAtlas.fetch(
                TFIM(; J=1.0, h=h), XXStructureFactor(), OBC(8); beta=β, q=q
            )
            S_yy = QAtlas.fetch(
                TFIM(; J=1.0, h=h), YYStructureFactor(), OBC(8); beta=β, q=q
            )
            @test S_xx ≥ -1e-10
            @test S_yy ≥ -1e-10
        end
    end
end

@testset "Sum rule: ∫dq S_αα(q)/(2π) = ⟨(σᵅ)²⟩ = 1 (per-site, OBC bulk)" begin
    # In the bulk of an OBC chain Σᵢⱼ ⟨σᵅᵢσᵅⱼ⟩ = N (diagonal) +
    # off-diagonal cross-terms.  The diagonal gives ⟨(σᵅ)²⟩ = 1 per
    # site, and the off-diagonal averages out under q-integration.
    # ∫₀^{2π} S(q) dq / (2π) = (1/N) Σᵢⱼ δᵢⱼ ⟨...⟩ = (1/N) · N = 1.
    for h in (0.5, 1.5), β in (1.0,)
        N = 8
        # 12-point uniform grid over [0, 2π]
        qs = range(0; step=2π / 12, length=12)
        sxx_avg = sum(
            QAtlas.fetch(TFIM(; J=1.0, h=h), XXStructureFactor(), OBC(N); beta=β, q=q)
            for q in qs
        ) / length(qs)
        syy_avg = sum(
            QAtlas.fetch(TFIM(; J=1.0, h=h), YYStructureFactor(), OBC(N); beta=β, q=q)
            for q in qs
        ) / length(qs)
        @test sxx_avg ≈ 1.0 atol = 1e-10
        @test syy_avg ≈ 1.0 atol = 1e-10
    end
end

@testset "Infinite proxy = OBC at matching N_proxy" begin
    h, β, q = 0.7, 1.0, π / 3
    for N_proxy in (40, 80)
        S_inf = QAtlas.fetch(
            TFIM(; J=1.0, h=h), XXStructureFactor(), Infinite();
            beta=β, q=q, N_proxy=N_proxy,
        )
        S_obc = QAtlas.fetch(
            TFIM(; J=1.0, h=h), XXStructureFactor(), OBC(N_proxy); beta=β, q=q
        )
        @test S_inf == S_obc
    end
end
