# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/XXZ/test_xxz1d_grand_canonical.jl
#
# Grand-canonical (`mu`) extension of XXZ1D OBC Energy / MagnetizationZ:
#
#   ⟨A⟩_{β,μ} = Tr(A exp(-β(Ĥ - μ N̂))) / Tr(exp(-β(Ĥ - μ N̂))),
#   N̂ = Ŝᶻ_tot = Σᵢ Sᶻᵢ = (1/2) Σᵢ σᶻᵢ   (spin-½ units).
#
# Independent verification route: a fully self-contained dense-ED built here
# from bare kron Pauli strings (NO QAtlas internals), plus analytic limits
# (μ = 0 anchor, μ ↔ -μ symmetry, μ → ∞ saturation). Reference for grand-
# canonical TPQ (gTPQ, Hyuga et al., PRB 90, 121110(R) (2014)). Refs #737.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra

# --- Independent dense-ED grand-canonical helpers (do NOT use QAtlas internals) ---

const _SXg = ComplexF64[0 1; 1 0]
const _SYg = ComplexF64[0 -im; im 0]
const _SZg = ComplexF64[1 0; 0 -1]
const _I2g = ComplexF64[1 0; 0 1]

_g_pauli(N, i, σ) = reduce(kron, [j == i ? σ : _I2g for j in 1:N])
_g_pair(N, i, j, σ) = reduce(kron, [(k == i || k == j) ? σ : _I2g for k in 1:N])

function _g_xxz_H(N, J, Δ)
    D = 2^N
    H = zeros(ComplexF64, D, D)
    for i in 1:(N - 1)
        H .+= (J / 4) .* _g_pair(N, i, i + 1, _SXg)
        H .+= (J / 4) .* _g_pair(N, i, i + 1, _SYg)
        H .+= (J * Δ / 4) .* _g_pair(N, i, i + 1, _SZg)
    end
    return H
end

_g_Sz_tot(N) = 0.5 .* sum(_g_pauli(N, i, _SZg) for i in 1:N)   # N̂ = Σ Sᶻ = ½ Σ σᶻ
_g_Mz(N) = sum(_g_pauli(N, i, _SZg) for i in 1:N)              # Pauli Σ σᶻ

# ⟨A⟩ under the grand ensemble exp(-β(H - μ N̂))
function _g_expect(A, H, Sz, β, μ)
    F = eigen(Hermitian(H .- μ .* Sz))
    emin = minimum(F.values)
    ws = exp.(-β .* (F.values .- emin))
    ws ./= sum(ws)
    Ad = real.(diag(F.vectors' * A * F.vectors))
    return sum(ws .* Ad)
end

@testset "XXZ1D grand-canonical Energy/MagnetizationZ vs independent ED" begin
    # Exact to ED round-off: same operators, theoretically distinct route (bare
    # kron ED here vs QAtlas's own eigen kernel). atol tracks the dense-ED noise
    # floor, matching the canonical XXZ1D observable tests.
    for (J, Δ) in ((1.0, 0.5), (1.0, 1.0), (1.0, -0.7)), N in (4, 6)
        m = XXZ1D(; J=J, Δ=Δ)
        H = _g_xxz_H(N, J, Δ)
        Sz = _g_Sz_tot(N)
        Mz = _g_Mz(N)
        for β in (0.5, 2.0), μ in (0.3, 1.5, -0.8)
            E_ref = _g_expect(H, H, Sz, β, μ)
            Mz_ref = _g_expect(Mz, H, Sz, β, μ) / N
            @test QAtlas.fetch(m, Energy(), OBC(N); beta=β, mu=μ) ≈ E_ref atol = 1e-10
            @test QAtlas.fetch(m, MagnetizationZ(), OBC(N); beta=β, mu=μ) ≈ Mz_ref atol =
                1e-10
        end
    end
end

@testset "XXZ1D grand-canonical: μ = 0 anchor, μ↔-μ symmetry, saturation" begin
    for (J, Δ) in ((1.0, 0.5), (1.0, 1.0)), N in (4, 6)
        m = XXZ1D(; J=J, Δ=Δ)

        # μ = 0 reduces to the canonical quantities (default-mu path).
        for β in (0.5, 2.0)
            @test QAtlas.fetch(m, Energy(), OBC(N); beta=β, mu=0.0) ≈
                QAtlas.fetch(m, Energy(), OBC(N); beta=β) atol = 1e-12
            @test abs(QAtlas.fetch(m, MagnetizationZ(), OBC(N); beta=β, mu=0.0)) < 1e-10
        end

        # Spin-flip symmetry of the XXZ Hamiltonian ⇒ E even in μ, M_z odd in μ.
        β = 1.3
        for μ in (0.4, 1.1)
            @test QAtlas.fetch(m, Energy(), OBC(N); beta=β, mu=μ) ≈
                QAtlas.fetch(m, Energy(), OBC(N); beta=β, mu=(-μ)) atol = 1e-10
            @test QAtlas.fetch(m, MagnetizationZ(), OBC(N); beta=β, mu=μ) ≈
                -QAtlas.fetch(m, MagnetizationZ(), OBC(N); beta=β, mu=(-μ)) atol = 1e-10
        end

        # Monotone polarisation: M_z increases with μ at fixed β.
        β = 1.0
        mz = [
            QAtlas.fetch(m, MagnetizationZ(), OBC(N); beta=β, mu=μ) for
            μ in (0.0, 0.5, 1.0, 2.0)
        ]
        @test issorted(mz)
        @test all(diff(mz) .> 0)

        # Saturation (μ, β → large): the field selects the all-up state, so
        # ⟨σᶻ⟩/site → 1 and ⟨Ĥ⟩ → the all-up energy J·Δ·(N-1)/4 (only the σᶻσᶻ
        # bond term survives on |↑…↑⟩). Purely analytic, no ED.
        E_allup = J * Δ * (N - 1) / 4
        @test QAtlas.fetch(m, MagnetizationZ(), OBC(N); beta=50.0, mu=10.0) ≈ 1.0 atol =
            1e-8
        @test QAtlas.fetch(m, Energy(), OBC(N); beta=50.0, mu=10.0) ≈ E_allup atol = 1e-8
    end
end
