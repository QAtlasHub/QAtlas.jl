# =============================================================================
# Tests for TFIM XXStructureFactor / YYStructureFactor on OBC.
#
# Migrated from pure-legacy @test pattern to verify()-first (PR #449 phase B,
# zero-legacy end-state). The ED cross-checks now run through verify() so
# each (J, h, N, β, q) point becomes a structural INVENTORY card via the
# :ed_finite_size route. The structural-property @testsets (positive
# definiteness, q-integral sum rule) stay as raw @test because they are
# multi-point invariants of the implementation rather than single-value
# hub pins — verify() is scalar by design.
# =============================================================================

using QAtlas, Test, LinearAlgebra

const _SY = ComplexF64[0 -im; im 0]

# Independent ED structure factor: build the dense H, partial-trace the
# thermal state, and compute `S_αα(q) = (1/N) Σ_{ij} e^{-iq(i-j)} ⟨σᵅ_i σᵅ_j⟩`.
function _ed_structure_factor(
    σ_local::AbstractMatrix, N::Int, J::Real, h::Real, β::Real, q::Real
)
    H = _build_tfim_dense(N, J, h)
    E, V = eigen(Hermitian(H))
    ws = exp.(-β .* (E .- E[1]))
    ws ./= sum(ws)
    ρ = V * (Diagonal(ComplexF64.(ws))) * V'
    s = 0.0 + 0.0im
    for i in 1:N, j in 1:N
        op = _op_site(σ_local, i, N) * _op_site(σ_local, j, N)
        s += exp(-im * q * (i - j)) * tr(ρ * op)
    end
    return real(s) / N
end

@testset "TFIM XXStructureFactor OBC — ED cross-check (verify migration)" begin
    for h in (0.5, 1.0, 1.5), β in (0.5, 1.5)
        N = 4
        for q in (0.0, π / 2, π)
            verify(
                TFIM(; J=1.0, h=h),
                XXStructureFactor(),
                OBC(N);
                route=:ed_finite_size,
                independent=_ed_structure_factor(_SX, N, 1.0, h, β, q),
                agree_within=1e-10,
                at=["h=$(h)", "β=$(β)", "q=$(round(q; digits=4))"],
                refs=[
                    "ED black-box: _build_tfim_dense + dense thermal trace ⟨σˣᵢσˣⱼ⟩ Fourier-summed at q",
                ],
                fetch_kw=(; beta=β, q=q),
            )
        end
    end
end

@testset "TFIM YYStructureFactor OBC — ED cross-check (verify migration)" begin
    for h in (0.5, 1.0, 1.5), β in (0.5, 1.5)
        N = 4
        for q in (0.0, π / 2, π)
            verify(
                TFIM(; J=1.0, h=h),
                YYStructureFactor(),
                OBC(N);
                route=:ed_finite_size,
                independent=_ed_structure_factor(_SY, N, 1.0, h, β, q),
                agree_within=1e-10,
                at=["h=$(h)", "β=$(β)", "q=$(round(q; digits=4))"],
                refs=[
                    "ED black-box: _build_tfim_dense + dense thermal trace ⟨σʸᵢσʸⱼ⟩ Fourier-summed at q",
                ],
                fetch_kw=(; beta=β, q=q),
            )
        end
    end
end

@testset "TFIM XXStructureFactor Infinite — N_proxy pass-through (verify migration)" begin
    h, β, q = 0.7, 1.0, π / 3
    for N_proxy in (40, 80)
        verify(
            TFIM(; J=1.0, h=h),
            XXStructureFactor(),
            Infinite();
            route=:delegation_invariant,
            independent=QAtlas.fetch(
                TFIM(; J=1.0, h=h), XXStructureFactor(), OBC(N_proxy); beta=β, q=q
            ),
            agree_within=0.0,
            at=["N_proxy=$(N_proxy)"],
            refs=[
                "Infinite() proxy returns exactly the OBC(N_proxy) value by construction (pass-through identity, == compared)",
            ],
            fetch_kw=(; beta=β, q=q, N_proxy=N_proxy),
        )
    end
end

@testset "TFIM YYStructureFactor Infinite — N_proxy pass-through (verify migration)" begin
    h, β, q = 0.7, 1.0, π / 3
    for N_proxy in (40, 80)
        verify(
            TFIM(; J=1.0, h=h),
            YYStructureFactor(),
            Infinite();
            route=:delegation_invariant,
            independent=QAtlas.fetch(
                TFIM(; J=1.0, h=h), YYStructureFactor(), OBC(N_proxy); beta=β, q=q
            ),
            agree_within=0.0,
            at=["N_proxy=$(N_proxy)"],
            refs=[
                "Infinite() proxy returns exactly the OBC(N_proxy) value by construction (pass-through identity, == compared)",
            ],
            fetch_kw=(; beta=β, q=q, N_proxy=N_proxy),
        )
    end
end

# ── Structural invariants (NOT single-value hub pins; kept as raw @test) ───
#
# verify() is scalar-by-design (one fetch ↔ one independent). The blocks
# below assert multi-point properties of the implementation:
#   (1) Positive definiteness S_αα(q) ≥ 0 — Cauchy-Schwarz on |M_α(q)|²
#   (2) Sum rule ∫dq S(q)/(2π) = ⟨(σᵅ)²⟩ = 1 — diagonal piece per-site
# Neither maps to a single (model, quantity, bc) hub. The verify-card
# coverage above already pins the per-(q,β,h) value; these guards just
# enforce derived global properties.

@testset "TFIM SF — positive definiteness S_αα(q) ≥ 0 (structural)" begin
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

@testset "TFIM SF — q-integral sum rule ∫dq S(q)/(2π) = 1 (structural)" begin
    # ∫₀^{2π} S(q) dq / (2π) = (1/N) Σᵢⱼ δᵢⱼ ⟨(σᵅ)²⟩ = (1/N) · N · 1 = 1
    for h in (0.5, 1.5), β in (1.0,)
        N = 8
        qs = range(0; step=2π / 12, length=12)
        sxx_avg =
            sum(
                QAtlas.fetch(TFIM(; J=1.0, h=h), XXStructureFactor(), OBC(N); beta=β, q=q)
                for q in qs
            ) / length(qs)
        syy_avg =
            sum(
                QAtlas.fetch(TFIM(; J=1.0, h=h), YYStructureFactor(), OBC(N); beta=β, q=q)
                for q in qs
            ) / length(qs)
        @test sxx_avg ≈ 1.0 atol = 1e-10
        @test syy_avg ≈ 1.0 atol = 1e-10
    end
end
