# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: ChernSimons3D — Sugawara central charge for SU(N)_k.
#
# Verifies:
#   * Canonical c(SU(2)_k) = 3k/(k+2) for k = 1, 2, 3, 4
#   * c(SU(3)_1) = 2, c(SU(N)_1) = N - 1 (free-fermion realisation)
#   * Large-k limit c → N² - 1 (dim of SU(N))
#   * DomainError on N < 2 or k < 1
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "ChernSimons3D — SU(2)_k central charge equals 3k/(k+2)" begin
    for k in 1:6
        c = QAtlas.fetch(ChernSimons3D(; N=2, k=k), CentralCharge(), Infinite())
        @test c == Rational(3 * k, k + 2)
    end
end

@testset "ChernSimons3D — SU(N)_1 central charge equals N - 1" begin
    for N in 2:6
        c = QAtlas.fetch(ChernSimons3D(; N=N, k=1), CentralCharge(), Infinite())
        @test c == N - 1
    end
end

@testset "ChernSimons3D — Large-k limit c → N² - 1" begin
    N = 3
    for k in (10, 100, 1000)
        c = QAtlas.fetch(ChernSimons3D(; N=N, k=k), CentralCharge(), Infinite())
        @test c < N^2 - 1
        @test float(c) ≈ N^2 - 1 atol = (N^2 - 1) * N / (k + N)   # 1st-order correction
    end
end

@testset "ChernSimons3D — DomainError on invalid (N, k)" begin
    @test_throws DomainError ChernSimons3D(; N=1, k=1)
    @test_throws DomainError ChernSimons3D(; N=2, k=0)
    @test_throws DomainError ChernSimons3D(; N=2, k=-1)
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 2: PartitionFunction Z(S³) (Witten 1989) — modular S_{0,0}
# ─────────────────────────────────────────────────────────────────────────────

@testset "ChernSimons3D — PartitionFunction Z(S³) Witten 1989 (Phase 2)" begin
    # SU(2)_1: Z = √(2/3) sin(π/3) = 1/√2
    z21 = QAtlas.fetch(ChernSimons3D(; N=2, k=1), PartitionFunction(), Infinite())
    @test z21 ≈ 1 / sqrt(2)
    @test z21 ≈ 0.7071067811865476

    # SU(2)_2: Z = √(2/4) sin(π/4) = 1/2
    z22 = QAtlas.fetch(ChernSimons3D(; N=2, k=2), PartitionFunction(), Infinite())
    @test z22 ≈ 0.5

    # SU(2)_3: Z = √(2/5) sin(π/5)
    z23 = QAtlas.fetch(ChernSimons3D(; N=2, k=3), PartitionFunction(), Infinite())
    @test z23 ≈ sqrt(2 / 5) * sin(π / 5)
    @test z23 ≈ 0.3717480344601845

    # General SU(2)_k vs the SU(2) shortcut √(2/(k+2)) sin(π/(k+2)).
    for k in 3:6
        z_general = QAtlas.fetch(ChernSimons3D(; N=2, k=k), PartitionFunction(), Infinite())
        z_su2 = sqrt(2 / (k + 2)) * sin(π / (k + 2))
        @test z_general ≈ z_su2
    end

    # SU(3)_1: positive roots (j<l) are (1,2),(1,3),(2,3) with p = k+N = 4.
    # ∏ 2 sin(π m/4) for m=1,2,1 = (2 sin π/4)(2 sin π/2)(2 sin π/4)
    #                            = √2 · 2 · √2 = 4.
    # Prefactor: N^{-1/2} p^{-(N-1)/2} = (1/√3)(1/4).
    # ⇒ Z = (1/√3)(1/4)·4 = 1/√3.
    z31 = QAtlas.fetch(ChernSimons3D(; N=3, k=1), PartitionFunction(), Infinite())
    @test z31 ≈ 1 / sqrt(3)
    @test z31 ≈ 0.5773502691896258
end

@testset "ChernSimons3D — level-rank duality √N·Z(SU(N)_k) = √k·Z(SU(k)_N) (Phase 2)" begin
    # Level-rank duality on S³ in the bare SU(N)_k = S_{0,0} normalisation
    # holds up to a √(N/k) factor (cf. Naculich-Riggs-Schnitzer 1990,
    # Mlawer-Naculich-Riggs-Schnitzer 1991): the dual partition
    # functions agree after rescaling by √rank, which removes the
    # SU(N) ↔ U(N) normalisation mismatch.  This catches N↔k swap
    # symmetry of the product-over-positive-roots formula, which the
    # SU(2)/SU(3) numerical specialisations alone do not exercise.
    for (N, k) in [(2, 3), (2, 4), (3, 5), (4, 3), (5, 7)]
        z_Nk = QAtlas.fetch(ChernSimons3D(; N=N, k=k), PartitionFunction(), Infinite())
        z_kN = QAtlas.fetch(ChernSimons3D(; N=k, k=N), PartitionFunction(), Infinite())
        @test sqrt(N) * z_Nk ≈ sqrt(k) * z_kN
    end
end

@testset "ChernSimons3D — PartitionFunction rejects N<2 / k<1 via kwargs (Phase 2)" begin
    # Constructor rejects invalid (N,k) up-front, so to exercise the
    # in-fetch guard we override via kwargs on a valid base model.
    m = ChernSimons3D(; N=2, k=1)
    @test_throws DomainError QAtlas.fetch(m, PartitionFunction(), Infinite(); N=1)
    @test_throws DomainError QAtlas.fetch(m, PartitionFunction(), Infinite(); k=0)
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "ChernSimons3D — verification cards" begin
    # SU(2)_k WZW central charge c = 3k/(k+2) (Knizhnik-Zamolodchikov).
    for k in (1, 2, 3, 6)
        verify(
            ChernSimons3D(; N=2, k=k),
            CentralCharge(),
            Infinite();
            route=:second_closed_form,
            independent=3k // (k + 2),
            agree_within=1e-10,
            refs=["SU(2)_k WZW: c = 3k/(k+2) (Knizhnik-Zamolodchikov 1984)"],
        )
    end
end

# ── additional verification cards (#381 batch) ─────────────────────────────
@testset "ChernSimons3D — PartitionFunction Z(S³) closed forms (#381 batch)" begin
    # Witten 1989 / Verlinde S_{0,0}: Z(S³; SU(N)_k) = N^{-1/2} (k+N)^{-(N-1)/2}
    # ∏_{1≤j<l≤N} 2 sin(π(l-j)/(k+N)). Reproduced low-rank closed forms below.
    # SU(2)_1: Z = 1/√2
    verify(
        ChernSimons3D(; N=2, k=1),
        PartitionFunction(),
        Infinite();
        route=:second_closed_form,
        independent=1/sqrt(2),
        agree_within=1e-14,
        refs=["Witten 1989: Z(S³; SU(2)_1) = sqrt(2/3)·sin(π/3) = 1/√2"],
    )
    # SU(2)_2: Z = 1/2
    verify(
        ChernSimons3D(; N=2, k=2),
        PartitionFunction(),
        Infinite();
        route=:second_closed_form,
        independent=0.5,
        agree_within=1e-14,
        refs=["Witten 1989: Z(S³; SU(2)_2) = sqrt(2/4)·sin(π/4) = 1/2"],
    )
    # SU(2)_3: Z = √(2/5) · sin(π/5)
    verify(
        ChernSimons3D(; N=2, k=3),
        PartitionFunction(),
        Infinite();
        route=:second_closed_form,
        independent=sqrt(2/5) * sin(π/5),
        agree_within=1e-14,
        refs=["Witten 1989: Z(S³; SU(2)_3) = sqrt(2/5)·sin(π/5)"],
    )
    # SU(3)_1: Z = 1/√3
    verify(
        ChernSimons3D(; N=3, k=1),
        PartitionFunction(),
        Infinite();
        route=:second_closed_form,
        independent=1/sqrt(3),
        agree_within=1e-14,
        refs=["Witten 1989: Z(S³; SU(3)_1) = 1/√3 (Verlinde S₀₀ closed form)"],
    )
end

