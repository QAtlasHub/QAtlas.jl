using QAtlas, Test

@testset "TFIM MassGap — infinite chain closed form" begin
    # Δ = 2|h − J|
    @test QAtlas.fetch(TFIM(; J=1.0, h=0.0), MassGap(), Infinite()) ≈ 2.0
    @test QAtlas.fetch(TFIM(; J=1.0, h=0.5), MassGap(), Infinite()) ≈ 1.0
    @test QAtlas.fetch(TFIM(; J=1.0, h=1.0), MassGap(), Infinite()) == 0.0  # critical
    @test QAtlas.fetch(TFIM(; J=1.0, h=2.0), MassGap(), Infinite()) ≈ 2.0
    @test QAtlas.fetch(TFIM(; J=0.5, h=1.5), MassGap(), Infinite()) ≈ 2.0
end

@testset "TFIM MassGap — OBC BdG" begin
    # Large-N OBC gap converges to the infinite-chain closed form
    # Δ = 2|h − J|.  The BdG matrix is cheap (N×N), so N=200 costs nothing.
    for (J, h) in ((1.0, 0.3), (1.0, 3.0), (2.0, 0.5))
        Δ_obc = QAtlas.fetch(TFIM(; J=J, h=h), MassGap(), OBC(200))
        Δ_inf = QAtlas.fetch(TFIM(; J=J, h=h), MassGap(), Infinite())
        @test Δ_obc ≈ Δ_inf rtol = 1e-4
    end

    # Gap is always positive (lowest positive BdG eigenvalue).
    @test QAtlas.fetch(TFIM(; J=1.0, h=0.5), MassGap(), OBC(32)) > 0

    # At the critical point h = J the OBC gap closes as ~ π J / N.
    gap_8 = QAtlas.fetch(TFIM(; J=1.0, h=1.0), MassGap(), OBC(8))
    gap_16 = QAtlas.fetch(TFIM(; J=1.0, h=1.0), MassGap(), OBC(16))
    gap_32 = QAtlas.fetch(TFIM(; J=1.0, h=1.0), MassGap(), OBC(32))
    @test gap_8 > gap_16 > gap_32 > 0
    # Asymptotically gap_2N / gap_N → 1/2 for the CFT finite-size gap.
    @test gap_32 / gap_16 < 0.7

    # J scaling: gap is linear in the overall energy scale.
    @test QAtlas.fetch(TFIM(; J=3.0, h=0.0), MassGap(), Infinite()) ≈ 6.0
    @test QAtlas.fetch(TFIM(; J=3.0, h=0.0), MassGap(), OBC(32)) ≈
        3 * QAtlas.fetch(TFIM(; J=1.0, h=0.0), MassGap(), OBC(32)) rtol = 1e-12
end

@testset "TFIM MassGap — legacy Symbol dispatch" begin
    # :mass_gap + aliases (:gap, :MassGap, :excitation_gap) all route to the
    # analytic value. The shim emits one `@info` per (model, quantity) pair;
    # `@test_logs` catches those so the deprecation notices stay out of CI.
    @test (@test_logs (:info, r"symbol-dispatch") QAtlas.fetch(
        :TFIM, :mass_gap, Infinite(); J=1.0, h=2.0
    )) ≈ 2.0
    @test (@test_logs (:info, r"symbol-dispatch") QAtlas.fetch(
        :TFIM, :gap, Infinite(); J=1.0, h=2.0
    )) ≈ 2.0
    @test (@test_logs (:info, r"symbol-dispatch") QAtlas.fetch(
        :TFIM, :MassGap, Infinite(); J=1.0, h=2.0
    )) ≈ 2.0
    @test (@test_logs (:info, r"symbol-dispatch") QAtlas.fetch(
        :TFIM, :excitation_gap, Infinite(); J=1.0, h=2.0
    )) ≈ 2.0

    # OBC Symbol dispatch with legacy N kwarg
    Δ_legacy = @test_logs (:info, r"symbol-dispatch") QAtlas.fetch(
        :TFIM, :mass_gap, OBC(); N=24, J=1.0, h=3.0
    )
    Δ_new = QAtlas.fetch(TFIM(; J=1.0, h=3.0), MassGap(), OBC(24))
    @test Δ_legacy ≈ Δ_new
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "TFIM MassGap — verification cards" begin
    # Pfeuty 1970: the TFIM mass gap is exactly 2|h - J| in the
    # thermodynamic limit (independent closed form from the Bogoliubov
    # dispersion min, not read from src).
    for (J, h) in ((1.0, 0.5), (1.0, 2.0), (1.0, 0.0), (2.0, 0.7))
        verify(
            TFIM(; J=J, h=h),
            MassGap(),
            Infinite();
            route=:second_closed_form,
            independent=2 * abs(h - J),
            agree_within=1e-10,
            refs=["Pfeuty 1970: Delta = 2|h - J| (Bogoliubov dispersion minimum)"],
        )
    end

    # Independent OBC dense-ED gap at small N (first excitation above GS)
    let J = 1.0, h = 2.0, N = 8
        sp = dense_spectrum(_build_tfim_dense(N, J, h))
        verify(
            TFIM(; J=J, h=h),
            MassGap(),
            OBC(N);
            route=:ed_finite_size,
            independent=sp[2] - sp[1],
            agree_within=1e-9,
            refs=["Direct OBC dense ED first excitation via _build_tfim_dense"],
        )
    end
end
