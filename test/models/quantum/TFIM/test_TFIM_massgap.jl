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
