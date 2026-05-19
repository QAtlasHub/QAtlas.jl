using Test
using QAtlas
using LinearAlgebra: norm

# Time-domain symmetry identities for TFIM dynamic correlators.
#
# Universal identities for *Hermitian* operators A, B at thermal
# equilibrium with a *real* Hamiltonian:
#
#   (1)  Re ⟨A(t) A⟩  is even in t      (autocorrelator reality)
#        Im ⟨A(t) A⟩  is odd  in t
#
#   (2)  ⟨A(t) B⟩^* = ⟨B(-t) A⟩         (Hermiticity + cyclic +
#                                         time-translation in equilibrium)
#
#   (3)  ⟨A(0) B⟩ = ⟨A B⟩ (real)        (t = 0 reduction to static)
#
#   (4)  Lieb-Robinson lightcone: |⟨A_x(t) B_y⟩| decays exponentially
#        when `v_max · t < |x − y|` with `v_max = max_k |dΛ/dk|` =
#        `2 max(J, h)` for TFIM.
#
# These are independent of the BdG / Pfaffian implementation details
# in TFIM_dynamics.jl — they probe whether the dynamic σᶻσᶻ / σˣσˣ
# correlators respect the basic Hermitian time-evolution structure.

# ────────────────────────────────────────────────────────────────────
# (1) Reflection symmetry of the auto-correlator
# ────────────────────────────────────────────────────────────────────

@testset "Auto-correlator: Re even in t, Im odd in t" begin
    for h in (0.5, 1.0, 1.5), β in (Inf, 2.0)
        model = TFIM(; J=1.0, h=h)
        N = 8
        for i in (3, 4, 5), t in (0.4, 1.2, 2.5)
            c_pos = QAtlas.fetch(
                model, ZZCorrelation{:dynamic}(), OBC(N); i=i, j=i, t=t, beta=β
            )
            c_neg = QAtlas.fetch(
                model, ZZCorrelation{:dynamic}(), OBC(N); i=i, j=i, t=(-t), beta=β
            )
            @test real(c_pos) ≈ real(c_neg) atol = 1e-10
            @test imag(c_pos) ≈ -imag(c_neg) atol = 1e-10
        end
    end
end

@testset "XX auto-correlator: Re even, Im odd" begin
    h, N, β = 0.7, 8, Inf
    model = TFIM(; J=1.0, h=h)
    for i in (3, 4, 5), t in (0.4, 1.5)
        c_pos = QAtlas.fetch(
            model, XXCorrelation{:dynamic}(), OBC(N); i=i, j=i, t=t, beta=β
        )
        c_neg = QAtlas.fetch(
            model, XXCorrelation{:dynamic}(), OBC(N); i=i, j=i, t=(-t), beta=β
        )
        @test real(c_pos) ≈ real(c_neg) atol = 1e-10
        @test imag(c_pos) ≈ -imag(c_neg) atol = 1e-10
    end
end

# ────────────────────────────────────────────────────────────────────
# (2) Hermiticity + cyclic + time-translation
#     ⟨A(t) B⟩^* = ⟨B(-t) A⟩
# ────────────────────────────────────────────────────────────────────

@testset "Hermitian-equilibrium identity: ⟨σᶻ_i(t)σᶻ_j⟩* = ⟨σᶻ_j(-t)σᶻ_i⟩" begin
    h, N, β = 0.7, 8, Inf
    model = TFIM(; J=1.0, h=h)
    for (i, j) in ((3, 5), (2, 6), (3, 7)), t in (0.5, 1.7)
        c_ij_pos = QAtlas.fetch(
            model, ZZCorrelation{:dynamic}(), OBC(N); i=i, j=j, t=t, beta=β
        )
        c_ji_neg = QAtlas.fetch(
            model, ZZCorrelation{:dynamic}(), OBC(N); i=j, j=i, t=(-t), beta=β
        )
        @test conj(c_ij_pos) ≈ c_ji_neg atol = 1e-10
    end
end

@testset "Hermitian identity at finite β: ⟨σᶻ_i(t)σᶻ_j⟩* = ⟨σᶻ_j(-t)σᶻ_i⟩" begin
    # KMS-extension at finite β: equilibrium thermal expectation still
    # commutes with H, so the identity still holds.
    h, N, β = 0.5, 8, 1.5
    model = TFIM(; J=1.0, h=h)
    for (i, j) in ((3, 5), (4, 7)), t in (0.6, 1.4)
        c_ij_pos = QAtlas.fetch(
            model, ZZCorrelation{:dynamic}(), OBC(N); i=i, j=j, t=t, beta=β
        )
        c_ji_neg = QAtlas.fetch(
            model, ZZCorrelation{:dynamic}(), OBC(N); i=j, j=i, t=(-t), beta=β
        )
        @test conj(c_ij_pos) ≈ c_ji_neg atol = 1e-10
    end
end

# ────────────────────────────────────────────────────────────────────
# (3) t = 0 reduction:  dynamic(t=0) == static
# ────────────────────────────────────────────────────────────────────

@testset "Dynamic at t = 0 reduces to static correlator (ZZ)" begin
    for h in (0.5, 1.0, 1.5), β in (Inf, 2.0)
        model = TFIM(; J=1.0, h=h)
        N = 8
        for i in 2:(N - 1), j in i:(N - 1)
            c_dyn = QAtlas.fetch(
                model, ZZCorrelation{:dynamic}(), OBC(N); i=i, j=j, t=0.0, beta=β
            )
            c_stat = QAtlas.fetch(model, ZZCorrelation{:static}(), OBC(N); i=i, j=j, beta=β)
            @test imag(c_dyn) ≈ 0 atol = 1e-12
            @test real(c_dyn) ≈ c_stat atol = 1e-10
        end
    end
end

@testset "Dynamic at t = 0 reduces to static correlator (XX)" begin
    h, N, β = 0.7, 8, Inf
    model = TFIM(; J=1.0, h=h)
    for i in 2:(N - 1), j in i:(N - 1)
        c_dyn = QAtlas.fetch(
            model, XXCorrelation{:dynamic}(), OBC(N); i=i, j=j, t=0.0, beta=β
        )
        c_stat = QAtlas.fetch(model, XXCorrelation{:static}(), OBC(N); i=i, j=j, beta=β)
        @test imag(c_dyn) ≈ 0 atol = 1e-12
        @test real(c_dyn) ≈ c_stat atol = 1e-10
    end
end

# ────────────────────────────────────────────────────────────────────
# (4) Lieb-Robinson lightcone — connected-correlator change
# ────────────────────────────────────────────────────────────────────

@testset "Lieb-Robinson: |⟨σᶻ_i(t)σᶻ_j⟩ − ⟨σᶻ_i σᶻ_j⟩| small outside lightcone" begin
    # TFIM maximum group velocity for h ≤ J: v_max = 2 max(J, h) = 2J.
    # Lieb-Robinson bound applies to the *change* of the correlator:
    #   |⟨A_x(t) B_y⟩ − ⟨A_x B_y⟩| ≤ C exp(-(Δr − v_max t)/ξ)
    # The static piece itself need not be small (FM ordered phase has
    # long-range σᶻ correlation), so subtract it before bounding.
    J, h, N = 1.0, 0.5, 16
    model = TFIM(; J=J, h=h)
    v_max = 2 * max(J, h)
    i, j = 4, 12
    Δr = j - i  # = 8

    c_static = QAtlas.fetch(model, ZZCorrelation{:static}(), OBC(N); i=i, j=j, beta=Inf)

    # Just outside lightcone: t = Δr/(4 v_max) → exponential decay.
    t_out = Δr / (4 * v_max)
    c_out = abs(
        QAtlas.fetch(
            model, ZZCorrelation{:dynamic}(), OBC(N); i=i, j=j, t=t_out, beta=Inf
        ) - c_static,
    )

    # Well inside lightcone: t = 4 Δr / v_max → O(1) deviation.
    t_in = 4 * Δr / v_max
    c_in = abs(
        QAtlas.fetch(model, ZZCorrelation{:dynamic}(), OBC(N); i=i, j=j, t=t_in, beta=Inf) -
        c_static,
    )

    @test c_out < c_in    # outside-lightcone deviation is smaller
    @test c_out < 0.1     # exponentially small in `(Δr − v_max t)/ξ`
end

# ────────────────────────────────────────────────────────────────────
# (5) Heisenberg equation initial slope:
#     i ∂_t ⟨A(t) B⟩|_{t=0} = ⟨[H, A] B⟩  (operator commutator on σ_x)
# ────────────────────────────────────────────────────────────────────

@testset "Initial slope ↔ operator commutator (XX dynamic)" begin
    # Heisenberg equation: ∂_t A(t) = i [H, A(t)]; at t = 0:
    #   ∂_t ⟨A(t) B⟩|_{t=0} = i ⟨[H, A] B⟩.
    # For TFIM with A = σˣ_i, [H, A] generates a local operator string
    # supported on sites {i-1, i, i+1} only — so the initial slope of
    # ⟨σˣ_i(t) σˣ_j⟩ vanishes identically when |i − j| > 1 (locality
    # / Lieb-Robinson at t = 0).  We restrict to nearest-neighbour
    # pairs where the rhs commutator overlaps with σˣ_j.
    h, J, N, β = 0.7, 1.0, 8, Inf
    model = TFIM(; J=J, h=h)
    δ = 1e-4
    for i in (3, 5)
        j = i + 1   # nearest neighbour, well inside the bulk
        c_p = QAtlas.fetch(model, XXCorrelation{:dynamic}(), OBC(N); i=i, j=j, t=δ, beta=β)
        c_m = QAtlas.fetch(
            model, XXCorrelation{:dynamic}(), OBC(N); i=i, j=j, t=(-δ), beta=β
        )
        re_slope = (real(c_p) - real(c_m)) / (2δ)
        im_slope = (imag(c_p) - imag(c_m)) / (2δ)
        # Re part of ⟨A(t)B⟩ is even in t for this Hermitian-equilibrium
        # case, so its slope at t=0 vanishes; Im part is odd → its slope
        # is non-zero for nearest-neighbour pairs.  A Re ↔ Im swap
        # (sign error in the time-evolution code) would surface here.
        @test abs(re_slope) < 1e-4
        @test abs(im_slope) > 1e-4
    end
    # Locality check: at |i − j| = 3 the slope is identically 0.
    c_p = QAtlas.fetch(model, XXCorrelation{:dynamic}(), OBC(N); i=3, j=6, t=δ, beta=β)
    c_m = QAtlas.fetch(model, XXCorrelation{:dynamic}(), OBC(N); i=3, j=6, t=(-δ), beta=β)
    far_im_slope = (imag(c_p) - imag(c_m)) / (2δ)
    @test abs(far_im_slope) < 1e-10
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "TFIM dynamic symmetries — verification cards" begin
    # The dynamic-correlator symmetry tested here rests on the TFIM
    # Bogoliubov spectrum; its gap is the Pfeuty 1970 closed form
    # Δ = 2|h - J| (independent of src).
    for h in (0.5, 1.5, 2.0)
        verify(
            TFIM(; J=1.0, h=h),
            MassGap(),
            Infinite();
            route=:second_closed_form,
            independent=2 * abs(h - 1.0),
            agree_within=1e-10,
            refs=["Pfeuty 1970: Δ = 2|h - J| (Bogoliubov dispersion minimum)"],
        )
    end
end
