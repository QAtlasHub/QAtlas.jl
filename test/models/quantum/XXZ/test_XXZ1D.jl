using QAtlas, Test

@testset "XXZ1D — dispatch & construction" begin
    m = XXZ1D()
    @test m isa QAtlas.AbstractQAtlasModel
    @test m.J == 1.0
    @test m.Δ == 0.0

    m2 = XXZ1D(; J=2.5, Δ=0.7)
    @test m2.J == 2.5
    @test m2.Δ == 0.7

    # Symbol alias resolves to :XXZ1D canonical
    @test QAtlas.canonicalize_model(Val(:XXZ)) === :XXZ1D
    @test QAtlas.canonicalize_model(Val(:xxz)) === :XXZ1D
    @test QAtlas.canonicalize_model(Val(:xxz1d)) === :XXZ1D
end

@testset "XXZ1D — known closed-form points (J=1, Infinite)" begin
    # Δ = 0 (XX, free fermion): e₀/J = -1/π
    e_xx = QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.0), Energy(), Infinite())
    @test e_xx ≈ -1 / π atol = 1e-10

    # Δ = 1 (AF Heisenberg, Hulthén 1938): e₀/J = 1/4 - ln 2
    e_af = QAtlas.fetch(XXZ1D(; J=1.0, Δ=1.0), Energy(), Infinite())
    @test e_af ≈ 0.25 - log(2.0) atol = 1e-10

    # Δ = -1 (FM): e₀/J = -1/4
    e_fm = QAtlas.fetch(XXZ1D(; J=1.0, Δ=-1.0), Energy(), Infinite())
    @test e_fm ≈ -0.25 atol = 1e-14

    # GroundStateEnergyDensity alias returns the same value.
    e_gs = QAtlas.fetch(XXZ1D(; J=1.0, Δ=1.0), GroundStateEnergyDensity(), Infinite())
    @test e_gs ≈ e_af
end

@testset "XXZ1D — Energy at general -1 < Δ < 1 (Yang-Yang single integral)" begin
    # Yang-Yang single-integral form, evaluated by QuadGK.  Validation
    # touches the rational γ = π/3 point Δ = 1/2 where the literature
    # value is the elementary -3/8 (Yang-Yang II 1966), and a sweep of
    # generic γ across the interior of the gapless interval.
    @test QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.5), Energy(), Infinite()) ≈ -3 / 8 atol = 1e-10
    @test QAtlas.fetch(XXZ1D(; J=2.5, Δ=0.5), Energy(), Infinite()) ≈ -2.5 * 3 / 8 atol =
        1e-10

    # Smooth, monotone decreasing in Δ across (-1, 1) — more antiferro-
    # magnetic order means lower energy density.  Sample at non-canonical
    # points to exercise the integral branch end-to-end.
    Δs = -0.9:0.1:0.9
    es = [QAtlas.fetch(XXZ1D(; J=1, Δ=Δ), Energy(), Infinite()) for Δ in Δs]
    @test all(isfinite, es)
    @test all(diff(es) .< 0)

    # Boundary continuity: e₀ at Δ = ±0.99 sits between the Δ = ±1 closed
    # forms and the central XX value, with the right sign.
    e_near_p = QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.99), Energy(), Infinite())
    @test (0.25 - log(2.0)) ≤ e_near_p ≤ -1 / π
    e_near_m = QAtlas.fetch(XXZ1D(; J=1.0, Δ=-0.99), Energy(), Infinite())
    @test -1 / π ≤ e_near_m ≤ -0.25
end

@testset "XXZ1D — Luttinger parameter K" begin
    @test QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.0), LuttingerParameter(), Infinite()) ≈ 1.0 atol =
        1e-12
    @test QAtlas.fetch(XXZ1D(; J=1.0, Δ=1.0), LuttingerParameter(), Infinite()) ≈ 0.5 atol =
        1e-12

    # Monotone decreasing from Δ = -1 (K → ∞) to Δ = 1 (K = 1/2)
    ks = [
        QAtlas.fetch(XXZ1D(; J=1, Δ=Δ), LuttingerParameter(), Infinite()) for
        Δ in -0.9:0.2:0.9
    ]
    @test all(diff(ks) .< 0)
end

@testset "XXZ1D — NMRRelaxationExponent" begin
    # Δ = 0 (K = 1.0, XX free fermion) ⇒ θ_NMR = 1/(2K) - 1 = -1/2
    # (the textbook T^{-1/2} divergence of 1/T_1 in the XX chain)
    @test QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.0), NMRRelaxationExponent(), Infinite()) ≈ -0.5 atol =
        1e-12
    # Δ = 1 (K = 0.5, Heisenberg) ⇒ θ_NMR = 1/(2K) - 1 = 0.0 (constant, up to logs)
    @test QAtlas.fetch(XXZ1D(; J=1.0, Δ=1.0), NMRRelaxationExponent(), Infinite()) ≈ 0.0 atol =
        1e-12

    # Outside critical regime: returns NaN + warn.  The exponent delegates to
    # LuttingerParameter, which emits its own out-of-range :warn first, so use
    # match_mode=:any to assert our warning appears among the logs.
    @test isnan(
        @test_logs (:warn, r"critical Luttinger liquid regime") match_mode = :any QAtlas.fetch(
            XXZ1D(; J=1.0, Δ=1.5), NMRRelaxationExponent(), Infinite()
        )
    )
end

@testset "XXZ1D — LuttingerVelocity (+ SpinWaveVelocity alias)" begin
    # Canonical values at Δ=0 (XX) and Δ=1 (AF Heisenberg)
    @test QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.0), LuttingerVelocity(), Infinite()) ≈ 1.0 atol =
        1e-12
    @test QAtlas.fetch(XXZ1D(; J=1.0, Δ=1.0), LuttingerVelocity(), Infinite()) ≈ π / 2 atol =
        1e-12

    # J scaling linear
    @test QAtlas.fetch(XXZ1D(; J=3.0, Δ=0.0), LuttingerVelocity(), Infinite()) ≈ 3.0

    # SpinWaveVelocity() === LuttingerVelocity at the type level
    @test SpinWaveVelocity === LuttingerVelocity
    @test SpinWaveVelocity() isa LuttingerVelocity

    u1 = QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.5), LuttingerVelocity(), Infinite())
    u2 = QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.5), SpinWaveVelocity(), Infinite())
    @test u1 === u2
end

@testset "XXZ1D — central charge (critical regime only)" begin
    for Δ in (-0.9, -0.5, 0.0, 0.5, 0.99)
        @test QAtlas.fetch(XXZ1D(; J=1.0, Δ=Δ), CentralCharge(), Infinite()) == 1.0
    end
    # Outside the critical regime we return NaN + a warning.
    @test isnan(
        @test_logs (:warn, r"critical regime") QAtlas.fetch(
            XXZ1D(; J=1.0, Δ=1.5), CentralCharge(), Infinite()
        )
    )
end

@testset "XXZ1D (Δ=1) ↔ Heisenberg1D — dictionary cross-consistency" begin
    # The AF Heisenberg chain is the Δ = 1 point of the XXZ chain. QAtlas
    # exposes the ground-state energy density at this point through two
    # independent code paths — the dedicated `Heisenberg1D` model (Hulthén
    # closed form) and the generic `XXZ1D` model (Bethe-ansatz branch at
    # Δ = 1). The two must return identical values for every J.
    #
    # This catches any future divergence between the two dispatch paths
    # (e.g. a formula typo, a sign convention drift, an accidental J²
    # factor). Without this test the two are independent dictionary
    # entries that happen to compute the same physics.
    for J in (1.0, 2.5, -0.7)
        e_heis = QAtlas.fetch(Heisenberg1D(), GroundStateEnergyDensity(); J=J)
        e_xxz_energy = QAtlas.fetch(XXZ1D(; J=J, Δ=1.0), Energy(), Infinite())
        e_xxz_gs = QAtlas.fetch(XXZ1D(; J=J, Δ=1.0), GroundStateEnergyDensity(), Infinite())
        @test e_xxz_energy ≈ e_heis atol = 1e-12
        @test e_xxz_gs ≈ e_heis atol = 1e-12
    end

    # Sanity: Heisenberg1D value is exactly J·(1/4 − ln 2), so the XXZ
    # Δ=1 path also reproduces the literature closed form.
    @test QAtlas.fetch(XXZ1D(; J=1.0, Δ=1.0), Energy(), Infinite()) ≈ 0.25 - log(2.0) atol =
        1e-12
end

@testset "XXZ1D — gapped regime |Δ| > 1 is NaN (Orbach series deferred)" begin
    # |Δ| > 1 (Néel-like AF for Δ > 1, Ising-like FM for Δ < -1) is the
    # gapped regime; the Bethe ansatz takes a different series form
    # (Orbach 1958 / Walker 1959 / Yang-Yang III 1966) which is not yet
    # implemented.  The closed-form path warns and returns NaN there;
    # callers can use OBC dense ED for a finite-N reference.
    for Δ in (-2.0, -1.5, 1.5, 2.0)
        e = @test_logs (:warn, r"gapped regime") QAtlas.fetch(
            XXZ1D(; J=1.0, Δ=Δ), Energy(), Infinite()
        )
        @test isnan(e)
    end
end

@testset "XXZ1D — legacy Symbol dispatch routes to new API" begin
    # This testset deliberately reaches through the legacy Symbol shim in
    # `src/deprecate/legacy_xxz.jl` to confirm that the forwarding to the
    # concrete-struct fetch methods still produces the canonical values.
    # The shim emits a one-shot `@info` per (model, quantity) pair the
    # first time it's hit; `@test_logs` captures those infos so the
    # deprecation noise does not leak into CI output.
    e_sym = @test_logs (:info, r"symbol-dispatch") QAtlas.fetch(
        :XXZ, :energy, Infinite(); J=1.0, Δ=0.0
    )
    @test e_sym ≈ -1 / π atol = 1e-10

    u_sym = @test_logs (:info, r"symbol-dispatch") QAtlas.fetch(
        :XXZ, :spin_wave_velocity, Infinite(); J=1.0, Δ=1.0
    )
    @test u_sym ≈ π / 2 atol = 1e-12

    u_fv = @test_logs (:info, r"symbol-dispatch") QAtlas.fetch(
        :XXZ, :fermi_velocity, Infinite(); J=1.0, Δ=0.0
    )
    @test u_fv ≈ 1.0 atol = 1e-12

    K_sym = @test_logs (:info, r"symbol-dispatch") QAtlas.fetch(
        :XXZ, :luttinger_parameter, Infinite(); J=1.0, Δ=0.0
    )
    @test K_sym ≈ 1.0 atol = 1e-12
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "XXZ1D — verification cards" begin
    Sx, Sy, Sz = spin_ops(1 // 2)

    function xxz_bond(J, Δ)
        return J * (kron(Sx, Sx) + kron(Sy, Sy) + Δ * kron(Sz, Sz))
    end

    function xxz_e0_ed(J, Δ, N)
        return dense_spectrum(chain_hamiltonian(2, N, xxz_bond(J, Δ)))[1] / (N - 1)
    end

    Ns = verify_profile_Ns(; fast=(6, 8), full=(6, 8, 10, 12), nightly=(6, 8, 10, 12, 14))

    verify(
        XXZ1D(; J=1.0, Δ=0.0),
        Energy(),
        Infinite();
        route=:ed_finite_size,
        independent=[xxz_e0_ed(1.0, 0.0, N) for N in Ns],
        at=["N=$N" for N in Ns],
        agree_within=0.05,
        refs=["Yang-Yang 1966 I: e0 = -J/pi for Delta=0 (free fermion)"],
    )

    verify(
        XXZ1D(; J=1.0, Δ=1.0),
        Energy(),
        Infinite();
        route=:ed_finite_size,
        independent=[xxz_e0_ed(1.0, 1.0, N) for N in Ns],
        at=["N=$N" for N in Ns],
        agree_within=0.05,
        refs=["Hulthen 1938: e0 = J(1/4 - log 2) at Delta=1"],
    )

    verify(
        XXZ1D(; J=1.0, Δ=-1.0),
        Energy(),
        Infinite();
        route=:limiting_case,
        independent=-0.25,
        agree_within=1e-14,
        refs=["FM saturation: all-aligned state is exact GS, e0 = -J/4"],
    )

    verify(
        XXZ1D(; J=1.0, Δ=0.5),
        Energy(),
        Infinite();
        route=:ed_finite_size,
        independent=[xxz_e0_ed(1.0, 0.5, N) for N in Ns],
        at=["N=$N" for N in Ns],
        agree_within=0.05,
        refs=["Yang-Yang 1966 II: e0 = -3J/8 at Delta=1/2 (gamma=pi/3)"],
    )

    verify(
        XXZ1D(; J=1.0, Δ=0.0),
        LuttingerParameter(),
        Infinite();
        route=:second_closed_form,
        independent=1.0,
        agree_within=1e-12,
        refs=["Jordan-Wigner free fermion: K=1 at Delta=0"],
    )

    verify(
        XXZ1D(; J=1.0, Δ=1.0),
        LuttingerParameter(),
        Infinite();
        route=:limiting_case,
        independent=0.5,
        agree_within=1e-12,
        refs=["Luther-Peschel 1975: K=1/2 at the SU(2) isotropic point"],
    )

    verify(
        XXZ1D(; J=1.0, Δ=0.0),
        LuttingerVelocity(),
        Infinite();
        route=:second_closed_form,
        independent=1.0,
        agree_within=1e-12,
        refs=["Free fermion eps(k)=J cos k: v_F=J at k_F=pi/2"],
    )

    verify(
        XXZ1D(; J=1.0, Δ=1.0),
        LuttingerVelocity(),
        Infinite();
        route=:limiting_case,
        independent=π / 2,
        agree_within=1e-12,
        refs=["des Cloizeaux-Pearson 1962: u=piJ/2 at SU(2) isotropic point"],
    )

    verify(
        XXZ1D(; J=1.0, Δ=1.0),
        GroundStateEnergyDensity(),
        Infinite();
        route=:delegation_invariant,
        independent=QAtlas.fetch(Heisenberg1D(), GroundStateEnergyDensity(), Infinite()),
        agree_within=1e-12,
        refs=["XXZ1D at Delta=1 === Heisenberg1D: two independent code paths must agree"],
    )

    verify(
        XXZ1D(; J=1.0, Δ=0.0),
        CentralCharge(),
        Infinite();
        route=:second_closed_form,
        independent=1.0,
        agree_within=1e-14,
        refs=["Luttinger liquid: c=1 free compact boson CFT for |Delta| < 1"],
    )

    verify(
        XXZ1D(; J=1.0, Δ=0.0),
        NMRRelaxationExponent(),
        Infinite();
        route=:second_closed_form,
        independent=-0.5,
        agree_within=1e-12,
        refs=[
            "Sachdev 1994 (doi:10.1103/PhysRevB.50.13006) / Chitra & Giamarchi 1997 (doi:10.1103/PhysRevB.55.5816) Eq.27: leading 1/T_1 ∝ T^{1/(2K)-1}; at the XX free-fermion point K=1 ⇒ θ_NMR = -1/2 (T^{-1/2} divergence).",
        ],
    )
end
