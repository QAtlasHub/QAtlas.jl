# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: HeisenbergXYZ — axis-aligned XXZ-delegation.
#
# Verifies:
#   * Jx = Jy = J, Jz = J (isotropic AF Heisenberg):
#         E0/N = J (1/4 - log 2)   (Hulthén 1938 via XXZ1D Δ = 1 path)
#   * Jx = Jy = J, Jz = 0 (XX free fermion): E0/N = -J/π
#   * Jx = Jy = J, Jz = -J: E0/N = -J/4 (FM saturation)
#   * Jx = Jy = J, Jz = J/2: matches XXZ1D Yang-Yang Δ = 1/2 → -3J/8
#   * Jx ≠ Jy: DomainError (general XYZ deferred to Phase 2)
#   * Jx = Jy = 0: DomainError (Ising-like reduction)
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "HeisenbergXYZ — isotropic Heisenberg AF limit (Hulthén)" begin
    m = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=1.0)
    E0 = QAtlas.fetch(m, Energy(:per_site), Infinite())
    @test E0 ≈ 0.25 - log(2.0) atol = 1e-12
end

@testset "HeisenbergXYZ — XX free-fermion limit" begin
    m = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=0.0)
    E0 = QAtlas.fetch(m, Energy(:per_site), Infinite())
    @test E0 ≈ -1 / π atol = 1e-12
end

@testset "HeisenbergXYZ — isotropic FM saturated limit" begin
    m = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=-1.0)
    E0 = QAtlas.fetch(m, Energy(:per_site), Infinite())
    @test E0 ≈ -1 / 4 atol = 1e-12
end

@testset "HeisenbergXYZ — Yang-Yang single integral at Δ = 1/2" begin
    m = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=0.5)
    E0 = QAtlas.fetch(m, Energy(:per_site), Infinite())
    @test E0 ≈ -3 / 8 atol = 1e-10
end

@testset "HeisenbergXYZ — J-scaling delegated to XXZ1D" begin
    m1 = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=0.5)
    m2 = HeisenbergXYZ(; Jx=2.5, Jy=2.5, Jz=1.25)
    e1 = QAtlas.fetch(m1, Energy(:per_site), Infinite())
    e2 = QAtlas.fetch(m2, Energy(:per_site), Infinite())
    # XXZ1D returns J × (Δ-only function), so e2 = 2.5 e1.
    @test e2 ≈ 2.5 * e1 atol = 1e-10
end

@testset "HeisenbergXYZ — Jx ≠ Jy raises DomainError" begin
    @test_throws DomainError QAtlas.fetch(
        HeisenbergXYZ(; Jx=1.0, Jy=0.5, Jz=0.3), Energy(:per_site), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        HeisenbergXYZ(; Jx=1.0, Jy=2.0, Jz=1.0), Energy(:per_site), Infinite()
    )
end

@testset "HeisenbergXYZ — Jx = Jy = 0 raises DomainError" begin
    @test_throws DomainError QAtlas.fetch(
        HeisenbergXYZ(; Jx=0.0, Jy=0.0, Jz=1.0), Energy(:per_site), Infinite()
    )
end

@testset "HeisenbergXYZ — LuttingerParameter at isotropic point (Phase 2)" begin
    # Isotropic SU(2) point: K = 1/2 (Luther-Peschel 1975).
    # Strict ==: acos(1.0) == 0.0 in IEEE → π/(2π) == 0.5 exactly.
    m = HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=1.0)
    @test QAtlas.fetch(m, LuttingerParameter(), Infinite()) == 0.5
    # Delegation invariant: bit-identical to XXZ1D(Δ=1) directly.
    @test QAtlas.fetch(m, LuttingerParameter(), Infinite()) ===
        QAtlas.fetch(QAtlas.XXZ1D(; J=1.0, Δ=1.0), LuttingerParameter(), Infinite())
end

@testset "HeisenbergXYZ — LuttingerParameter non-isotropic throws DomainError (Phase 2 deferral)" begin
    # Jx = Jy ≠ Jz: axial anisotropy (use XXZ1D directly)
    @test_throws DomainError QAtlas.fetch(
        HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=0.5), LuttingerParameter(), Infinite()
    )
    # Generic XYZ: Baxter elliptic regime, deferred
    @test_throws DomainError QAtlas.fetch(
        HeisenbergXYZ(; Jx=1.0, Jy=0.7, Jz=0.5), LuttingerParameter(), Infinite()
    )
    # Regression: strict `Jx == Jy == Jz` — tiny deviations must NOT silently
    # delegate to XXZ1D(Δ=1). K is not constant in a neighbourhood of the
    # isotropic point (XXZ K(Δ) varies with Δ; XYZ requires Baxter elliptic).
    @test_throws DomainError QAtlas.fetch(
        HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=1.0 + 1e-13), LuttingerParameter(), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        HeisenbergXYZ(; Jx=1.0, Jy=1.0 + 1e-13, Jz=1.0), LuttingerParameter(), Infinite()
    )
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "HeisenbergXYZ — verification cards" begin
    Sx, Sy, Sz = spin_ops(1 // 2)
    Ns = verify_profile_Ns(; fast=(6, 8), full=(6, 8, 10, 12), nightly=(6, 8, 10, 12, 14))

    # Isotropic Heisenberg AF: E/site converges to J(1/4 - log 2)
    let J = 1.0
        bond = J * (kron(Sx, Sx) + kron(Sy, Sy) + kron(Sz, Sz))
        ind = [dense_spectrum(chain_hamiltonian(2, N, bond))[1] / (N - 1) for N in Ns]
        verify(
            HeisenbergXYZ(; Jx=J, Jy=J, Jz=J),
            Energy(:per_site),
            Infinite();
            route=:ed_finite_size,
            independent=ind,
            at=["N=$N" for N in Ns],
            agree_within=0.05,
            refs=["Hulthen 1938 via delegation to XXZ1D at Delta=1"],
        )
    end

    # XX limit (Jz=0): e0 = -J/pi (free fermion, Jordan-Wigner)
    verify(
        HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=0.0),
        Energy(:per_site),
        Infinite();
        route=:limiting_case,
        independent=-1 / π,
        agree_within=1e-12,
        refs=["Jordan-Wigner free fermion: XX limit Jz=0 gives e0 = -J/pi"],
    )

    # FM limit (Jz=-J): e0 = -J/4 (saturated ferromagnet, exact)
    verify(
        HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=-1.0),
        Energy(:per_site),
        Infinite();
        route=:limiting_case,
        independent=-1 / 4,
        agree_within=1e-12,
        refs=["FM saturation: Jz=-J gives e0 = -J/4"],
    )

    # LuttingerParameter delegation invariant at isotropic SU(2) point
    verify(
        HeisenbergXYZ(; Jx=1.0, Jy=1.0, Jz=1.0),
        LuttingerParameter(),
        Infinite();
        route=:delegation_invariant,
        independent=QAtlas.fetch(XXZ1D(; J=1.0, Δ=1.0), LuttingerParameter(), Infinite()),
        agree_within=1e-14,
        refs=["HeisenbergXYZ(isotropic) delegates to XXZ1D(Delta=1): K=1/2"],
    )
end
