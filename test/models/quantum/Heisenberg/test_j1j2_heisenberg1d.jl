# test/models/quantum/Heisenberg/test_j1j2_heisenberg1d.jl
#
# Phase-1 closed-form tests for J1J2Heisenberg1D (issue #297).  Only the
# two analytically-known points are exercised here:
#
#   - j = 0   → Bethe-Hulthén  (delegate Heisenberg1D)
#   - j = 1/2 → Majumdar-Ghosh (delegate MajumdarGhosh)
#
# Generic j (DMRG) is deferred to Phase 2 and must raise DomainError.

using QAtlas, Test

@testset "J1J2Heisenberg1D — j = 0 (pure Heisenberg, Bethe-Hulthén) (Phase 1)" begin
    e0 = QAtlas.fetch(J1J2Heisenberg1D(; J1=1.0, J2=0.0), Energy{:per_site}(), Infinite())
    @test e0 ≈ 0.25 - log(2)
    # Linear scaling with J₁
    e0_3 = QAtlas.fetch(J1J2Heisenberg1D(; J1=3.0, J2=0.0), Energy{:per_site}(), Infinite())
    @test e0_3 ≈ 3 * (0.25 - log(2))
    # Delegation matches Heisenberg1D directly (legacy GroundStateEnergyDensity API)
    @test e0 ≈ QAtlas.fetch(Heisenberg1D(), GroundStateEnergyDensity(), Infinite(); J=1.0)
end

@testset "J1J2Heisenberg1D — j = 1/2 (Majumdar-Ghosh) (Phase 1)" begin
    e0 = QAtlas.fetch(J1J2Heisenberg1D(; J1=1.0, J2=0.5), Energy{:per_site}(), Infinite())
    @test e0 ≈ -3 / 8
    # Default constructor IS the MG point
    @test QAtlas.fetch(J1J2Heisenberg1D(), Energy{:per_site}(), Infinite()) ≈ -3 / 8
    # Linear scaling — (J₁, J₂) = (2, 1) is still j = 1/2, so E/N = -3·2/8 = -3/4
    e0_2 = QAtlas.fetch(J1J2Heisenberg1D(; J1=2.0, J2=1.0), Energy{:per_site}(), Infinite())
    @test e0_2 ≈ -3 / 4
    # Delegation matches MajumdarGhosh directly (legacy GroundStateEnergyDensity API)
    @test e0 ≈ QAtlas.fetch(MajumdarGhosh(; J=1.0), GroundStateEnergyDensity(), Infinite())
end

@testset "J1J2Heisenberg1D — generic j throws DomainError (Phase 1)" begin
    # Various non-closed-form j values bracketing the MG point
    for (J1, J2) in [(1.0, 0.25), (1.0, 0.75), (2.0, 0.3), (1.5, 1.0)]
        @test_throws DomainError QAtlas.fetch(
            J1J2Heisenberg1D(; J1=J1, J2=J2), Energy{:per_site}(), Infinite()
        )
    end
end

@testset "J1J2Heisenberg1D — rejects J1 ≤ 0 or J2 < 0 (Phase 1)" begin
    @test_throws DomainError J1J2Heisenberg1D(; J1=0.0, J2=0.5)
    @test_throws DomainError J1J2Heisenberg1D(; J1=-1.0, J2=0.5)
    @test_throws DomainError J1J2Heisenberg1D(; J1=1.0, J2=-0.1)
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "J1J2Heisenberg1D — verification cards" begin
    # Independent J1-J2 PBC ground-state energy density (spin-1/2),
    # built black-box from site operators (never QAtlas internals).
    function j1j2_pbc_e0(N, J1, J2)
        Sx, Sy, Sz = spin_ops(1 // 2)
        SS(i, j) =
            site_op(Sx, 2, N, i) * site_op(Sx, 2, N, j) +
            site_op(Sy, 2, N, i) * site_op(Sy, 2, N, j) +
            site_op(Sz, 2, N, i) * site_op(Sz, 2, N, j)
        H = zeros(ComplexF64, 2^N, 2^N)
        for i in 1:N
            H .+= J1 * SS(i, mod1(i + 1, N))
            H .+= J2 * SS(i, mod1(i + 2, N))
        end
        return dense_spectrum(H)[1] / N
    end

    # j = 0 (pure Heisenberg): delegates to Heisenberg1D, e0 = J1(1/4 - log 2)
    verify(
        J1J2Heisenberg1D(; J1=1.0, J2=0.0),
        Energy(:per_site),
        Infinite();
        route=:delegation_invariant,
        independent=QAtlas.fetch(Heisenberg1D(), GroundStateEnergyDensity(), Infinite()),
        agree_within=1e-12,
        refs=["J1J2 at J2=0 delegates to Heisenberg1D (Hulthen 1938)"],
    )

    # j = 1/2 (Majumdar-Ghosh): delegates to MajumdarGhosh, e0 = -3 J1 / 8
    verify(
        J1J2Heisenberg1D(; J1=1.0, J2=0.5),
        Energy(:per_site),
        Infinite();
        route=:delegation_invariant,
        independent=QAtlas.fetch(MajumdarGhosh(; J=1.0), GroundStateEnergyDensity(), Infinite()),
        agree_within=1e-12,
        refs=["J1J2 at J2=J1/2 delegates to MajumdarGhosh (exact dimer, -3J/8)"],
    )

    # j = 1/2 MG point: independent PBC ED -> exact dimer -3 J1 / 8 (even N)
    let Ns = verify_profile_Ns(; fast=(6, 8), full=(6, 8, 10, 12), nightly=(6, 8, 10, 12))
        verify(
            J1J2Heisenberg1D(; J1=1.0, J2=0.5),
            Energy(:per_site),
            Infinite();
            route=:ed_finite_size,
            independent=[j1j2_pbc_e0(N, 1.0, 0.5) for N in Ns],
            at=["N=$N" for N in Ns],
            agree_within=1e-6,
            refs=["Majumdar-Ghosh 1969: exact dimer GS, E0/N = -3J/8 size-independent (PBC even N)"],
        )
    end
end
