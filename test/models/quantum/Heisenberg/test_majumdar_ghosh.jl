# ─────────────────────────────────────────────────────────────────────────────
# Test: Majumdar–Ghosh chain (S = 1/2 J₁–J₂ at J₂/J₁ = 1/2)
#
# Values are verified by the verify() cards below (the new system). This
# file retains ONLY the structural / error / identity / relational guards
# that verify() architecturally cannot express (`@test_throws`,
# constructor invariants, registry sanity, deprecation-alias behaviour,
# SU(2) cross-checks, strict bounds). Legacy hand-rolled value @testsets
# are deleted — superseded by the verification cards.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test
using Logging: with_logger, NullLogger

@testset "MajumdarGhosh — structural / error / identity guards" begin
    @testset "Constructor variants agree" begin
        @test MajumdarGhosh(1.0).J == MajumdarGhosh(; J=1.0).J
    end

    @testset "GSED PBC — odd N rejected, N-kwarg form callable" begin
        m = MajumdarGhosh(; J=1.0)
        @test_throws DomainError QAtlas.fetch(m, GroundStateEnergyDensity(), PBC(5))
        # N-kwarg API form reachable (value covered by verify cards below).
        @test QAtlas.fetch(m, GroundStateEnergyDensity(), PBC(); N=8) isa Real
    end

    @testset "MassGap method dispatch & deprecation" begin
        m = MajumdarGhosh(; J=1.0)
        Δ_num = QAtlas.fetch(m, MassGap(), Infinite(); method=:numerical)
        Δ_trimer = QAtlas.fetch(m, MassGap(), Infinite(); method=:trimer_bound)
        # Default routes to :numerical (no fetch_kw == same path).
        @test QAtlas.fetch(m, MassGap(), Infinite()) == Δ_num
        # Strict relational: SS-1981 trimer-sector bound exceeds the
        # actual gap (cannot be expressed by a single-value verify card).
        @test Δ_trimer > Δ_num
        # Legacy :lower_bound alias resolves to J/4 with a one-shot @warn
        # (deprecation behaviour; not a verify-card concern).
        Δ_legacy = with_logger(NullLogger()) do
            QAtlas.fetch(m, MassGap(), Infinite(); method=:lower_bound)
        end
        @test Δ_legacy ≈ 0.25 atol = 1e-14
        # Unsupported method symbols raise DomainError.
        @test_throws DomainError QAtlas.fetch(m, MassGap(), Infinite(); method=:dmrg_strict)
        @test_throws DomainError QAtlas.fetch(m, MassGap(), Infinite(); method=:bogus)
    end

    @testset "Registry knows about MajumdarGhosh" begin
        rows = QAtlas.implementation_status(MajumdarGhosh)
        @test !isempty(rows)
        quantities = unique(r.quantity for r in rows)
        @test GroundStateEnergyDensity in quantities
        @test MassGap in quantities
    end

    @testset "SpinGap — strict bounds, DomainError, SU(2) identity" begin
        # Strict bounds (relational; not single-value verify cards):
        #   Δ < J/4  (Shastry-Sutherland 1981 trimer-sector bound)
        #   Δ > 0.117 J  (Magnus 1991 strict absolute-gap bound)
        Δ = QAtlas.fetch(MajumdarGhosh(), SpinGap(), Infinite())
        @test Δ < 0.25
        @test Δ > 0.117

        # SpinGap rejects J ≤ 0.
        @test_throws DomainError QAtlas.fetch(MajumdarGhosh(), SpinGap(), Infinite(); J=0.0)
        @test_throws DomainError QAtlas.fetch(
            MajumdarGhosh(), SpinGap(), Infinite(); J=-1.5
        )

        # SU(2) identity: MG is SU(2)-symmetric, so the spectral gap
        # equals the spin gap (S=0 → S=1 triplet excitation).
        for J in (0.5, 1.0, 3.0)
            m = MajumdarGhosh(; J=J)
            @test QAtlas.fetch(m, SpinGap(), Infinite()) ≈
                QAtlas.fetch(m, MassGap(), Infinite(); method=:numerical)
        end
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "MajumdarGhosh — verification cards" begin
    # Independent J1-J2 PBC ground-state energy density (spin-1/2),
    # built black-box from site operators (J2 locked to J/2 by MG).
    function mg_pbc_e0(N, J)
        Sx, Sy, Sz = spin_ops(1 // 2)
        SS(i, j) =
            site_op(Sx, 2, N, i) * site_op(Sx, 2, N, j) +
            site_op(Sy, 2, N, i) * site_op(Sy, 2, N, j) +
            site_op(Sz, 2, N, i) * site_op(Sz, 2, N, j)
        H = zeros(ComplexF64, 2^N, 2^N)
        for i in 1:N
            H .+= J * SS(i, mod1(i + 1, N))
            H .+= (J / 2) * SS(i, mod1(i + 2, N))
        end
        return dense_spectrum(H)[1] / N
    end

    # GSED Infinite: closed-form -3J/8 (each NN singlet ⟨S·S⟩ = -3/4; the
    # orthogonal-dimer covering gives e0 = -3J/8). J-scaling linear.
    for J in (0.5, 1.0, 2.0, 3.7)
        verify(
            MajumdarGhosh(; J=J),
            GroundStateEnergyDensity(),
            Infinite();
            route=:second_closed_form,
            independent=-3 * J / 8,
            agree_within=1e-14,
            refs=["Majumdar-Ghosh 1969: exact orthogonal-dimer product state, e0 = -3J/8"],
        )
    end

    # GSED Infinite ed_finite_size: independent J1-J2 ring ED via mg_pbc_e0,
    # exact ∀ even N (dimer state is the exact GS); non-circular vs the
    # closed form above.
    let Ns = verify_profile_Ns(; fast=(6, 8), full=(6, 8, 10, 12), nightly=(6, 8, 10, 12))
        verify(
            MajumdarGhosh(; J=1.0),
            GroundStateEnergyDensity(),
            Infinite();
            route=:ed_finite_size,
            independent=[mg_pbc_e0(N, 1.0) for N in Ns],
            at=["N=$N" for N in Ns],
            agree_within=1e-6,
            refs=["Exact MG dimer GS of the J1-J2 ring at J2=J/2 (even N), -3J/8"],
        )
    end

    # GSED PBC ed_finite_size: same independent ED with bc=PBC, so the
    # MajumdarGhosh/GroundStateEnergyDensity/PBC hub is corroborated.
    let Ns = verify_profile_Ns(; fast=(6, 8), full=(6, 8, 10, 12), nightly=(6, 8, 10, 12))
        verify(
            MajumdarGhosh(; J=1.0),
            GroundStateEnergyDensity(),
            PBC(8);
            route=:ed_finite_size,
            independent=[mg_pbc_e0(N, 1.0) for N in Ns],
            at=["N=$N" for N in Ns],
            # Dimer product state is the exact GS of the MG ring (J2=J/2)
            # for every even N, so the size-by-size residual is at machine
            # precision; 1e-12 is tight enough to catch real bugs but stays
            # well clear of accumulated float round-off in the ED solver.
            agree_within=1e-12,
            refs=[
                "Exact MG dimer GS of the J1-J2 ring at J2=J/2 " *
                "(even N), e0 = -3J/8 size-independent",
            ],
        )
    end

    # MassGap Infinite (:numerical default = White-Affleck DMRG ≈ 0.234 J).
    # J-scaling linear; src returns the literal 0.234*J to machine precision.
    # Replaces legacy J-loop value assertions.
    for J in (0.5, 1.0, 2.0, 3.7)
        verify(
            MajumdarGhosh(; J=J),
            MassGap(),
            Infinite();
            route=:literature_value,
            independent=0.234 * J,
            agree_within=1e-14,
            refs=["White-Affleck 1996 DMRG; Eggert 1996: Δ ≈ 0.234 J (J-linear)"],
        )
    end

    # MassGap Infinite (:trimer_bound = J/4, Shastry-Sutherland 1981).
    # Replaces legacy :trimer_bound value + J-scale tests.
    for J in (0.5, 1.0, 2.0, 3.7)
        verify(
            MajumdarGhosh(; J=J),
            MassGap(),
            Infinite();
            route=:literature_value,
            independent=J / 4,
            agree_within=1e-14,
            refs=["Shastry-Sutherland 1981: trimer-sector bound Δ ≥ J/4"],
            fetch_kw=(; method=:trimer_bound),
        )
    end

    # SpinGap Infinite (White-Affleck/Eggert DMRG ≈ 0.234 J), J-linear.
    for J in (0.5, 1.0, 3.0)
        verify(
            MajumdarGhosh(; J=J),
            SpinGap(),
            Infinite();
            route=:literature_value,
            independent=0.234 * J,
            agree_within=1e-14,
            refs=["White-Affleck 1996 DMRG; Eggert 1996: spin gap Δ ≈ 0.234 J"],
        )
    end
end
# ── additional verification cards (#381 batch 6) ─────────────────────────
@testset "MajumdarGhosh — SpinGap White-Affleck (#381 batch 6)" begin
    # Spin-1/2 J1-J2 chain at MG point J2 = J1/2: White-Affleck 1996 DMRG
    # singlet-triplet spin gap Δ_S ≈ 0.234 J.
    # Single J point: solver returns J * 0.234 from a stored constant;
    # a multi-J sweep collapses to identical constant-multiplication
    # residuals on the same code path and adds no discriminating power.
    let J = 1.0
        verify(
            MajumdarGhosh(; J=J),
            SpinGap(),
            Infinite();
            route=:literature_value,
            independent=0.234 * J,
            agree_within=5e-3,
            refs=["White-Affleck 1996 PRB 54 9862: MG point singlet-triplet spin gap Δ_S ≈ 0.234 J (DMRG)"],
        )
    end
end

