# ─────────────────────────────────────────────────────────────────────────────
# Test: Majumdar–Ghosh chain — literature / closed-form verify cards.
#
# Split out of test_majumdar_ghosh.jl (10.4 min on s01). This file contains
# ONLY the constant-time verify cards (closed-form GSED, literature MassGap
# default + trimer-bound, SpinGap White-Affleck) plus the #381 PBC closed-
# form batch and the #381 batch-6 SpinGap card. No ED — those live in the
# sibling test_majumdar_ghosh_verify_ed_{infinite,pbc}.jl.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "MajumdarGhosh — verify cards (literature / closed-form)" begin
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

    # MassGap Infinite (:numerical default = White-Affleck DMRG ≈ 0.234 J).
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

# ── additional verification card (#381 batch) ─────────────────────────────
@testset "MajumdarGhosh — GroundStateEnergyDensity/PBC closed-form (#381 batch)" begin
    # The dimer-product state is an exact eigenstate of the J1-J2 ring at
    # J2 = J/2 for any even N (Majumdar-Ghosh 1969), so e0 = -3J/8 is
    # size-independent on PBC.
    for N in (6, 8, 10, 12)
        for J in (0.5, 1.0, 2.0)
            verify(
                MajumdarGhosh(; J=J),
                GroundStateEnergyDensity(),
                PBC(N);
                route=:second_closed_form,
                independent=-3 * J / 8,
                agree_within=1e-14,
                refs=[
                    "Majumdar-Ghosh 1969: dimer GS exact for any even N, e0 = -3J/8 (BC- and size-independent)",
                ],
            )
        end
    end
end

# ── additional verification cards (#381 batch 6) ─────────────────────────
@testset "MajumdarGhosh — SpinGap White-Affleck (#381 batch 6)" begin
    # Spin-1/2 J1-J2 chain at MG point J2 = J1/2: White-Affleck 1996 DMRG
    # singlet-triplet spin gap Δ_S ≈ 0.234 J.
    let J = 1.0
        verify(
            MajumdarGhosh(; J=J),
            SpinGap(),
            Infinite();
            route=:literature_value,
            independent=0.234 * J,
            agree_within=5e-3,
            refs=[
                "White-Affleck 1996 PRB 54 9862: MG point singlet-triplet spin gap Δ_S ≈ 0.234 J (DMRG)",
            ],
        )
    end
end
