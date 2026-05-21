# ─────────────────────────────────────────────────────────────────────────────
# AKLT1D — closed-form / DMRG cards (#381 batch) (split, formerly testset 3 of 3)
#
# Split out of test/models/quantum/misc/test_aklt.jl (5.9 min on s02) so
# the three top-level testsets each run on their own shard. Helpers
# spin_ops, chain_hamiltonian, two_point, verify_profile_Ns come from
# test/util/{generic_ed,verify}.jl via runtests.jl ambient include.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra

@testset "AKLT1D — closed-form / DMRG cards (#381 batch)" begin
    # GroundStateEnergyDensity/Infinite: AKLT 1988 exact VBS GS energy
    # density e₀ = -2J/3, J-linear, J-independent of the wavefunction.
    for J in (0.5, 1.0, 2.0, 3.7)
        verify(
            AKLT1D(; J=J),
            GroundStateEnergyDensity(),
            Infinite();
            route=:second_closed_form,
            independent=-2 * J / 3,
            agree_within=1e-14,
            refs=[
                "AKLT 1988: VBS ground state is the exact null space of every bond P₂ projector ⇒ e₀ = -2J/3",
            ],
        )
    end

    # CorrelationLength/Infinite: ξ = 1/log 3 (AKLT 1988); J-independent
    # because the VBS wavefunction does not depend on J > 0.
    for J in (0.5, 1.0, 2.0)
        verify(
            AKLT1D(; J=J),
            CorrelationLength(),
            Infinite();
            route=:second_closed_form,
            independent=1 / log(3),
            agree_within=1e-14,
            refs=["AKLT 1988: ⟨S^z_0 S^z_r⟩ = (-1)^r (4/3) 3^{-r} ⇒ ξ = 1/log 3"],
        )
    end

    # StringOrderParameter/Infinite: O_str = 4/9 (AKLT 1988, Kennedy-Tasaki
    # 1992 hidden Z2×Z2 symmetry-breaking order parameter); J-independent.
    for J in (0.5, 1.0, 2.0)
        verify(
            AKLT1D(; J=J),
            StringOrderParameter(),
            Infinite();
            route=:second_closed_form,
            independent=4 / 9,
            agree_within=1e-14,
            refs=[
                "Kennedy-Tasaki 1992 on AKLT VBS: O_str = -⟨S^z_i e^{iπ Σ S^z_k} S^z_j⟩ → 4/9 at r → ∞",
            ],
        )
    end

    # MassGap/Infinite: Haldane gap Δ ≈ 0.350 J — DMRG literature value
    # (García-Saez/Murg/Verstraete 2013, PRB 88, 245118); no closed form.
    #
    # IMPORTANT — agree_within=5e-3 is the DMRG literature uncertainty floor
    # (García-Saez–Murg–Verstraete report Δ/J ≈ 0.3502 ± 0.0001; ~25× their
    # quoted uncertainty gives headroom for hub-stored precision choices
    # between 0.350 and 0.35048). This tolerance is INTENTIONALLY loose and
    # should NOT be tightened by future maintainers without a new high-
    # precision DMRG/MERA result superseding GMV 2013.
    for J in (1.0, 2.0)
        verify(
            AKLT1D(; J=J),
            MassGap(),
            Infinite();
            route=:literature_value,
            independent=0.350 * J,
            agree_within=5e-3,
            refs=[
                "García-Saez–Murg–Verstraete 2013 (PRB 88, 245118): AKLT Haldane gap Δ ≈ 0.350 J (DMRG)",
            ],
        )
    end
end
