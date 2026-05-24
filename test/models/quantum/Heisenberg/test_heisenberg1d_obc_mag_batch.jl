# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_heisenberg1d_obc_mag_batch.jl
#
# SU(2)-symmetry verification cards for the OBC Heisenberg chain: the
# total magnetisations <S^α>/N (α ∈ {x, y, z}) are identically zero in
# the thermal ensemble at any temperature and any J, N, because the
# Hamiltonian commutes with each S^α_total generator, so the thermal
# trace over the symmetric basis annihilates linear magnetisation.
#
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Heisenberg1D — Magnetization{X,Y,Z}/OBC = 0 (SU(2)) (#381 batch)" begin
    for J in (0.5, 1.0, 2.0)
        for N in (4, 6, 8)
            for β in (1.0, 10.0)
                verify(
                    Heisenberg1D(),
                    MagnetizationX(),
                    OBC(N);
                    route=:limiting_case,
                    independent=0.0,
                    agree_within=1e-10,
                    refs=[
                        "SU(2) symmetry of the Heisenberg Hamiltonian: <S^α>_β = 0 for α ∈ {x,y,z} in the unbroken thermal ensemble, exact for any J, N, β",
                    ],
                    fetch_kw=(; J=J, beta=β),
                )
                verify(
                    Heisenberg1D(),
                    MagnetizationY(),
                    OBC(N);
                    route=:limiting_case,
                    independent=0.0,
                    agree_within=1e-10,
                    refs=[
                        "SU(2) symmetry of the Heisenberg Hamiltonian: <S^α>_β = 0 for α ∈ {x,y,z} in the unbroken thermal ensemble, exact for any J, N, β",
                    ],
                    fetch_kw=(; J=J, beta=β),
                )
                verify(
                    Heisenberg1D(),
                    MagnetizationZ(),
                    OBC(N);
                    route=:limiting_case,
                    independent=0.0,
                    agree_within=1e-10,
                    refs=[
                        "SU(2) symmetry of the Heisenberg Hamiltonian: <S^α>_β = 0 for α ∈ {x,y,z} in the unbroken thermal ensemble, exact for any J, N, β",
                    ],
                    fetch_kw=(; J=J, beta=β),
                )
            end
        end
    end
end
