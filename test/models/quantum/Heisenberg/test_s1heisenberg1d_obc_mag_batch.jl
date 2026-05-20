# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_s1heisenberg1d_obc_mag_batch.jl
#
# SU(2)-symmetry verification cards for the OBC spin-1 Heisenberg chain:
# total magnetisations <S^α>/N (α ∈ {x,y,z}) vanish in the thermal
# ensemble at any J, N, β, by SU(2) symmetry of the Hamiltonian.
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "S1Heisenberg1D — Magnetization{X,Y,Z}/OBC = 0 (SU(2)) (#381 batch)" begin
    for J in (0.5, 1.0, 2.0)
        for N in (3, 4, 5, 6)
            for β in (1.0, 10.0)
                for q in (MagnetizationX(), MagnetizationY(), MagnetizationZ())
                    verify(
                        S1Heisenberg1D(),
                        q,
                        OBC(N);
                        route=:second_closed_form,
                        independent=0.0,
                        agree_within=1e-10,
                        refs=["SU(2) symmetry of S=1 Heisenberg: <S^α>_β = 0 for α ∈ {x,y,z} in the unbroken thermal ensemble"],
                        fetch_kw=(; J=J, beta=β),
                    )
                end
            end
        end
    end
end
