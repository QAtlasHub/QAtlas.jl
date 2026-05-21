# ─────────────────────────────────────────────────────────────────────────────
# Test: Majumdar–Ghosh chain — ED-finite-size verify for PBC.
#
# Split out of test_majumdar_ghosh.jl (10.4 min on s01). Parallel to the
# sibling test_majumdar_ghosh_verify_ed_infinite.jl. Each sweep is
# dominated by the dense-eigvals call at N=12 (4096-dim Hermitian).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "MajumdarGhosh — verify (PBC ED N=6,8,10,12)" begin
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

    let Ns = verify_profile_Ns(; fast=(6, 8), full=(6, 8, 10, 12), nightly=(6, 8, 10, 12))
        verify(
            MajumdarGhosh(; J=1.0),
            GroundStateEnergyDensity(),
            PBC(8);
            route=:ed_finite_size,
            independent=[mg_pbc_e0(N, 1.0) for N in Ns],
            at=["N=$N" for N in Ns],
            # Dimer product state is the exact GS of the MG ring (J2=J/2)
            # for every even N — machine-precision tolerance.
            agree_within=1e-12,
            refs=[
                "Exact MG dimer GS of the J1-J2 ring at J2=J/2 " *
                "(even N), e0 = -3J/8 size-independent",
            ],
        )
    end
end
