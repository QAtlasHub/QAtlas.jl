# ─────────────────────────────────────────────────────────────────────────────
# Test: Majumdar–Ghosh chain — ED-finite-size verify for the Infinite BC.
#
# Split out of test_majumdar_ghosh.jl (10.4 min on s01). The previous file
# ran TWO ED sweeps (Infinite + PBC) in one shard for ~10 min total; this
# file isolates the Infinite sweep so it parallelises with the PBC sweep
# in test_majumdar_ghosh_verify_ed_pbc.jl. Each sweep is dominated by the
# dense-eigvals call at N=12 (4096-dim Hermitian).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "MajumdarGhosh — verify (Infinite ED N=6,8,10,12)" begin
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
end
