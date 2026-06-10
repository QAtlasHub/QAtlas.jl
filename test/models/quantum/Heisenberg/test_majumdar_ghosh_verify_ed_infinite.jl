# ─────────────────────────────────────────────────────────────────────────────
# Test: Majumdar–Ghosh chain — ED-finite-size verify for the Infinite BC.
#
# Split out of test_majumdar_ghosh.jl (10.4 min on s01). The previous file
# ran TWO ED sweeps (Infinite + PBC) in one shard for ~10 min total; this
# file isolates the Infinite sweep so it parallelises with the PBC sweep
# in test_majumdar_ghosh_verify_ed_pbc.jl. Each sweep is dominated by the
# dense-eigvals call at N=12 (4096-dim Hermitian).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, KrylovKit, Random

@testset "MajumdarGhosh — verify (Infinite ED N=6,8,10,12)" begin
    # Independent J1-J2 PBC ground-state energy density (spin-1/2),
    # built black-box from site operators (J2 locked to J/2 by MG).
    # Sparse (Lanczos) GS: only the lowest eigenpair is needed, so we
    # avoid the O(D³) `eigvals(Hermitian(H))` cost — at N=12 that drops
    # the run from ~4 min (full 4096-dim eigvals) to ~1 s.
    function mg_pbc_e0(N, J)
        Sx, Sy, Sz = spin_ops(1 // 2)
        SS(i, j) =
            embed_two_site_sparse(Sx, Sx, i, j, N) +
            embed_two_site_sparse(Sy, Sy, i, j, N) +
            embed_two_site_sparse(Sz, Sz, i, j, N)
        H = sum(J * SS(i, mod1(i + 1, N)) + (J / 2) * SS(i, mod1(i + 2, N)) for i in 1:N)
        D = 2^N
        vals, _, _ = eigsolve(
            H,
            randn(MersenneTwister(0), ComplexF64, D),
            1,
            :SR;
            ishermitian=true,
            krylovdim=30,
            tol=1e-12,
        )
        return real(vals[1]) / N
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
