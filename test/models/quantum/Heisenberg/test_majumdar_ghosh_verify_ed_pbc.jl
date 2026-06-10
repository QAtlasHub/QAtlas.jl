# ─────────────────────────────────────────────────────────────────────────────
# Test: Majumdar–Ghosh chain — ED-finite-size verify for PBC.
#
# Split out of test_majumdar_ghosh.jl (10.4 min on s01). Parallel to the
# sibling test_majumdar_ghosh_verify_ed_infinite.jl. Each sweep is
# dominated by the dense-eigvals call at N=12 (4096-dim Hermitian).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, KrylovKit, Random

@testset "MajumdarGhosh — verify (PBC ED N=6,8,10,12)" begin
    # Sparse Lanczos GS: build H as a SparseMatrixCSC (embed_two_site_sparse,
    # O(2^N) nnz per bond) instead of dense site_op products, then take only
    # the lowest eigenpair. This is what the O(D³)→O(nnz·k) note always meant;
    # the dense build was the real cost (N=12: ~8 min → ~1 s). Eigenvalues are
    # basis-independent, so the exact MG dimer energy -3J/8 is reproduced.
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
