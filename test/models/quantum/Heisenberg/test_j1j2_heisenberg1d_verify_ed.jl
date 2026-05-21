# test/models/quantum/Heisenberg/test_j1j2_heisenberg1d_verify_ed.jl
#
# Split out of test_j1j2_heisenberg1d.jl (9.9 min on s02). Heavy ED sweep
# at the MG point (PBC, N=6,8,10,12 with dense Hermitian eigvals; dominated
# by N=12 = 4096-dim). Isolated so it parallelises in its own shard.

using QAtlas, Test, KrylovKit, Random

@testset "J1J2Heisenberg1D — verify (MG point PBC ED N=6,8,10,12)" begin
    # Independent J1-J2 PBC ground-state energy density (spin-1/2),
    # built black-box from site operators (never QAtlas internals).
    # Sparse Lanczos GS — only lowest eigenpair needed, drops N=12 from
    # ~4 min (full eigvals) to ~1 s.
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
            refs=[
                "Majumdar-Ghosh 1969: exact dimer GS, E0/N = -3J/8 size-independent (PBC even N)",
            ],
        )
    end
end
