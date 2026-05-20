# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_s1heisenberg1d_obc_thermo_ED_batch.jl
#
# ED-independent structural corroboration for S1Heisenberg1D OBC
# thermodynamic observables at arbitrary β. Builds H from scratch with
# spin-1 primitives, takes full spectrum (N≤5; 3^5=243 is still fast for full
# spectrum), computes per-site (E, F, S, C). Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra

include(joinpath(@__DIR__, "..", "..", "..", "util", "generic_ed.jl"))

let Sx = spin_ops(1)[1], Sy = spin_ops(1)[2], Sz = spin_ops(1)[3]
    function ed_s1h_thermo_per_site(N::Int, J::Real, beta::Real)
        bond = J * (kron(Sx, Sx) + kron(Sy, Sy) + kron(Sz, Sz))
        H = chain_hamiltonian(3, N, bond)
        evals = dense_spectrum(H)
        E, F, S, C = thermo_from_spectrum(evals, beta)
        return E/N, F/N, S/N, C/N
    end

    @testset "S1Heisenberg1D — Energy + thermo/OBC vs ED at arbitrary β (#381 batch)" begin
        for J in (0.5, 1.0, 2.0)
            for N in (3, 4, 5)  # 3^5 = 243; still fast for full spectrum
                for beta in (0.5, 2.0, 5.0)
                    ed_E, ed_F, ed_S, ed_C = ed_s1h_thermo_per_site(N, J, beta)
                    verify(
                        S1Heisenberg1D(),
                        Energy(:per_site),
                        OBC(N);
                        route=:ed_finite_size,
                        independent=ed_E,
                        at=["N=$(N)"],
                        agree_within=1e-9,
                        refs=[
                            "ED black-box (spin-1): chain_hamiltonian(3,N, J·(SxSx+SySy+SzSz)), thermo_from_spectrum",
                        ],
                        fetch_kw=(; J=J, beta=beta),
                    )
                    verify(
                        S1Heisenberg1D(),
                        FreeEnergy(),
                        OBC(N);
                        route=:ed_finite_size,
                        independent=ed_F,
                        at=["N=$(N)"],
                        agree_within=1e-9,
                        refs=["ED black-box (spin-1): full-spectrum log-sum-exp F/N"],
                        fetch_kw=(; J=J, beta=beta),
                    )
                    verify(
                        S1Heisenberg1D(),
                        ThermalEntropy(),
                        OBC(N);
                        route=:ed_finite_size,
                        independent=ed_S,
                        at=["N=$(N)"],
                        agree_within=1e-9,
                        refs=["ED black-box (spin-1): S = β·(E - F) from full spectrum"],
                        fetch_kw=(; J=J, beta=beta),
                    )
                    verify(
                        S1Heisenberg1D(),
                        SpecificHeat(),
                        OBC(N);
                        route=:ed_finite_size,
                        independent=ed_C,
                        at=["N=$(N)"],
                        agree_within=1e-9,
                        refs=["ED black-box (spin-1): C = β²·Var(E) from full spectrum"],
                        fetch_kw=(; J=J, beta=beta),
                    )
                end
            end
        end
    end
end
