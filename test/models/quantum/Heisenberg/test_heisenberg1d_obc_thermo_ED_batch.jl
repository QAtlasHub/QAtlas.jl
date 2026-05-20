# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_heisenberg1d_obc_thermo_ED_batch.jl
#
# ED-independent structural corroboration for Heisenberg1D OBC
# thermodynamic observables at arbitrary β.  Builds H_Heisenberg from
# scratch via chain_hamiltonian, takes full spectrum, computes
# (E, F, S, C)/N via the standard log-sum-exp canonical formulas in
# thermo_from_spectrum.  No QAtlas thermal kernel touched.
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra

include(joinpath(@__DIR__, "..", "..", "..", "util", "generic_ed.jl"))

let Sx = spin_ops(1//2)[1], Sy = spin_ops(1//2)[2], Sz = spin_ops(1//2)[3]
    function ed_h1d_thermo_per_site(N::Int, J::Real, beta::Real)
        bond = J * (kron(Sx, Sx) + kron(Sy, Sy) + kron(Sz, Sz))
        H = chain_hamiltonian(2, N, bond)
        evals = dense_spectrum(H)
        E, F, S, C = thermo_from_spectrum(evals, beta)
        return E/N, F/N, S/N, C/N
    end

    @testset "Heisenberg1D — Energy + thermo/OBC vs ED at arbitrary β (#381 batch)" begin
        for J in (0.5, 1.0, 2.0)
            for N in (4, 6, 8)
                for beta in (0.5, 2.0, 10.0)
                    ed_E, ed_F, ed_S, ed_C = ed_h1d_thermo_per_site(N, J, beta)
                    verify(
                        Heisenberg1D(),
                        Energy(:per_site),
                        OBC(N);
                        route=:ed_finite_size,
                        independent=ed_E,
                        at=["N=$(N)", "β=$(beta)"],
                        agree_within=1e-9,
                        refs=[
                            "ED black-box: chain_hamiltonian(2,N, J·(Sx⊗Sx+Sy⊗Sy+Sz⊗Sz)), thermo_from_spectrum",
                        ],
                        fetch_kw=(; J=J, beta=beta),
                    )
                    verify(
                        Heisenberg1D(),
                        FreeEnergy(),
                        OBC(N);
                        route=:ed_finite_size,
                        independent=ed_F,
                        at=["N=$(N)", "β=$(beta)"],
                        agree_within=1e-9,
                        refs=["ED black-box: full-spectrum log-sum-exp free energy F/N"],
                        fetch_kw=(; J=J, beta=beta),
                    )
                    verify(
                        Heisenberg1D(),
                        ThermalEntropy(),
                        OBC(N);
                        route=:ed_finite_size,
                        independent=ed_S,
                        at=["N=$(N)", "β=$(beta)"],
                        agree_within=1e-9,
                        refs=["ED black-box: S = β·(E - F) from full spectrum"],
                        fetch_kw=(; J=J, beta=beta),
                    )
                    verify(
                        Heisenberg1D(),
                        SpecificHeat(),
                        OBC(N);
                        route=:ed_finite_size,
                        independent=ed_C,
                        at=["N=$(N)", "β=$(beta)"],
                        agree_within=1e-9,
                        refs=["ED black-box: C = β²·(⟨E²⟩ - ⟨E⟩²) from full spectrum"],
                        fetch_kw=(; J=J, beta=beta),
                    )
                end
            end
        end
    end
end
