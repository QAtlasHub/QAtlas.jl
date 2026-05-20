# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/XXZ/test_xxz1d_obc_thermo_ED_batch.jl
#
# ED-independent verify cards for XXZ1D OBC thermodynamic observables.
#
# NOTE (2026-05-20, ED-verify-first policy): preliminary probe shows src
# returns values that DIFFER from ED for ALL (Δ, N, β) tested — including
# Δ=1 isotropic where XXZ1D should agree with Heisenberg1D. Heisenberg1D
# src matches ED exactly, XXZ1D(Δ=1) src disagrees → same bug class as
# the Susceptibility one (#428): the XXZ thermal kernel does not
# correctly propagate Δ (or even build the right H at Δ=1).
# Bug-surfacing card per ED-verify-first policy. Refs #381, #428.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra

include(joinpath(@__DIR__, "..", "..", "..", "util", "generic_ed.jl"))

let Sx = spin_ops(1//2)[1],
    Sy = spin_ops(1//2)[2],
    Sz = spin_ops(1//2)[3]

    function ed_xxz_thermo_per_site(N::Int, J::Real, dz::Real, beta::Real)
        bond = J * (kron(Sx, Sx) + kron(Sy, Sy) + dz * kron(Sz, Sz))
        H = chain_hamiltonian(2, N, bond)
        evals = dense_spectrum(H)
        E, F, S, C = thermo_from_spectrum(evals, beta)
        return E/N, F/N, S/N, C/N
    end

    @testset "XXZ1D — Energy + thermo/OBC vs ED (bug-surfacing) (#381 batch)" begin
        for (J, dz) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0))
            for N in (3, 4)
                for beta in (0.5, 5.0)
                    ed_E, ed_F, ed_S, ed_C = ed_xxz_thermo_per_site(N, J, dz, beta)
                    verify(
                        XXZ1D(),
                        Energy(:per_site),
                        OBC(N);
                        route=:ed_finite_size,
                        independent=ed_E,
                        at=["N=$(N)"],
                        agree_within=1e-9,
                        refs=["ED black-box: chain_hamiltonian(2,N, J·(SxSx+SySy+Δ·SzSz)), thermo_from_spectrum"],
                        fetch_kw=(; J=J, Δ=dz, beta=beta),
                    )
                    verify(
                        XXZ1D(),
                        FreeEnergy(),
                        OBC(N);
                        route=:ed_finite_size,
                        independent=ed_F,
                        at=["N=$(N)"],
                        agree_within=1e-9,
                        refs=["ED black-box: full-spectrum log-sum-exp F/N"],
                        fetch_kw=(; J=J, Δ=dz, beta=beta),
                    )
                    verify(
                        XXZ1D(),
                        ThermalEntropy(),
                        OBC(N);
                        route=:ed_finite_size,
                        independent=ed_S,
                        at=["N=$(N)"],
                        agree_within=1e-9,
                        refs=["ED black-box: S = β·(E - F) from full spectrum"],
                        fetch_kw=(; J=J, Δ=dz, beta=beta),
                    )
                    verify(
                        XXZ1D(),
                        SpecificHeat(),
                        OBC(N);
                        route=:ed_finite_size,
                        independent=ed_C,
                        at=["N=$(N)"],
                        agree_within=1e-9,
                        refs=["ED black-box: C = β²·Var(E) from full spectrum"],
                        fetch_kw=(; J=J, Δ=dz, beta=beta),
                    )
                end
            end
        end
    end
end
