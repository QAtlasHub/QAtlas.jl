# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_thermo_PBC_ED_batch.jl
#
# ED-independent verify card for TFIM/{Energy,FreeEnergy,ThermalEntropy,
# SpecificHeat}/PBC.
#
# NOTE (2026-05-20, ED-verify-first policy): preliminary probe shows src
# disagrees with ED for PBC in the disordered phase h > J — same bug
# class as TFIM/SusceptibilityXX/PBC (#431): likely a Jordan-Wigner
# fermion-parity sector handling issue in the BdG PBC path. Ordered
# phase h < J matches ED.
#
# Bug-surfacing card per ED-verify-first policy. Refs #381, #431.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra

include(joinpath(@__DIR__, "..", "..", "..", "util", "generic_ed.jl"))

let Sx = spin_ops(1//2)[1], Sz = spin_ops(1//2)[3]
    sigmax = 2 * Sx
    sigmaz = 2 * Sz

    function ed_tfim_thermo_pbc(N::Int, J::Real, h::Real, beta::Real)
        H = Matrix(chain_hamiltonian_pbc(2, N, [(-J * sigmaz, sigmaz)]))
        for i in 1:N
            H .+= -h * site_op(sigmax, 2, N, i)
        end
        evals = sort(real.(eigvals(Hermitian(H))))
        E, F, S, C = thermo_from_spectrum(evals, beta)
        return E/N, F/N, S/N, C/N
    end

    @testset "TFIM — Energy + thermo/PBC vs ED (bug-surfacing at h>J) (#381 batch)" begin
        for (J, h) in ((1.0, 0.5), (1.0, 2.0))  # ordered + disordered phases
            for N in (4, 6, 8)
                for beta in (1.0, 5.0)
                    ed_E, ed_F, ed_S, ed_C = ed_tfim_thermo_pbc(N, J, h, beta)
                    for (qty, ed_val, lab) in (
                        (Energy(:per_site), ed_E, "E"),
                        (FreeEnergy(), ed_F, "F"),
                        (ThermalEntropy(), ed_S, "S"),
                        (SpecificHeat(), ed_C, "C"),
                    )
                        verify(
                            TFIM(; J=J, h=h),
                            qty,
                            PBC(N);
                            route=:ed_finite_size,
                            independent=ed_val,
                            at=["N=$(N)"],
                            agree_within=1e-9,
                            expected_fail=true,  # tracker issue #444 — TFIM PBC parity-sector handling (same root cause as #444)
                            refs=[
                                "ED black-box (PBC ring): chain_hamiltonian_pbc + onsite -h σ_x, full spectrum, thermo_from_spectrum ($(lab))",
                            ],
                            fetch_kw=(; beta=beta),
                        )
                    end
                end
            end
        end
    end
end
