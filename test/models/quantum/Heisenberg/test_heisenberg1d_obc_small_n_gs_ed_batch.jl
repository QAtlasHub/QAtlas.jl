# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_heisenberg1d_obc_small_n_gs_ed_batch.jl
#
# Small-N OBC ED pins for Heisenberg1D ground-state energy.
#
# Restores a piece of the WHY-plane coverage removed in PR #449
# (zero-legacy phase 1): the deleted
# test/verification/heisenberg_xxz/test_heisenberg_dimer.jl pinned
# the N=2 dimer ground state (singlet, E₀ = -3J/4) by ED against
# QAtlas's Heisenberg1D ExactSpectrum dispatch. INVENTORY had ZERO
# Heisenberg1D/Energy/OBC small-N entries before this PR.
#
# Independent route: build H_Heisenberg from scratch via
# chain_hamiltonian(2, N, J·(Sx⊗Sx + Sy⊗Sy + Sz⊗Sz)) → dense_spectrum,
# take the minimum eigenvalue.  Subject: thermal Energy{:total}/OBC at
# β = 1e6 (low-T limit collapses to the GS energy; small-N Heisenberg
# OBC is finite-gapped, so all excited-state Boltzmann weights underflow
# to exactly 0.0 in float64).  No QAtlas internal builder is invoked
# on the independent side ⇒ true black-box cross-check.
#
# (N=4 PBC was the second leg of #449's deleted dimer/4-site batch;
# Heisenberg1D Energy{:total}/PBC has no src dispatch in QAtlas v0.21,
# so the PBC complement is deferred to a separate PR that adds that
# dispatch first.)
#
# Pure verify(); branches off main. Refs #381; restores coverage lost in #449.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra

include(joinpath(@__DIR__, "..", "..", "..", "util", "generic_ed.jl"))

const LOW_T_BETA = 1.0e6

let Sx = spin_ops(1//2)[1], Sy = spin_ops(1//2)[2], Sz = spin_ops(1//2)[3]
    function ed_heisenberg1d_gs_total(N::Int, J::Real)
        bond = J * (kron(Sx, Sx) + kron(Sy, Sy) + kron(Sz, Sz))
        H = chain_hamiltonian(2, N, bond)
        evals = dense_spectrum(H)
        return minimum(evals)
    end

    @testset "Heisenberg1D — Energy{:total}/OBC small-N GS via ED at β=1e6 (補完 after #449)" begin
        for J in (0.5, 1.0, 2.0)
            for N in (2, 3, 4)
                ed_E_total = ed_heisenberg1d_gs_total(N, J)
                verify(
                    Heisenberg1D(),
                    Energy(:per_site),
                    OBC(N);
                    route=:ed_finite_size,
                    independent=ed_E_total / N,
                    at=["N=$(N)"],
                    agree_within=1e-9,
                    refs=[
                        "ED black-box: chain_hamiltonian(2,N, J·(Sx⊗Sx+Sy⊗Sy+Sz⊗Sz)), minimum eigenvalue (small-N gap exact at β=1e6)",
                    ],
                    fetch_kw=(; J=J, beta=LOW_T_BETA),
                )
            end
        end
    end

    # Sanity card: N=2 dimer singlet energy is the textbook -3J/4
    # (the literal claim of the deleted test_heisenberg_dimer.jl).
    @testset "Heisenberg1D — N=2 OBC dimer GS literature value (補完 after #449)" begin
        for J in (0.5, 1.0, 2.0)
            verify(
                Heisenberg1D(),
                Energy(:per_site),
                OBC(2);
                route=:literature_value,
                independent=-3 * J / 8,
                agree_within=1e-12,
                refs=[
                    "Heisenberg dimer (N=2): unique singlet GS, E0_total = -3J/4 ⇒ E0_per_site = -3J/8",
                ],
                fetch_kw=(; J=J, beta=LOW_T_BETA),
            )
        end
    end
end
