# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_susc_pbc_ED_batch.jl
#
# ED-independent verify card for TFIM/SusceptibilityXX/PBC.
# Builds the PBC TFIM Hamiltonian on a ring of N sites from scratch
# (chain_hamiltonian_pbc + onsite -h σ_x), diagonalises, computes
# χ_xx = β·Var(M_x)/N.
#
# NOTE (2026-05-20, ED-verify-first policy): preliminary probe shows src
# returns values that DIFFER from ED for PBC at low T (β ≥ 1) — e.g. at
# N=8, J=1, h=0.5, β=5: src ≈ 0.56 vs ED ≈ 5.0 (an order-of-magnitude
# discrepancy). OBC matches ED perfectly; only PBC is off. Likely a
# Jordan-Wigner / fermion-parity sector handling issue in the BdG PBC
# path. Bug-surfacing card per ED-verify-first policy.
# Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra

include(joinpath(@__DIR__, "..", "..", "..", "util", "generic_ed.jl"))

let Sx = spin_ops(1//2)[1],
    Sz = spin_ops(1//2)[3]
    sigmax = 2 * Sx
    sigmaz = 2 * Sz

    function ed_tfim_pbc_chi_xx(N::Int, J::Real, h::Real, beta::Real)
        H = Matrix(chain_hamiltonian_pbc(2, N, [(-J * sigmaz, sigmaz)]))
        for i in 1:N
            H .+= -h * site_op(sigmax, 2, N, i)
        end
        M = sum(site_op(sigmax, 2, N, i) for i in 1:N)
        evals, evecs = eigen(Hermitian(H))
        emin = minimum(evals)
        w = exp.(-beta .* (evals .- emin))
        Z = sum(w)
        Md = evecs' * M * evecs
        diagM = real.(diag(Md))
        diagM2 = real.(diag(Md * Md))
        M1 = sum(diagM .* w) / Z
        M2 = sum(diagM2 .* w) / Z
        return beta * (M2 - M1^2) / N
    end

    @testset "TFIM — SusceptibilityXX/PBC vs ED (bug-surfacing) (#381 batch)" begin
        for (J, h) in ((1.0, 0.5), (1.0, 2.0))
            for N in (4, 6, 8)
                for beta in (1.0, 5.0)
                    ed_val = ed_tfim_pbc_chi_xx(N, J, h, beta)
                    verify(
                        TFIM(; J=J, h=h),
                        SusceptibilityXX(),
                        PBC(N);
                        route=:ed_finite_size,
                        independent=ed_val,
                        at=["N=$(N)"],
                        agree_within=1e-9,
                        refs=["ED black-box: build PBC H_TFIM with chain_hamiltonian_pbc + onsite -h σ_x, diagonalise, compute β·Var(M_x)/N"],
                        fetch_kw=(; beta=beta),
                    )
                end
            end
        end
    end
end
