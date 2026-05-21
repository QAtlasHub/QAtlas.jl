# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/XXZ/test_xxz1d_obc_susc_ed_batch.jl
#
# ED-independent verification cards for XXZ1D/Susceptibility{XX,YY,ZZ}/OBC.
# Builds H_XXZ = J Σ (S^x S^x + S^y S^y + Δ S^z S^z) from scratch with
# the canonical spin-1/2 primitives, diagonalises, computes
#   χ_αα = β · Var(M_α) / N,  M_α = Σ σ^α_i
# in the canonical ensemble.  No QAtlas thermal kernel touched.
#
# NOTE (2026-05-20): A preliminary probe shows the current src implementation
# returns values that DIFFER from the ED reference across every (Δ, N, β)
# tested — including a clear SU(2) violation at Δ=1.0 where ED gives
# χ_x = χ_y = χ_z but src gives χ_x = χ_y ≠ χ_z. This batch is filed as
# a "bug-surfacing" card: CI failure is the intended diagnostic signal.
# Refs: tracker issue #445 (XXZ1D thermal kernel discrepancy vs brute-force
# ED at finite β).
#
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra

include(joinpath(@__DIR__, "..", "..", "..", "util", "generic_ed.jl"))

let Sx_op = spin_ops(1//2)[1], Sy_op = spin_ops(1//2)[2], Sz_op = spin_ops(1//2)[3]
    sigmax = 2 * Sx_op
    sigmay = 2 * Sy_op
    sigmaz = 2 * Sz_op

    function ed_xxz_chi(N::Int, J::Real, dz::Real, beta::Real, sigma_alpha::AbstractMatrix)
        bond = J * (kron(Sx_op, Sx_op) + kron(Sy_op, Sy_op) + dz * kron(Sz_op, Sz_op))
        H = chain_hamiltonian(2, N, bond)
        M = sum(site_op(sigma_alpha, 2, N, i) for i in 1:N)
        evals, evecs = eigen(Hermitian(Matrix(H)))
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

    @testset "XXZ1D — Susceptibility{XX,YY,ZZ}/OBC vs ED (#381 batch)" begin
        for (J, dz) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0))
            for N in (4, 5)
                for beta in (1.0, 10.0)
                    for (sigma_alpha, qty, axis_name) in (
                        (sigmax, SusceptibilityXX(), "x"),
                        (sigmay, SusceptibilityYY(), "y"),
                        (sigmaz, SusceptibilityZZ(), "z"),
                    )
                        ed_val = ed_xxz_chi(N, J, dz, beta, sigma_alpha)
                        verify(
                            XXZ1D(),
                            qty,
                            OBC(N);
                            route=:ed_finite_size,
                            independent=ed_val,
                            at=["N=$(N)"],
                            agree_within=1e-9,
                            expected_fail=true,  # tracker issue #445 — XXZ thermal kernel discrepancy vs brute-force ED at finite β
                            refs=[
                                "ED black-box: build H_XXZ from scratch with spin_ops(1/2), diagonalise, compute beta*Var(M_alpha)/N (alpha=$(axis_name)). Refs: tracker issue #445 (XXZ1D thermal kernel discrepancy vs brute-force ED at finite β).",
                            ],
                            fetch_kw=(; J=J, Δ=dz, beta=beta),
                        )
                    end
                end
            end
        end
    end
end
