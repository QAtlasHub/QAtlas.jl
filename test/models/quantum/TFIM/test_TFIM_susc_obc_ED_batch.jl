# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_susc_obc_ED_batch.jl
#
# ED-independent structural corroboration for TFIM/Susceptibility{XX,YY,ZZ}/OBC.
# Builds H_TFIM = -J Σ σ^z σ^z - h Σ σ^x from scratch with the canonical
# spin-1/2 primitives, diagonalises, computes χ_αα = β·Var(M_α)/N in the
# canonical ensemble. No QAtlas thermal kernel touched.
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra

include(joinpath(@__DIR__, "..", "..", "..", "util", "generic_ed.jl"))

let Sx = spin_ops(1//2)[1], Sy = spin_ops(1//2)[2], Sz = spin_ops(1//2)[3]
    sigmax = 2 * Sx
    sigmay = 2 * Sy
    sigmaz = 2 * Sz

    function ed_tfim_chi(N::Int, J::Real, h::Real, beta::Real, sigma_alpha::AbstractMatrix)
        bond = -J * kron(sigmaz, sigmaz)
        onsite = -h * sigmax
        H = chain_hamiltonian(2, N, bond; onsite=onsite)
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

    @testset "TFIM — Susceptibility{XX,YY,ZZ}/OBC vs ED (#381 batch)" begin
        for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0))
            for N in (4, 6, 8)
                for beta in (1.0, 5.0)
                    for (sigma_alpha, qty, axis_name) in (
                        (sigmax, SusceptibilityXX(), "x"),
                        (sigmay, SusceptibilityYY(), "y"),
                        (sigmaz, SusceptibilityZZ(), "z"),
                    )
                        ed_val = ed_tfim_chi(N, J, h, beta, sigma_alpha)
                        verify(
                            TFIM(; J=J, h=h),
                            qty,
                            OBC(N);
                            route=:ed_finite_size,
                            independent=ed_val,
                            at=["N=$(N)"],
                            agree_within=1e-9,
                            refs=[
                                "ED black-box: build H_TFIM from scratch with spin_ops(1/2), diagonalise, compute beta*Var(M_alpha)/N (alpha=$(axis_name))",
                            ],
                            fetch_kw=(; beta=beta),
                        )
                    end
                end
            end
        end
    end
end
