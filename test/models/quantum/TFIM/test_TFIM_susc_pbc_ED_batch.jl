# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_susc_pbc_ED_batch.jl
#
# ED-independent verify card for TFIM/SusceptibilityXX/PBC.
# Builds the PBC TFIM Hamiltonian on a ring of N sites from scratch
# (chain_hamiltonian_pbc + onsite -h σ_x), diagonalises, computes
# χ_xx = β·Var(M_x)/N.
#
# NOTE (2026-06-03, post-fix): the original 0.56 vs 5.0 discrepancy at
# N=8, J=1, h=0.5, β=5 reflected an ED reference that used the FDT
# formula β·Var(M)/N — only valid for classical / commuting [H, M].
# TFIM has [H, σ_x] ≠ 0, so the correct quantum static susceptibility
# is the Kubo sum-over-eigenpairs form, which agrees with src to
# machine precision. The thermo (E/F/S/C) half of #444 was a genuine
# JW parity-sector bug fixed in TFIM_pbc_thermal.jl (sign of the
# (R, sinh) sector now depends on |h| vs |J|). See PR for #444.
# Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra

include(joinpath(@__DIR__, "..", "..", "..", "util", "generic_ed.jl"))

let Sx = spin_ops(1//2)[1], Sz = spin_ops(1//2)[3]
    sigmax = 2 * Sx
    sigmaz = 2 * Sz

    function ed_tfim_pbc_chi_xx(N::Int, J::Real, h::Real, beta::Real)
        # χ_xx = ∂⟨M⟩/∂h via Kubo on the full ED spectrum (correct for
        # quantum non-commuting [H, M] ≠ 0). The naive FDT form
        # β·Var(M)/N is only valid when [H, M] = 0, which TFIM violates,
        # so this card uses the proper sum-over-eigenpairs form.
        H = Matrix(chain_hamiltonian_pbc(2, N, [(-J * sigmaz, sigmaz)]))
        for i in 1:N
            H .+= -h * site_op(sigmax, 2, N, i)
        end
        M = sum(site_op(sigmax, 2, N, i) for i in 1:N)
        evals, evecs = eigen(Hermitian(H))
        emin = minimum(evals)
        w = exp.(-beta .* (evals .- emin))
        Z = sum(w)
        p = w ./ Z
        Mab = evecs' * M * evecs
        Mmean = sum(p[m] * real(Mab[m, m]) for m in eachindex(p))
        χ = 0.0
        for m in eachindex(evals), n in eachindex(evals)
            ΔE = evals[m] - evals[n]
            if abs(ΔE) > 1e-10
                χ += (p[n] - p[m]) / ΔE * abs(Mab[m, n])^2
            else
                χ += beta * p[m] * abs(Mab[m, n])^2
            end
        end
        return (χ - beta * Mmean^2) / N
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
                        refs=[
                            "ED black-box: build PBC H_TFIM with chain_hamiltonian_pbc + onsite -h σ_x, diagonalise, compute β·Var(M_x)/N",
                        ],
                        fetch_kw=(; beta=beta),
                    )
                end
            end
        end
    end
end
