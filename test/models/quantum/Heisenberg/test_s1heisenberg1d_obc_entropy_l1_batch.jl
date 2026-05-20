# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_s1heisenberg1d_obc_entropy_l1_batch.jl
#
# Single-site entanglement entropy of the OBC S=1 Heisenberg chain: by
# SU(2) symmetry ρ₁ = I/3 (3 m=±1,0 states maximally mixed) so
# S_vN(ℓ=1) = S_α(ℓ=1) = log 3 for ALL J, N, β. (The Haldane edge
# quartet on OBC is itself SU(2)-symmetric — a thermal mixture of
# the four degenerate GS components — so the bulk argument applies
# unchanged at any β.)
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "S1Heisenberg1D — VonNeumannEntropy + RenyiEntropy/OBC at ℓ=1 (#381 batch)" begin
    for J in (0.5, 1.0, 2.0)
        for N in (3, 4, 5)
            for β in (0.5, 10.0, 1e6)
                verify(
                    S1Heisenberg1D(),
                    VonNeumannEntropy(),
                    OBC(N);
                    route=:second_closed_form,
                    independent=log(3),
                    agree_within=1e-10,
                    refs=["S1Heisenberg1D SU(2) symmetry: ρ₁ = I/3 ⇒ S_vN(ℓ=1) = log 3 for all J, N, β"],
                    fetch_kw=(; J=J, ℓ=1, beta=β),
                )
                for α in (2, 3)
                    verify(
                        S1Heisenberg1D(),
                        RenyiEntropy(α),
                        OBC(N);
                        route=:second_closed_form,
                        independent=log(3),
                        agree_within=1e-10,
                        refs=["S1Heisenberg1D SU(2): ρ₁ = I/3 ⇒ S_α(ℓ=1) = log 3 for all Rényi index"],
                        fetch_kw=(; J=J, ℓ=1, beta=β),
                    )
                end
            end
        end
    end
end
