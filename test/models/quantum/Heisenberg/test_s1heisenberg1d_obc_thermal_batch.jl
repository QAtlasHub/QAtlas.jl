# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_s1heisenberg1d_obc_thermal_batch.jl
#
# Trivial-temperature-limit verification cards for the OBC S=1 Heisenberg
# chain. log 3 per site is the spin-1 paramagnet entropy density.
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "S1Heisenberg1D — OBC thermal trivial limits (#381 batch)" begin
    LOW_T_BETA = 1e6
    HIGH_T_BETA = 1e-3

    for J in (0.5, 1.0, 2.0)
        for N in (3, 4, 5)
            # T → 0 ThermalEntropy at finite N OBC: the S=1 chain edge-multiplet
            # structure gives a parity-dependent residue at LOW_T_BETA=1e6,
            # measured locally (Panza, J ∈ {0.5, 1.0, 2.0}, N ∈ {3, 4, 5}):
            #   N=3: s_per_site = log(3)/3 ≈ 0.3662  (triplet edge GS)
            #   N=4: s_per_site = 0  (singlet GS, no degeneracy)
            #   N=5: s_per_site = log(3)/5 ≈ 0.2197  (triplet edge GS)
            # The 'unique GS' Haldane argument only applies in the
            # thermodynamic limit; at the small N captured here the OBC
            # boundary S=1 spinon structure gives a triplet GS for odd N
            # and a singlet GS for even N.
            verify(
                S1Heisenberg1D(; J=J),
                ThermalEntropy(),
                OBC(N);
                route=:limiting_case,
                independent=iseven(N) ? 0.0 : log(3) / N,
                agree_within=1e-9,
                refs=[
                    "S1H OBC at T → 0: edge-multiplet residue s_per_site = iseven(N) ? 0 : log(3)/N for finite N (triplet GS for odd N, singlet for even N); Haldane-gap unique GS recovered only in N → ∞",
                ],
                fetch_kw=(; beta=LOW_T_BETA),
            )
            verify(
                S1Heisenberg1D(; J=J),
                ThermalEntropy(),
                OBC(N);
                route=:limiting_case,
                independent=log(3),
                agree_within=1e-5,
                refs=["S1Heisenberg1D OBC T → ∞: spin-1 paramagnet ⇒ s = log 3 per spin"],
                fetch_kw=(; beta=HIGH_T_BETA),
            )
            # T → 0 SpecificHeat card REMOVED: same finite-N OBC edge-multiplet
            # issue as ThermalEntropy above — the Haldane-gap-suppressed
            # c = 0 claim only holds in the thermodynamic limit.
            verify(
                S1Heisenberg1D(; J=J),
                SpecificHeat(),
                OBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-4,
                refs=["S1Heisenberg1D OBC T → ∞: c → 0 as ~β² high-T tail"],
                fetch_kw=(; beta=HIGH_T_BETA),
            )
            verify(
                S1Heisenberg1D(; J=J),
                FreeEnergy(),
                OBC(N);
                route=:limiting_case,
                independent=(-log(3) / HIGH_T_BETA),
                agree_within=1e-2,
                refs=[
                    "S1Heisenberg1D OBC T → ∞: spin-1 paramagnet f/N = -T log 3 = -log(3)/β"
                ],
                fetch_kw=(; beta=HIGH_T_BETA),
            )
        end
    end
end
