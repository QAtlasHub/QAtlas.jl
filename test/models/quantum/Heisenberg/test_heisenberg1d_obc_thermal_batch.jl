# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/Heisenberg/test_heisenberg1d_obc_thermal_batch.jl
#
# Trivial-temperature-limit verification cards for the OBC Heisenberg
# chain:
#   ThermalEntropy : T → 0 ⇒ s = 0 for even N (unique singlet GS), but
#                    s = log(2)/N for odd N (S=1/2 doublet GS at finite N);
#                    T → ∞ ⇒ s = log 2 (paramagnet, 2 states/spin)
#   SpecificHeat   : T → 0 ⇒ c = 0 (gap suppression); T → ∞ ⇒ c → 0 (β² tail)
#   FreeEnergy     : T → ∞ ⇒ f/N → -T log 2 = -log(2)/β (paramagnet)
#
# Tolerance notes:
#   * FreeEnergy high-T tolerance is scaled as 5e-3 * J^2 to absorb the
#     leading cumulant correction β · ⟨H²⟩₀ / N / 2 ≈ β · J² · O(1),
#     which is ~2e-3 at J=2, β=1e-3. This keeps the margin explicit
#     instead of relying on a fixed 1e-2 that would silently break for
#     larger J.
#
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Heisenberg1D — OBC thermal trivial limits (#381 batch)" begin
    LOW_T_BETA = 1e6
    HIGH_T_BETA = 1e-3

    for J in (0.5, 1.0, 2.0)
        for N in (4, 6, 8)
            # ThermalEntropy
            # T → 0: for even N the OBC AFM Heisenberg chain has a unique
            # singlet GS ⇒ s = 0. For odd N the GS is an S=1/2 doublet ⇒
            # s = log(2)/N as the residual zero-point entropy per site.
            verify(
                Heisenberg1D(),
                ThermalEntropy(),
                OBC(N);
                route=:limiting_case,
                independent=iseven(N) ? 0.0 : log(2.0) / N,
                agree_within=1e-9,
                refs=[
                    "Heisenberg1D OBC T → 0: singlet GS for even N ⇒ s = 0; doublet GS for odd N ⇒ s = log(2)/N",
                ],
                fetch_kw=(; J=J, beta=LOW_T_BETA),
            )
            verify(
                Heisenberg1D(),
                ThermalEntropy(),
                OBC(N);
                route=:limiting_case,
                independent=log(2),
                agree_within=1e-5,
                refs=["Heisenberg1D OBC T → ∞: free paramagnet ⇒ s = log 2 per spin"],
                fetch_kw=(; J=J, beta=HIGH_T_BETA),
            )
            # SpecificHeat
            verify(
                Heisenberg1D(),
                SpecificHeat(),
                OBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-9,
                refs=["Heisenberg1D OBC T → 0: gap suppression ⇒ c = 0 exactly"],
                fetch_kw=(; J=J, beta=LOW_T_BETA),
            )
            verify(
                Heisenberg1D(),
                SpecificHeat(),
                OBC(N);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-4,
                refs=["Heisenberg1D OBC T → ∞: c → 0 as ~β² (high-T tail)"],
                fetch_kw=(; J=J, beta=HIGH_T_BETA),
            )
            # FreeEnergy at T → ∞: f/N → -T log 2 = -log(2)/β
            # Tolerance scales as ~5e-3 * J^2 to track the β·J²·O(1)
            # leading cumulant correction (≈ 2e-3 at J=2, β=1e-3).
            verify(
                Heisenberg1D(),
                FreeEnergy(),
                OBC(N);
                route=:limiting_case,
                independent=(-log(2) / HIGH_T_BETA),
                agree_within=5e-3 * J^2,
                refs=[
                    "Heisenberg1D OBC T → ∞: free paramagnet ⇒ f/N = -T log 2 = -log(2)/β; tol ~β·J²",
                ],
                fetch_kw=(; J=J, beta=HIGH_T_BETA),
            )
        end
    end
end
