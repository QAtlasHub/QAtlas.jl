# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/misc/test_hubbard1d_gsed_U0_batch.jl
#
# Textbook U → 0 limit of Hubbard1D GSED at half-filling: e_0 = -4t/π
# (Essler et al. 2005 Lieb-Wu chapter Eq. 1.21; matches src docstring at
# Hubbard1D.jl:167). Previously bug-surfacing for issue #423 (src returned
# -4t²/π with an extra factor of t); root cause fixed in this PR by changing
# _hubbard1d_e0 prefactor from -4*t^2*integral to -4*t*integral. Refs #381,
# #390, #423.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Hubbard1D — GSED/Infinite at U → 0 textbook (bug-surfacing) (#381 batch)" begin
    for t in (0.5, 1.0, 2.0)
        U = 1e-4 * t
        verify(
            Hubbard1D(; t=t, U=U, μ=U/2),
            GroundStateEnergyDensity(),
            Infinite();
            route=:limiting_case,
            independent=-4 * t / π,
            agree_within=1e-3,
            refs=[
                "Essler et al. 2005 (Lieb-Wu Eq. 1.21): e_0 = -4t/π at U → 0 half-filling (two decoupled tight-binding chains)",
            ],
        )
    end
end

# Large-U asymptote: e_0 → -4 t² log 2 / U (Heisenberg-AFM reduction;
# Lieb-Wu 1968, Essler et al. 2005). t ≠ 1 is mandatory to distinguish
# the t prefactor from t² (at t=1, t = t² = 1 and the test is degenerate).
# agree_within = 1e-4 covers the leading sub-leading correction
# (relative ~ t/U ~ 1% at U/t = 100, abs ~ 5e-5 at the worst sweep point).
@testset "Hubbard1D — GSED/Infinite at U → ∞ Heisenberg asymptote (#381 batch)" begin
    for t in (0.5, 2.0)
        U = 100.0 * t
        verify(
            Hubbard1D(; t=t, U=U, μ=U/2),
            GroundStateEnergyDensity(),
            Infinite();
            route=:limiting_case,
            independent=-4 * t^2 * log(2) / U,
            agree_within=1e-4,
            refs=[
                "Lieb-Wu 1968 / Essler et al. 2005: e_0 → -4t² log(2) / U at U → ∞ (Heisenberg-AFM reduction); t ≠ 1 distinguishes the t prefactor from t²",
            ],
        )
    end
end
