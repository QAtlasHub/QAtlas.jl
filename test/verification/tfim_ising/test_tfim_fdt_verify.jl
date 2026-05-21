# ─────────────────────────────────────────────────────────────────────────────
# TFIM FDT — verify() cards (WHY-correct plane)
#
# Split out of test/verification/tfim_ising/test_tfim_fdt.jl (6.2 min on s03)
# so each beta slice runs on its own shard. Same FDT physics:
#     S_zz(q, omega; beta) = (2 / (1 - exp(-beta * omega))) * chi''_zz(q, omega; beta)
# Reference: Kubo 1957; Mahan, *Many-Particle Physics*, ch. 3.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM FDT — verification cards" begin
    # The fluctuation-dissipation setup uses TFIM(J=1, h=0.5), gapped with
    # Pfeuty gap Delta = 2|h - J| = 1 (independent closed form).
    verify(
        TFIM(; J=1.0, h=0.5),
        MassGap(),
        Infinite();
        route=:second_closed_form,
        independent=1.0,
        agree_within=1e-10,
        refs=["Pfeuty 1970: Delta = 2|h - J| = 1 at (J=1, h=0.5)"],
    )
end
