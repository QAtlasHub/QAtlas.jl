# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_energy_pfeuty_batch.jl
#
# Pfeuty/Bogoliubov closed-form verification card for the TFIM ground-
# state energy per site at the three exactly-evaluable points of the
# integral ε₀(J,h) = -(1/π) ∫₀^π √(J² + h² - 2Jh cos k) dk.
# Pure verify(); self-contained, branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM — Energy/Infinite Pfeuty special points (#381 batch)" begin
    # h = 0 (pure ferromagnet): Λ(k) = 2|J|, ε₀ = -|J|.
    # Include a negative-J point so the card distinguishes ε₀ = -|J| from
    # a hub that mistakenly returns -J.
    for J in (-1.0, 0.5, 1.0, 2.0)
        verify(
            TFIM(; J=J, h=0.0),
            Energy(:per_site),
            Infinite();
            route=:second_closed_form,
            independent=(-abs(J)),
            agree_within=1e-9,
            refs=["Pfeuty 1970: at h=0 dispersion is flat Λ(k)=2|J| ⇒ ε₀ = -|J|"],
        )
    end

    # J = 0 (pure transverse field): Λ(k) = 2|h|, ε₀ = -|h|.
    for h in (0.5, 1.0, 2.0)
        verify(
            TFIM(; J=0.0, h=h),
            Energy(:per_site),
            Infinite();
            route=:second_closed_form,
            independent=(-h),
            agree_within=1e-9,
            refs=["Pfeuty 1970: at J=0 dispersion is flat Λ(k)=2h ⇒ ε₀ = -h"],
        )
    end

    # Critical point h = J. The integrand simplifies via the half-angle
    # identity √(2(1−cos k)) = 2|sin(k/2)| (k ∈ [0,π] ⇒ sin(k/2) ≥ 0), so
    #   √(J²+h²−2Jh cos k)|_{h=J} = √(2J²(1−cos k)) = 2|J|·sin(k/2).
    # Then ∫₀^π 2J sin(k/2) dk = 2J·[−2 cos(k/2)]₀^π = 2J·2 = 4J,
    # giving ε₀ = −(1/π)·4J = −4J/π.
    for J in (0.5, 1.0, 2.0)
        verify(
            TFIM(; J=J, h=J),
            Energy(:per_site),
            Infinite();
            route=:second_closed_form,
            independent=-4 * J / π,
            agree_within=1e-9,
            refs=["Pfeuty 1970: at h=J critical point ε₀ = -4J/π (closed-form integral)"],
        )
    end
end
