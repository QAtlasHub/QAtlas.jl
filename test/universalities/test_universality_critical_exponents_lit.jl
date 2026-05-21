# ─────────────────────────────────────────────────────────────────────────────
# test/universalities/test_universality_critical_exponents_lit.jl
#
# Literature-value pins for Universality{X}/CriticalExponents/Infinite at d=3
# (3D Ising, XY, Heisenberg conformal-bootstrap exponents).
#
# Restores a piece of the WHY-plane coverage removed in PR #449:
# the deleted test/verification/universality/test_universality_literature_values.jl
# pinned the stored exponent decimals in src/universalities/Ising2D.jl and
# src/universalities/ONModel.jl by raw @test against the cited paper.
# This file does the same via verify() using the new subject_extract hook
# (added in the same PR), so each field becomes its own structural
# literature-pin card with route :literature_value.
#
# Tolerance choice: agree_within=0 — we're pinning the exact decimal
# stored in src, so any drift past the typed precision must surface.
# (The deleted test used atol=1e-5; agree_within=0 is strictly tighter.)
#
# Hubs added: Universality(:Ising/:XY/:Heisenberg)/CriticalExponents/Infinite
# (3 hubs × 6 fields = 18 verify cards).
#
# References:
#   Ising d=3:      Kos, Poland, Simmons-Duffin, Vichi, JHEP 08, 036 (2016).
#   XY d=3:         Chester et al., JHEP 02, 098 (2020).
#   Heisenberg d=3: Chester et al., Phys. Rev. D 104, 105013 (2021).
#
# Pure verify(); branches off main. Refs #381; restores coverage lost in #449.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Universality d=3 CriticalExponents — literature pins (補完 after #449)" begin
    # ─── 3D Ising (Kos-Poland-Simmons-Duffin-Vichi 2016) ───
    let lit = (α=0.11009, β=0.32642, γ=1.23708, δ=4.78984, ν=0.62997, η=0.03630)
        for (field, value) in pairs(lit)
            verify(
                Universality(:Ising),
                CriticalExponents(),
                Infinite();
                route=:literature_value,
                independent=value,
                agree_within=0,
                refs=[
                    "Kos-Poland-Simmons-Duffin-Vichi 2016 (JHEP 08, 036): 3D Ising conformal bootstrap, $(field) = $(value)",
                ],
                fetch_kw=(; d=3),
                subject_extract=e -> getproperty(e, field),
            )
        end
    end

    # ─── 3D XY / O(2) (Chester et al. 2020) ───
    let lit = (α=-0.01526, β=0.34869, γ=1.3179, δ=4.77937, ν=0.67175, η=0.038176)
        for (field, value) in pairs(lit)
            verify(
                Universality(:XY),
                CriticalExponents(),
                Infinite();
                route=:literature_value,
                independent=value,
                agree_within=0,
                refs=[
                    "Chester-Landry-Liu-Poland-Simmons-Duffin-Su-Vichi 2020 (JHEP 02, 098): 3D O(2) conformal bootstrap, $(field) = $(value)",
                ],
                fetch_kw=(; d=3),
                subject_extract=e -> getproperty(e, field),
            )
        end
    end

    # ─── 3D Heisenberg / O(3) (Chester et al. 2021) ───
    let lit = (α=-0.1336, β=0.3689, γ=1.3960, δ=4.783, ν=0.7112, η=0.0375)
        for (field, value) in pairs(lit)
            verify(
                Universality(:Heisenberg),
                CriticalExponents(),
                Infinite();
                route=:literature_value,
                independent=value,
                agree_within=0,
                refs=[
                    "Chester-Landry-Liu-Poland-Simmons-Duffin-Su-Vichi 2021 (Phys. Rev. D 104, 105013): 3D O(3) conformal bootstrap, $(field) = $(value)",
                ],
                fetch_kw=(; d=3),
                subject_extract=e -> getproperty(e, field),
            )
        end
    end
end
