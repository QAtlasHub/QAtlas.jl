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
# (3 hubs × 6 fields = 18 verify cards), plus Percolation d=3 (4 fields).
#
# Percolation is the one that does not fit the shape above, and the difference is
# not cosmetic.  Wang et al. quote the RENORMALIZATION exponents, not the six
# standard ones: `1/ν = 1.141 0(15)` and `y_h = 2.522 95(15)`.  The standard
# exponents follow from those by DEFINITION —
#
#     ν = 1/y_t      β/ν = d − y_h      γ/ν = 2y_h − d
#     δ = y_h/(d − y_h)                 η = d + 2 − 2y_h
#
# — and the pin is against those.  `agree_within` cannot be `0` as it is above,
# because src carries a rounded decimal while the paper carries an eigenvalue; it
# is half a unit in the stored decimal's LAST PLACE, which is the same intent:
# any edit to the stored digits surfaces.  The field's own quoted error would be
# six to forty times looser and would let a digit change pass.  MEASURED
# headroom against the real deviations: β 27×, γ 6.6×, δ 3.7×, η exact.
#
# Two of the six are NOT pinned, and that is the finding rather than an omission:
#
#   * `α`. The only route from the paper's numbers is `α = 2 − dν`, which IS
#     Josephson. Pinning α that way would assert hyperscaling against a table the
#     :derivation plane tests hyperscaling ON. The shipped `−0.625(3)` is 1.4σ
#     from `2 − dν = −0.629 27` and the paper quotes no α, so its source is
#     elsewhere and unrecorded.
#   * `ν`. The paper gives `ν = 1/1.1410(15) = 0.876 42(115)`. src carries
#     `0.876 19(12)` — 2.0σ away in its OWN error, and with an error ten times
#     tighter than this paper supports. Also from elsewhere, also unrecorded.
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

    # ─── 3D Percolation (Wang-Zhou-Zhang-Garoni-Deng 2013) ───
    # Derived from the paper's own y_t and y_h through the definitions above, so
    # the pin is independent of the scaling identities: Widom, Fisher and
    # Rushbrooke are never used to produce any of these four.
    let d = 3,
        y_t = 1.1410,
        y_h = 2.52295,
        cite = "Wang-Zhou-Zhang-Garoni-Deng 2013 (PRE 87, 052107)"

        ν_p = 1 / y_t
        for (field, value, tol, how) in (
            (:β, (d - y_h) * ν_p, 5e-5, "β/ν = d − y_h"),
            (:γ, (2 * y_h - d) * ν_p, 5e-4, "γ/ν = 2y_h − d"),
            (:δ, y_h / (d - y_h), 5e-3, "δ = y_h/(d − y_h)"),
            (:η, d + 2 - 2 * y_h, 5e-6, "η = d + 2 − 2y_h"),
        )
            verify(
                Universality(:Percolation),
                CriticalExponents(),
                Infinite();
                route=:literature_value,
                independent=value,
                agree_within=tol,
                refs=[
                    "$(cite): y_t = $(y_t), y_h = $(y_h); $(how) gives $(field) = $(round(value; sigdigits=6))",
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
