# ─────────────────────────────────────────────────────────────────────────────
# test/universalities/test_universality_critical_exponents_lit.jl
#
# Literature-value pins for Universality{X}/CriticalExponents/Infinite: the d=3
# conformal-bootstrap exponents (Ising, XY, Heisenberg), the d=3 percolation
# Monte-Carlo table, and the d=2 Coulomb-gas family (percolation, Potts q=3,
# Potts q=4).
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
# Percolation d=3 is pinned against DEFINITIONS of the paper's y_t and y_h, and
# leaves α and ν unpinned; Percolation.jl's own header says why.  Its
# `agree_within` is half a unit in the stored decimal's last place rather than
# `0`, because src rounds where the paper does not.
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

# ─────────────────────────────────────────────────────────────────────────────
# d=2 Coulomb gas: percolation (q=1), Potts q=3, Potts q=4
#
# One closed form covers all three ([Xu2025](@cite) Eqs. 20, 22a, 23a, coded
# below), with g read off that paper's Table I.  The pin is EXACT rational
# equality: these tables are exact, so there is no rounding to hide in.  alpha is
# excluded as in the d=3 block, its only route being Josephson.
const _CG_G = Dict(
    :percolation => 2 // 3, :ising => 3 // 4, :potts3 => 5 // 6, :potts4 => 1 // 1
)
_cg_yt(g) = 3 * (2g - 1) // (2g)
_cg_yh(g) = (2g + 1) * (2g + 3) // (8g)
function _cg_exponents(g; d=2)
    y_t, y_h = _cg_yt(g), _cg_yh(g)
    return (
        β=(d - y_h) // y_t,
        γ=(2y_h - d) // y_t,
        δ=y_h // (d - y_h),
        ν=1 // y_t,
        η=d + 2 - 2y_h,
        α=2 - d // y_t,
    )
end

@testset "d=2 Coulomb gas reproduces every exact exponent this atlas ships" begin
    # The paper's own Table I, so the formulas are checked before they are used.
    for (q, g, y_t, y_h) in (
        (1, 2 // 3, 3 // 4, 91 // 48),
        (2, 3 // 4, 1 // 1, 15 // 8),
        (3, 5 // 6, 6 // 5, 28 // 15),
        (4, 1 // 1, 3 // 2, 15 // 8),
    )
        @test _CG_G[(:percolation, :ising, :potts3, :potts4)[q]] == g
        @test 4 * cos(pi * float(g))^2 ≈ q atol = 1e-12
        @test _cg_yt(g) == y_t
        @test _cg_yh(g) == y_h
    end

    # alpha is included here and excluded from the pins below: the closed form
    # produces it, so a route that also produces it cannot check it.
    for (key, M, kw) in (
        (:percolation, Universality(:Percolation), (; d=2)),
        (:ising, Universality(:Ising), (; d=2)),
        (:potts3, Universality(:Potts3), (; d=2)),
        (:potts4, Universality(:Potts4), (; d=2)),
    )
        shipped = QAtlas.fetch(M, CriticalExponents(); kw...)
        closed = _cg_exponents(_CG_G[key])
        for f in (:α, :β, :γ, :δ, :ν, :η)
            @test shipped[f] == closed[f]
        end
    end

    # Positive control: the equality is a claim about THESE g, not something the
    # formulas satisfy for any argument.  One step along the g axis breaks it.
    off = _cg_exponents(_CG_G[:potts3] + 1 // 12)
    @test any(
        off[f] != QAtlas.fetch(Universality(:Potts3), CriticalExponents(); d=2)[f] for
        f in (:α, :β, :γ, :δ, :ν, :η)
    )
end

@testset "d=2 Coulomb-gas literature pins" begin
    for (key, M, label) in (
        (:percolation, Universality(:Percolation), "q=1 (percolation)"),
        (:potts3, Universality(:Potts3), "q=3"),
        (:potts4, Universality(:Potts4), "q=4"),
    )
        g = _CG_G[key]
        closed = _cg_exponents(g)
        for f in (:β, :γ, :δ, :ν, :η)
            verify(
                M,
                CriticalExponents(),
                Infinite();
                route=:literature_value,
                independent=float(closed[f]),
                agree_within=0,
                at=["d=2", "field=$(f)"],
                refs=[
                    "Xu-Salas-Deng 2025 (Entropy 27, 418) Eqs. (20), (22a), (23a) and Table I: \
                     d=2 Potts $(label) at g = $(g), y_t = $(_cg_yt(g)), y_h = $(_cg_yh(g)) \
                     gives $(f) = $(closed[f])",
                ],
                fetch_kw=(; d=2),
                subject_extract=e -> float(getproperty(e, f)),
            )
        end
    end
end
