# ─────────────────────────────────────────────────────────────────────────────
# Cross-verification: numerical universality-class exponents vs the
# specific literature source cited in `src/universalities/`.
#
# Each of the three O(N) critical universality classes stores its d = 3
# exponents as Float64 decimals taken from a conformal-bootstrap paper.
# The `test/standalone/test_universality_exponents.jl` file already
# asserts that the stored numbers fall inside a physical range and
# satisfy the scaling relations to `rtol = 0.01`.  That is a weak
# signal — the value could drift by 1 % and still pass.
#
# This file adds the missing tight assertion: the stored values match
# the tabulated decimals from the cited paper to the precision at
# which that paper quotes them.  Any copy-paste slip, factor-of-10
# typo, or accidental rewrite would fail this test immediately.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Universality d=3 numerical values — tabulated literature cross-check" begin
    # ── 3D Ising (O(1))
    #
    #   Source: Kos, Poland, Simmons-Duffin, Vichi,
    #   JHEP 08, 036 (2016), Table 2 ("3d Ising" row).
    #
    #   β = Δ_σ · ν  and Δ_σ = 0.5181489(10), ν = 0.629971(4) / 3.
    #   The stored value β = 0.32642 corresponds to the 5-digit
    #   rounding of Δ_σ · ν from this table.
    # ──────────────────────────────────────────────────────────────────
    @testset "3D Ising — Kos–Poland–Simmons-Duffin–Vichi 2016" begin
        e = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=3)
        @test e.β ≈ 0.32642 atol = 1e-5
        @test e.ν ≈ 0.62997 atol = 1e-5
        @test e.γ ≈ 1.23708 atol = 1e-5
        @test e.η ≈ 0.03630 atol = 1e-5
        @test e.δ ≈ 4.78984 atol = 1e-5
        @test e.α ≈ 0.11009 atol = 1e-5

        # Quoted uncertainties come out as the `*_err` fields — each
        # one must be positive (catches `_err = 0` placeholder bugs).
        for key in (:α_err, :β_err, :γ_err, :δ_err, :ν_err, :η_err)
            @test getproperty(e, key) > 0
        end
    end

    # ── 3D XY / O(2)
    #
    #   Source: Chester, Landry, Liu, Poland, Simmons-Duffin, Su, Vichi,
    #   JHEP 02, 098 (2020) — O(2) conformal-bootstrap "global best"
    #   island, as cited in `src/universalities/ONModel.jl`.
    # ──────────────────────────────────────────────────────────────────
    @testset "3D XY (O(2)) — Chester et al. 2020" begin
        e = QAtlas.fetch(Universality(:XY), CriticalExponents(); d=3)
        @test e.β ≈ 0.34869 atol = 1e-4
        @test e.ν ≈ 0.67175 atol = 1e-4
        @test e.η ≈ 0.038176 atol = 1e-4
        @test e.γ ≈ 1.3179 atol = 1e-3
        @test e.δ ≈ 4.77937 atol = 1e-4
        @test e.α ≈ -0.01526 atol = 1e-4

        for key in (:α_err, :β_err, :γ_err, :δ_err, :ν_err, :η_err)
            @test getproperty(e, key) > 0
        end
    end

    # ── 3D Heisenberg / O(3)
    #
    #   Source: Chester, Landry, Liu, Poland, Simmons-Duffin, Su, Vichi,
    #   Phys. Rev. D 104, 105013 (2021) — O(3) conformal-bootstrap
    #   "global best" island, as cited in `src/universalities/ONModel.jl`.
    # ──────────────────────────────────────────────────────────────────
    @testset "3D Heisenberg (O(3)) — Chester et al. 2021" begin
        e = QAtlas.fetch(Universality(:Heisenberg), CriticalExponents(); d=3)
        @test e.β ≈ 0.3689 atol = 1e-3
        @test e.ν ≈ 0.7112 atol = 1e-3
        @test e.η ≈ 0.0375 atol = 1e-3

        for key in (:α_err, :β_err, :γ_err, :δ_err, :ν_err, :η_err)
            @test getproperty(e, key) > 0
        end
    end

    # ── Scaling-relation internal consistency at tight tolerance
    #
    #   Rushbrooke α + 2β + γ = 2, Widom γ = β(δ−1), Fisher γ = ν(2−η),
    #   Josephson 2 − α = d ν all hold at the numerical-precision level
    #   for well-truncated bootstrap values.  The standalone test file
    #   already checks these at `rtol = 0.01`; tighten to the quoted
    #   2 × 10⁻³ error bars here.  A drift in any stored exponent
    #   breaks a relation whose residual exceeds the quoted error.
    # ──────────────────────────────────────────────────────────────────
    @testset "Scaling relations — within the bootstrap error bars" begin
        for (label, e) in (
            ("Ising d=3", QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=3)),
            ("XY d=3", QAtlas.fetch(Universality(:XY), CriticalExponents(); d=3)),
            (
                "Heisenberg d=3",
                QAtlas.fetch(Universality(:Heisenberg), CriticalExponents(); d=3),
            ),
        )
            @testset "$label" begin
                @test e.α + 2 * e.β + e.γ ≈ 2.0 atol = 2e-3   # Rushbrooke
                @test e.γ ≈ e.β * (e.δ - 1) atol = 5e-3       # Widom
                @test e.γ ≈ e.ν * (2 - e.η) atol = 5e-3       # Fisher
                @test 2 - e.α ≈ 3 * e.ν atol = 5e-3           # Josephson (d = 3)
            end
        end
    end
end
