# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: Bethe ansatz ground-state energy density
#
# Verify the Hulthén (1938) result e₀ = J(1/4 − ln 2) for the
# thermodynamic-limit Heisenberg chain, and its consistency with the
# finite-size ED spectra already stored in QAtlas.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Bethe ansatz: e₀ = J(1/4 − ln 2)" begin
    J = 1.0
    e0 = QAtlas.fetch(Heisenberg1D(), GroundStateEnergyDensity(); J=J)

    # Numerical value
    @test e0 ≈ J * (0.25 - log(2)) atol = 1e-14
    @test e0 ≈ -0.44314718055994530 atol = 1e-12

    # J scaling
    for Jval in (0.5, 2.0, 3.7)
        @test QAtlas.fetch(Heisenberg1D(), GroundStateEnergyDensity(); J=Jval) ≈ Jval * e0 rtol =
            1e-14
    end

    # Consistency with finite-size spectra:
    # E₀(N)/N should overshoot e₀ for finite PBC chains (finite-size
    # correction is negative), i.e., E₀(N)/N < e₀.
    λ4 = QAtlas.fetch(Heisenberg1D(), ExactSpectrum(); N=4, J=J, bc=:PBC)
    @test λ4[1] / 4 < e0  # N=4: E₀/N ≈ -0.5 < -0.443

    # The finite-size error |E₀(N)/N − e₀| should be O(1/N²)
    # For N=4: error ≈ |-0.5 − (-0.443)| ≈ 0.057
    @test abs(λ4[1] / 4 - e0) < 0.1
end
