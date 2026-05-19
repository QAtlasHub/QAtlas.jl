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

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "Heisenberg1D GroundStateEnergyDensity — verification cards" begin
    Sx, Sy, Sz = spin_ops(1 // 2)
    bond_heis = kron(Sx, Sx) + kron(Sy, Sy) + kron(Sz, Sz)
    Ns = verify_profile_Ns(; fast=(6, 8), full=(6, 8, 10, 12), nightly=(6, 8, 10, 12, 14))
    ind = [dense_spectrum(chain_hamiltonian(2, N, bond_heis))[1] / (N - 1) for N in Ns]

    verify(
        Heisenberg1D(),
        GroundStateEnergyDensity(),
        Infinite();
        route=:ed_finite_size,
        independent=ind,
        at=["N=$N" for N in Ns],
        agree_within=0.05,
        refs=["Hulthen 1938: e0 = J(1/4 - log 2), verified by OBC ED convergence"],
    )
end
