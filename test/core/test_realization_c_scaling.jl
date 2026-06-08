# INDEPENDENT verification of the model ↔ universality correspondence: a model's
# OWN entanglement-entropy scaling must reproduce the central charge of the class
# it realizes — computed by a method that has NO knowledge of c.
#
# The honest case here is TFIM: at criticality (h = J) its von Neumann entropy
# comes from Peschel's free-fermion correlation-matrix method (O(ℓ³), exact, and
# entirely c-agnostic). Extracting c from its OBC scaling and matching the
# realized class (:Ising, c = 1/2) genuinely verifies the correspondence — it is
# NOT the Calabrese–Cardy formula compared against itself.
#
# The interacting realizations (XXZ Δ=1, Heisenberg1D, HaldaneShastry) only have
# a closed-form CC entanglement at `Infinite` (which *is* the universal formula),
# so an independent c there needs ED/DMRG — deferred to the verification-density
# work. We verify what can be verified independently, and say so.

using QAtlas, Test
using QAtlas: TFIM, Universality, realized_class

@testset "realization c-scaling — independent (TFIM Peschel)" begin
    m = TFIM(; J=1.0, h=1.0)                       # Ising critical point
    @test realized_class(m) === :Ising

    # OBC half-cut CC scaling: S(ℓ) = (c/6) log[(2N/π) sin(πℓ/N)] + const.
    # Extract c from two subsystem sizes via Peschel (independent of any c input).
    N = 400
    chord(ℓ) = (2N / π) * sin(π * ℓ / N)
    ℓ1, ℓ2 = 100, 200
    s1 = QAtlas.fetch(m, VonNeumannEntropy(), OBC(); ℓ=ℓ1, N=N)
    s2 = QAtlas.fetch(m, VonNeumannEntropy(), OBC(); ℓ=ℓ2, N=N)
    c_est = 6 * (s2 - s1) / (log(chord(ℓ2)) - log(chord(ℓ1)))

    c_class = Float64(QAtlas.fetch(Universality(:Ising), CentralCharge()))
    @test isapprox(c_est, c_class; rtol=0.05)      # ≈ 0.50, matches :Ising
    # crucially, the check actually DISCRIMINATES: it must reject the other
    # common 1+1D value c = 1 (XY / Heisenberg), else it could not tell classes apart.
    @test !isapprox(c_est, 1.0; rtol=0.05)
end
