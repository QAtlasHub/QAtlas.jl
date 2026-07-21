# Derived-input suppliers (src/core/derivative.jl) + the AD package extensions.
#
# The point of this file is not that ForwardDiff can differentiate x^2.  It is
# that a derivative taken THROUGH a real `fetch` reproduces a thermodynamic
# identity the atlas does not otherwise check — `S = -∂F/∂T` — and that the
# accuracy class of the answer is carried honestly by the backend rather than
# assumed by the caller.

using QAtlas, Test, ForwardDiff
using QAtlas:
    FiniteDifference,
    ForwardDiffBackend,
    ZygoteBackend,
    derivative,
    backend_available,
    preferred_backend,
    default_rtol

@testset "backends: availability and the accuracy class they carry" begin
    # FD needs nothing; ForwardDiff is live because this file loaded it.
    @test backend_available(FiniteDifference())
    @test backend_available(ForwardDiffBackend())
    # Zygote is NOT in the test target, so its extension must stay dormant —
    # and a dormant backend must refuse loudly rather than quietly degrade to a
    # finite difference at a different accuracy class.
    @test !backend_available(ZygoteBackend())
    @test_throws ErrorException derivative(ZygoteBackend(), x -> x^2, 1.0)

    # The tolerance travels with the method: AD earns the tight one, a central
    # difference does not.  (Measured: FD reproduces the analytic entropy to
    # ~5e-5 relative on CurieWeissIsing — 1e-6 there would be a false failure.)
    @test default_rtol(ForwardDiffBackend()) < default_rtol(FiniteDifference())

    # With ForwardDiff loaded it must win the preference order over FD.
    @test preferred_backend() isa ForwardDiffBackend

    @test_throws ArgumentError FiniteDifference(; h=0.0)
end

@testset "both backends agree on a closed form" begin
    f = x -> sin(3x) * exp(-x / 2)
    exact = 3cos(3.0) * exp(-0.5) - 0.5 * sin(3.0) * exp(-0.5)
    @test derivative(ForwardDiffBackend(), f, 1.0) ≈ exact rtol = 1e-12
    @test derivative(FiniteDifference(), f, 1.0) ≈ exact rtol = 1e-6
end

# ── The reason this machinery exists ──────────────────────────────────
# S = -∂F/∂T through a real fetch.  TFIM at Infinite is the reference hub: it
# implements both participants and has no validity-window guard on this sweep.
@testset "S = -∂F/∂T through a real fetch (the derived input relations need)" begin
    m = TFIM(; J=1.0, h=0.5)
    bc = Infinite()
    for T in (0.7, 1.0, 1.6)
        F(t) = QAtlas.fetch(m, FreeEnergy(), bc; beta=1 / t)
        S = QAtlas.fetch(m, ThermalEntropy(), bc; beta=1 / T)

        ad = -derivative(ForwardDiffBackend(), F, T)
        fd = -derivative(FiniteDifference(), F, T)

        # AD is held to its own tolerance, FD to its own — the whole point of
        # `default_rtol` being a property of the backend.
        @test ad ≈ S rtol = default_rtol(ForwardDiffBackend())
        @test fd ≈ S rtol = default_rtol(FiniteDifference())
        # ...and AD is strictly the better estimate, which is why it is preferred.
        @test abs(ad - S) ≤ abs(fd - S) + 1e-12
    end
end
