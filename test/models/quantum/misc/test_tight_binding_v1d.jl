# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: TightBindingV1D (1D spinless fermion t-V chain).
#
# Targeted run (skips Pkg.test()):
#   julia --project=test test/models/quantum/misc/test_tight_binding_v1d.jl
#
# Phase 1 coverage (V = 0 free-fermion closed forms only):
#   • MassGap          = max(0, |μ| - 2t)
#   • FermiVelocity    = 2t · sin(arccos(-μ/(2t)))   for |μ| < 2t
#   • |μ| ≥ 2t         ⇒ FermiVelocity raises DomainError
#   • V ≠ 0            ⇒ MassGap / FermiVelocity raise DomainError (Phase 2)
#   • t ≤ 0            ⇒ constructor raises DomainError
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TightBindingV1D — V=0 free-fermion MassGap (Phase 1)" begin
    @test QAtlas.fetch(TightBindingV1D(), MassGap(), Infinite()) == 0.0      # default μ=0, gapless
    @test QAtlas.fetch(TightBindingV1D(; μ=1.5), MassGap(), Infinite()) == 0.0  # inside band
    @test QAtlas.fetch(TightBindingV1D(; μ=2.0), MassGap(), Infinite()) == 0.0  # band edge
    @test QAtlas.fetch(TightBindingV1D(; μ=3.0), MassGap(), Infinite()) == 1.0  # insulating
    @test QAtlas.fetch(TightBindingV1D(; μ=-5.0), MassGap(), Infinite()) == 3.0
    @test QAtlas.fetch(TightBindingV1D(; t=2.0, μ=5.0), MassGap(), Infinite()) == 1.0
end

@testset "TightBindingV1D — V=0 free-fermion FermiVelocity (Phase 1)" begin
    @test QAtlas.fetch(TightBindingV1D(), FermiVelocity(), Infinite()) ≈ 2.0          # μ=0 → v_F = 2t
    @test QAtlas.fetch(TightBindingV1D(; t=3.0), FermiVelocity(), Infinite()) ≈ 6.0   # v_F = 2·3
    @test QAtlas.fetch(TightBindingV1D(; μ=1.0), FermiVelocity(), Infinite()) ≈ sqrt(3)
    @test_throws DomainError QAtlas.fetch(
        TightBindingV1D(; μ=2.5), FermiVelocity(), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        TightBindingV1D(; μ=-3.0), FermiVelocity(), Infinite()
    )
end

@testset "TightBindingV1D — V≠0 throws DomainError (Phase 2 deferral)" begin
    @test_throws DomainError QAtlas.fetch(TightBindingV1D(; V=0.5), MassGap(), Infinite())
    @test_throws DomainError QAtlas.fetch(
        TightBindingV1D(; V=-1.0), FermiVelocity(), Infinite()
    )
end

@testset "TightBindingV1D — tiny V (1e-13) is non-zero (iszero strictness)" begin
    # Regression: V boundary must be exact iszero(V), not isapprox(...; atol=1e-12).
    # Any non-zero V puts us on the interacting JW-XXZ branch (Phase 2).
    @test_throws DomainError QAtlas.fetch(TightBindingV1D(; V=1e-13), MassGap(), Infinite())
    @test_throws DomainError QAtlas.fetch(
        TightBindingV1D(; V=-1e-13), FermiVelocity(), Infinite()
    )
end

@testset "TightBindingV1D — rejects t ≤ 0" begin
    @test_throws DomainError TightBindingV1D(; t=0.0)
    @test_throws DomainError TightBindingV1D(; t=-1.0)
end

@testset "TightBindingV1D — V=0 free-fermion Energy{:per_site} (Phase 1)" begin
    # μ=0 (half-filling, J=t=1): e₀ = -2/π
    e_hf = QAtlas.fetch(TightBindingV1D(), Energy(:per_site), Infinite())
    @test isapprox(e_hf, -2 / pi; atol=1e-12)
    @test isapprox(e_hf, -0.6366197723675814; atol=1e-12)

    # Band edges: empty (μ=-2t) → 0; filled (μ=+2t) → -μ = -2
    @test QAtlas.fetch(TightBindingV1D(; μ=-2.0), Energy(:per_site), Infinite()) == 0.0
    @test QAtlas.fetch(TightBindingV1D(; μ=2.0), Energy(:per_site), Infinite()) == -2.0
    # Deep insulating limits
    @test QAtlas.fetch(TightBindingV1D(; μ=-5.0), Energy(:per_site), Infinite()) == 0.0
    @test QAtlas.fetch(TightBindingV1D(; μ=5.0), Energy(:per_site), Infinite()) == -5.0

    # t-linearity at μ=0: e₀(t, 0) = -2t/π
    @test isapprox(
        QAtlas.fetch(TightBindingV1D(; t=3.0), Energy(:per_site), Infinite()),
        -2 * 3 / pi;
        atol=1e-12,
    )

    # Continuity across band edges (one-sided)
    e_just_inside = QAtlas.fetch(
        TightBindingV1D(; μ=2.0 - 1e-9), Energy(:per_site), Infinite()
    )
    @test isapprox(e_just_inside, -2.0; atol=1e-6)
    e_just_empty = QAtlas.fetch(
        TightBindingV1D(; μ=-2.0 + 1e-9), Energy(:per_site), Infinite()
    )
    @test isapprox(e_just_empty, 0.0; atol=1e-6)
end

@testset "TightBindingV1D — Energy{:per_site} V≠0 / tiny-V DomainError (Phase 2 gate)" begin
    @test_throws DomainError QAtlas.fetch(
        TightBindingV1D(; V=0.5), Energy(:per_site), Infinite()
    )
    # iszero(V) strictness regression
    @test_throws DomainError QAtlas.fetch(
        TightBindingV1D(; V=1e-13), Energy(:per_site), Infinite()
    )
    @test_throws DomainError QAtlas.fetch(
        TightBindingV1D(; V=-1e-13), Energy(:per_site), Infinite()
    )
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "TightBindingV1D — verification cards" begin
    # V = 0 free-fermion subspace: half-filling is gapless.
    verify(
        TightBindingV1D(),
        MassGap(),
        Infinite();
        route=:second_closed_form,
        independent=0.0,
        agree_within=1e-10,
        refs=["V=0 half-filled chain: gapless Fermi surface"],
    )
    # μ > 2t band insulator: gap = |μ| - 2t
    verify(
        TightBindingV1D(; μ=3.0),
        MassGap(),
        Infinite();
        route=:second_closed_form,
        independent=1.0,
        agree_within=1e-9,
        refs=["V=0, μ>2t: band-insulator gap = |μ| - 2t"],
    )
end
# ── additional verification cards (#381 batch) ─────────────────────────────
@testset "TightBindingV1D — additional verification cards (#381 batch)" begin
    # V=0 free-fermion v_F = 2 t sin(k_F); default t=1, μ=0 ⇒ v_F = 2.
    verify(
        TightBindingV1D(),
        FermiVelocity(),
        Infinite();
        route=:second_closed_form,
        independent=2.0,
        agree_within=1e-12,
        refs=["Ashcroft-Mermin 1976: v_F = 2 t sin(k_F); V=0 half-filling"],
    )
    # V=0 free-fermion e₀ = -(2t/π) sin(k_F) - (μ/π) k_F; default ⇒ -2/π.
    verify(
        TightBindingV1D(),
        Energy(),
        Infinite();
        route=:second_closed_form,
        independent=-2 / π,
        agree_within=1e-12,
        refs=["Mahan 2000; Ashcroft-Mermin 1976: e₀ = -(2t/π) sin(k_F) - (μ/π) k_F"],
    )
end

@testset "TightBindingV1D — finite-T thermodynamics" begin
    # V = 0 free-fermion finite-T at Infinite.  Cards mirror those for
    # TightBinding1D (same integrand); V ≠ 0 / β ≤ 0 traps stay as
    # @test_throws since they probe exception shape rather than values.

    for (t, β) in [(1.0, 1e-3), (2.5, 5e-4)]
        verify(
            QAtlas.TightBindingV1D(; t=t, V=0.0, μ=0.0),
            QAtlas.FreeEnergy(),
            QAtlas.Infinite();
            route=:limiting_case,
            fetch_kw=(; beta=β),
            independent=(-log(2) / β),
            agree_within=abs(log(2) / β) * 2e-3,
            refs=[
                "Mahan, Many-Particle Physics §1.3: V=0 free-fermion β → 0⁺ limit ω → -T log 2",
            ],
        )
    end

    for (t, μ, β) in [(1.0, 0.0, 1e-3), (1.0, 0.5, 1e-3)]
        verify(
            QAtlas.TightBindingV1D(; t=t, V=0.0, μ=μ),
            QAtlas.ThermalEntropy(),
            QAtlas.Infinite();
            route=:limiting_case,
            fetch_kw=(; beta=β),
            independent=log(2),
            agree_within=log(2) * 4e-3,
            refs=[
                "Mahan, Many-Particle Physics §1.3: V=0 free-fermion β → 0⁺ entropy → log 2 per mode",
            ],
        )
    end

    let β = 1e-2, t = 1.0
        bound = (β * 2 * t)^2
        verify(
            QAtlas.TightBindingV1D(; t=t, V=0.0, μ=0.0),
            QAtlas.SpecificHeat(),
            QAtlas.Infinite();
            route=:limiting_case,
            fetch_kw=(; beta=β),
            independent=0.0,
            agree_within=bound + 1e-12,
            refs=["Mahan, Many-Particle Physics §1.3: V=0 c_μ ~ β² · ⟨ε²⟩/4 → 0 at β → 0⁺"],
        )
    end

    verify(
        QAtlas.TightBindingV1D(),
        QAtlas.FreeEnergy(),
        QAtlas.Infinite();
        route=:limiting_case,
        fetch_kw=(; beta=200.0),
        independent=-2 / π,
        agree_within=5e-3,
        refs=[
            "Ashcroft-Mermin (1976) Ch 9: half-filling 1D free-fermion E/N = -2/π = lim_{β→∞} ω(β) (V=0)",
        ],
    )

    @testset "V ≠ 0 raises DomainError on finite-T quantities" begin
        m = QAtlas.TightBindingV1D(; t=1.0, V=1.0, μ=0.0)
        @test_throws DomainError QAtlas.fetch(
            m, QAtlas.FreeEnergy(), QAtlas.Infinite(); beta=1.0
        )
        @test_throws DomainError QAtlas.fetch(
            m, QAtlas.ThermalEntropy(), QAtlas.Infinite(); beta=1.0
        )
        @test_throws DomainError QAtlas.fetch(
            m, QAtlas.SpecificHeat(), QAtlas.Infinite(); beta=1.0
        )
    end

    @testset "DomainError on β ≤ 0" begin
        @test_throws DomainError QAtlas.fetch(
            QAtlas.TightBindingV1D(), QAtlas.FreeEnergy(), QAtlas.Infinite(); beta=0.0
        )
        @test_throws DomainError QAtlas.fetch(
            QAtlas.TightBindingV1D(), QAtlas.ThermalEntropy(), QAtlas.Infinite(); beta=-1.0
        )
        @test_throws DomainError QAtlas.fetch(
            QAtlas.TightBindingV1D(), QAtlas.SpecificHeat(), QAtlas.Infinite(); beta=0.0
        )
    end
end
