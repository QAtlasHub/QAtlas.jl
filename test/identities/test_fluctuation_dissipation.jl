# ─────────────────────────────────────────────────────────────────────────────
# Test: the fluctuation–dissipation foundation (test/util/fluctuation_dissipation.jl).
#
# Three layers of self-consistency, each comparing two GENUINELY DIFFERENT
# computations of the same quantity — never the same closed form re-typed:
#
#   LAYER 1 — arbitrary spectra.  For random + physical spectra the ensemble
#   moments (⟨E⟩, Var(E)) equal the β-derivatives of lnZ (ForwardDiff) to
#   machine precision, and the static linear-response FDT ∂⟨O⟩/∂λ = β·Var(O)
#   holds for a diagonal observable: the thermodynamic + FDT relations are
#   internally consistent for ANY spectrum.
#
#   LAYER 2 — a real quantum many-body spectrum.  The TFIM many-body spectrum
#   (dense ED of build_tfim) obeys the energy FDT C = β²Var(E) = -β²∂⟨E⟩/∂β.
#
#   LAYER 3 — atlas tie-in.  QAtlas's registered IsingChain1D SpecificHeat and
#   SusceptibilityZZ equal β²·Var(E)/N and β·Var(M)/N from a brute-force config
#   enumeration (directly, and through the FLUCTUATION_DISSIPATION_IDENTITIES /
#   ThermoIdentity harness) — fluctuation (independent) == registered response.
#
# `thermo_from_spectrum`, `log_partition`, `mean_energy`, `gibbs_moments`, the
# variance providers, the FDT identities and `verify_thermodynamic_identities`
# all come from test/util/ (included by runtests before this file).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Lattice2D, LinearAlgebra, Test, ForwardDiff, Random
using QAtlas: IsingChain1D, SpecificHeat, SusceptibilityZZ, Infinite, fetch

# ── LAYER 1a — arbitrary spectra: ensemble moments == lnZ derivatives ─────────
@testset "FDT layer 1 — arbitrary spectra: thermodynamic + energy FDT self-consistent" begin
    rng = MersenneTwister(0x5eed)
    spectra = [
        ("two-level", [0.0, 1.0]),
        ("equally-spaced (5)", collect(0.0:1.0:4.0)),
        ("harmonic n+½ (12)", [n + 0.5 for n in 0:11]),
        ("degenerate", [0.0, 0.0, 1.0, 1.0, 1.0, 3.5]),
        ("random uniform (40)", sort(8 .* rand(rng, 40))),
        ("random + offset (30)", sort(-4 .+ 10 .* rand(rng, 30))),
        ("wide-scale (20)", sort(40 .* rand(rng, 20))),
    ]
    βs = (0.1, 0.5, 1.0, 2.0)
    for (name, levels) in spectra, β in βs
        th = thermo_from_spectrum(levels, β)

        # (1) mean energy: ensemble average  ⟨E⟩  ==  -∂lnZ/∂β  (derivative)
        dlnZ = ForwardDiff.derivative(b -> log_partition(levels, b), β)
        @test isapprox(th.E, -dlnZ; rtol=1e-9, atol=1e-10)

        # (2) energy FDT (core):  Var(E)  ==  ∂²lnZ/∂β²  ==  -∂⟨E⟩/∂β
        d2lnZ = ForwardDiff.derivative(
            b -> ForwardDiff.derivative(bb -> log_partition(levels, bb), b), β
        )
        dE = ForwardDiff.derivative(b -> mean_energy(levels, b), β)
        @test isapprox(th.varE, d2lnZ; rtol=1e-6, atol=1e-8)
        @test isapprox(th.varE, -dE; rtol=1e-6, atol=1e-8)

        # (3) specific heat from fluctuations  ==  from temperature response
        @test isapprox(th.C, -β^2 * dE; rtol=1e-6, atol=1e-8)

        # (4) entropy as a free-energy response:  S  ==  -∂F/∂T  =  β²·∂F/∂β
        dF = ForwardDiff.derivative(b -> -log_partition(levels, b) / b, β)
        @test isapprox(th.S, β^2 * dF; rtol=1e-6, atol=1e-8)

        # (5) Gibbs ⟨E⟩ = F + T·S, with S from the DERIVATIVE route (not β(E-F)):
        #     three independent computations (ensemble E, -lnZ/β, ∂F/∂β) reconcile.
        S_deriv = β^2 * dF
        @test isapprox(th.E, th.F + S_deriv / β; rtol=1e-7, atol=1e-9)
    end
end

# ── LAYER 1b — static linear-response FDT for a commuting observable ──────────
@testset "FDT layer 1 — linear response: ∂⟨O⟩/∂λ = β·Var(O)  for  H(λ)=H₀-λO" begin
    rng = MersenneTwister(0xf17)
    # ⟨O⟩ of a diagonal observable `obs` under H(λ) = H₀ - λ·O, at field λ.
    mean_obs(E0, obs, β, λ) = gibbs_moments(E0 .- λ .* obs, obs, β).mean
    cases = [
        ("spin-1 (m=-1,0,1)", [0.0, 0.0, 0.0], [-1.0, 0.0, 1.0]),
        ("two paramagnetic spins", [0.0, 0.5, 0.5, 1.0], [-2.0, 0.0, 0.0, 2.0]),
        ("random levels + random obs (12)", sort(5 .* rand(rng, 12)), randn(rng, 12)),
    ]
    for (name, E0, obs) in cases, β in (0.3, 1.0, 2.5), λ in (-0.7, 0.0, 0.4)
        # response side: ∂⟨O⟩/∂λ via ForwardDiff
        dO = ForwardDiff.derivative(l -> mean_obs(E0, obs, β, l), λ)
        # fluctuation side: β·Var(O) at the same field
        varO = gibbs_moments(E0 .- λ .* obs, obs, β).var
        @test isapprox(dO, β * varO; rtol=1e-7, atol=1e-9)
    end
end

# ── LAYER 2 — real quantum many-body spectrum (TFIM dense ED) obeys energy FDT ─
@testset "FDT layer 2 — TFIM many-body spectrum obeys the energy FDT" begin
    for N in (6, 8), (J, h) in ((1.0, 0.5), (1.0, 1.0), (0.7, 1.3))
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())
        levels = eigvals(Symmetric(build_tfim(lat, J, h)))   # 2ᴺ many-body energies
        @test length(levels) == 1 << N
        for β in (0.3, 0.8, 1.5)
            th = thermo_from_spectrum(levels, β)
            # variance route (2nd moment) vs derivative route — the FDT theorem
            # on a genuine interacting quantum spectrum.
            d2lnZ = ForwardDiff.derivative(
                b -> ForwardDiff.derivative(bb -> log_partition(levels, bb), b), β
            )
            dE = ForwardDiff.derivative(b -> mean_energy(levels, b), β)
            @test isapprox(th.varE, d2lnZ; rtol=1e-7, atol=1e-9)
            @test isapprox(th.C, -β^2 * dE; rtol=1e-7, atol=1e-9)
            @test th.C ≥ 0          # heat capacity is a variance ⇒ non-negative
        end
    end
end

# ── LAYER 3 — atlas tie-in: registered IsingChain1D C, χ == brute-force flux. ──
@testset "FDT layer 3 — IsingChain1D registered C, χ == brute-force fluctuations" begin
    m = IsingChain1D(; J=1.0)
    # High-T β: the energy-variance finite-size error scales as O(N²·tanh(βJ)^N)
    # (N² from differentiating the (λ₋/λ₊)^N transfer-matrix gap twice in β), so
    # at N=16 it is < 2e-5 for βJ ≤ 0.4 but already ~2e-4 by βJ = 0.5.  The FDT
    # itself is proven at ALL temperatures by layers 1–2 (exact AutoDiff); this
    # layer ties the registered (N→∞) closed forms to a finite-N enumeration,
    # which is only tight where finite-size corrections are negligible.
    βs = [0.2, 0.3, 0.4]

    # (a) Direct, transparent: fetched closed form vs brute-force fluctuation.
    for β in βs
        c_v = fetch(m, SpecificHeat(), Infinite(); beta=β)
        χ = fetch(m, SusceptibilityZZ(), Infinite(); beta=β)
        @test isapprox(
            c_v,
            β^2 * independent_energy_variance_per_site(m, Infinite(); beta=β);
            rtol=1e-4,
        )
        @test isapprox(
            χ,
            β * independent_magnetization_variance_per_site(m, Infinite(); beta=β);
            rtol=1e-4,
        )
    end

    # (b) Same statements through the ThermoIdentity harness (exercises the
    # FLUCTUATION_DISSIPATION_IDENTITIES integration).  All must RUN (not skip)
    # and pass — a skip would mean the trait/dispatch wiring broke.
    results = verify_thermodynamic_identities(
        m,
        Infinite();
        βs=βs,
        identities=FLUCTUATION_DISSIPATION_IDENTITIES,
        rtol=1e-4,
        atol=1e-8,
    )
    @test length(results) == length(FLUCTUATION_DISSIPATION_IDENTITIES) * length(βs)
    @test all(r -> r.status === :pass, results)
    ran = unique(r.identity for r in results)
    @test any(x -> occursin("Var(E)", x), ran)   # energy FDT ran
    @test any(x -> occursin("Var(M)", x), ran)   # magnetisation FDT ran
end

# ── BONUS — the existing thermodynamic identities also hold for IsingChain1D ──
# (No IsingChain1D identity file existed; this fills that gap and confirms the
# new util integrates with the standard DEFAULT_IDENTITIES path.)
@testset "IsingChain1D — standard thermodynamic identities (Gibbs, c_v) via harness" begin
    m = IsingChain1D(; J=1.0)
    results = verify_thermodynamic_identities(m, Infinite(); βs=[0.3, 0.7, 1.2])
    @test all(r -> r.status !== :fail, results)
    @test any(r -> occursin("Gibbs", r.identity) && r.status === :pass, results)
    @test any(r -> occursin("∂ε/∂β", r.identity) && r.status === :pass, results)
    @test any(r -> occursin("∂s/∂β", r.identity) && r.status === :pass, results)
end
