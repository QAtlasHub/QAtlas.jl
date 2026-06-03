# Self-test for the black-box verification harness (test/util/verify.jl
# + test/util/generic_ed.jl).  Proves: (1) the subject is fetched, not
# re-typed; (2) a model-agnostic ED route (built from the *physics*, not
# any QAtlas internal builder) reproduces a src closed form; (3) route
# validation rejects anything outside the independent vocabulary.

using QAtlas, Test
using LinearAlgebra: kron

@testset "verify harness — black-box independence (self-test)" begin
    # AKLT bond from the WRITTEN Hamiltonian H = Σ S·S + 1/3 (S·S)²,
    # spin-1, J = 1 — reconstructed purely from generic_ed primitives.
    Sx, Sy, Sz = spin_ops(1)
    SS = kron(Sx, Sx) + kron(Sy, Sy) + kron(Sz, Sz)
    bond = SS + (1 / 3) * (SS * SS)

    Ns = (4, 6)
    ed_eps = Float64[]
    for N in Ns
        H = chain_hamiltonian(3, N, bond)          # OBC, no src builder
        E0, _ = ground_state(H)
        push!(ed_eps, E0 / (N - 1))                # frustration-free ⇒ -2/3
    end

    # `verify` fetches the subject itself (AKLT1D Energy{:per_site},
    # Infinite = closed form -2J/3); we only hand it the independent ED.
    s = verify(
        AKLT1D(),
        Energy(:per_site),
        Infinite();
        route=:ed_finite_size,
        independent=ed_eps,
        at=["N=$N" for N in Ns],
        agree_within=1e-9,
        refs=["Affleck-Kennedy-Lieb-Tasaki 1988"],
    )
    @test s ≈ -2 / 3 atol = 1e-12          # subject came from fetch
    @test all(e -> isapprox(e, -2 / 3; atol=1e-10), ed_eps)  # ED reproduced it

    # Route vocabulary is closed — "retype the formula" is unexpressible.
    @test_throws ErrorException verify(
        AKLT1D(),
        Energy(:per_site),
        Infinite();
        route=:retype_formula,
        independent=[-2 / 3],
        agree_within=1e-9,
        refs=["x"],
    )
end

@testset "verify_bound harness — variational one-sided witness" begin
    # Rayleigh–Ritz: the |→…→⟩ product-state energy density (-h exactly,
    # since ⟨X⟩=1 and ⟨ZZ⟩=0 there) is an UPPER bound for the TFIM
    # ground-state energy density. So the fetched GS (subject) sits at or
    # below -h, i.e. the measured witness -h is ≥ subject. This is a real
    # one-sided check, not an equality dressed up as a bound.
    J, h = 1.0, 2.0
    m = TFIM(; J=J, h=h)
    trial_density = -h                       # ⟨H⟩/N for |→…→⟩, analytic
    s = verify_bound(
        m,
        Energy(:per_site),
        Infinite();
        route=:variational_state,
        measured=[trial_density],
        relation=:geq,
        refs=["Rayleigh–Ritz variational principle"],
    )
    @test s == QAtlas.fetch(m, Energy(:per_site), Infinite())   # subject came from fetch
    @test trial_density >= s                              # the bound genuinely holds

    # Both vocabularies are closed — a typo can't silently weaken the card.
    @test_throws ErrorException verify_bound(
        m,
        Energy(:per_site),
        Infinite();
        route=:retype,
        measured=[trial_density],
        relation=:geq,
        refs=["x"],
    )
    @test_throws ErrorException verify_bound(
        m,
        Energy(:per_site),
        Infinite();
        route=:variational_state,
        measured=[trial_density],
        relation=:above,
        refs=["x"],
    )
end

@testset "verify_approx harness — domain-limited high-T tail" begin
    # TFIM specific heat decays as ~(βΔ)² at high T, so c → 0 as β → 0.
    # In-domain (βJ ≪ 1) the fetched c agrees with the leading value 0 to
    # the stated tolerance; the domain + error order ride along on the card.
    m = TFIM(; J=1.0, h=1.0)
    s = verify_approx(
        m,
        SpecificHeat(),
        Infinite();
        route=:high_temperature,
        reference=0.0,
        agree_within=1e-3,
        valid_domain="betaJ << 1",
        error_order="O((betaJ)^2)",
        refs=["high-T tail: c ~ (beta*Delta)^2"],
        fetch_kw=(; beta=1e-3),
    )
    @test isapprox(s, 0.0; atol=1e-3)        # in-domain: matches leading value

    @test_throws ErrorException verify_approx(
        m,
        SpecificHeat(),
        Infinite();
        route=:bogus,
        reference=0.0,
        agree_within=1e-3,
        valid_domain="x",
        error_order="y",
        refs=["z"],
        fetch_kw=(; beta=1e-3),
    )
end
