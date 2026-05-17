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
