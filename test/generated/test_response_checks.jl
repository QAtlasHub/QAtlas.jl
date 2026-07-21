# Generated RESPONSE checks (#734 Phase B) — the derivative-supplied slice.
#
# These are the AbstractQAtlas relations that were reachable on QAtlas hubs but
# unusable because one slot is a derivative rather than a fetched value.
#
# Note the backend: `runtests.jl` loads ForwardDiff, so the extension is active
# here and these run on AD, at AD tolerance.  Without it the suite still works,
# on finite differences at `default_rtol(FiniteDifference())` — that fallback is
# covered in test/core/test_derivative.jl.

include("util_run_checks.jl")
using QAtlas: generated_checks, RESPONSES, preferred_backend, ForwardDiffBackend

@testset "generated response checks" begin
    # The test environment loads ForwardDiff, so the AD extension must be the
    # one that ran — otherwise these would silently be FD results reported at
    # an AD tolerance.
    @test preferred_backend() isa ForwardDiffBackend

    checks = generated_checks(; kinds=(:response,))
    @test !isempty(checks)
    ids = [c.id for c in checks]
    @test any(startswith("response/entropy_response/"), ids)
    @test any(startswith("response/specific_heat_from_entropy/"), ids)
    @test any(startswith("response/gibbs_helmholtz/"), ids)
    @test any(startswith("response/specific_heat_fdt/"), ids)
    @test length(unique(ids)) == length(ids)

    run_generated_suite(checks; label="generated response checks")
end

@testset "response edges resolve their slots from the relation" begin
    @test !isempty(RESPONSES)
    for e in RESPONSES
        # The subject is derived, never declared; it must be a real slot name.
        @test e.subject in first.(QAtlas.variable_slots(e.relation))
        # Every untyped slot is supplied — an unsupplied one could not run.
        untyped = [n for (n, T) in QAtlas.variable_slots(e.relation) if T === nothing]
        @test Set(untyped) == Set(keys(e.derived))
    end
end

@testset "response! refuses what it cannot materialize" begin
    # A field derivative needs a model-parameter mechanism that does not exist.
    @test_throws ArgumentError QAtlas.∂(QAtlas.FreeEnergy, :h)
    # An untyped slot left unsupplied would be a check that never runs.
    @test_throws ArgumentError QAtlas.response!(
        :_probe_unsupplied; relation=QAtlas.EntropyResponse, derived=NamedTuple()
    )
    @test_throws ArgumentError QAtlas.response!(
        :entropy_response;
        relation=QAtlas.EntropyResponse,
        derived=(dF_dT=QAtlas.∂(QAtlas.FreeEnergy, :T),),
    )
end

@testset "derived-input transforms do what they say" begin
    # `of` picks WHAT is differentiated: d(βF)/dβ, not dF/dβ.
    dβF = QAtlas.∂(QAtlas.FreeEnergy, :β; of=(F, β) -> β * F)
    @test dβF.of(2.0, 3.0) == 6.0
    @test dβF.then(5.0) == 5.0
    # `then` post-processes: Var(E) = -∂⟨E⟩/∂β.
    var_E = QAtlas.∂(QAtlas.Energy{:per_site}, :β; then=d -> -d)
    @test var_E.then(5.0) == -5.0
    @test var_E.of(2.0, 3.0) == 2.0
    # Defaults stay the identity, so an untransformed edge is unaffected.
    plain = QAtlas.∂(QAtlas.FreeEnergy, :T)
    @test plain.of(2.0, 3.0) == 2.0 && plain.then(5.0) == 5.0
end

@testset "the two specific-heat routes are independent" begin
    # :specific_heat_from_entropy derives C from S, :specific_heat_fdt from U.
    # A model computing C by one formula and S or U by another disagrees with
    # exactly one of them — which is the point of running both.
    ids = [c.id for c in generated_checks(; kinds=(:response,))]
    hubs(pre) = Set(join(split(i, "/")[3:4], "/") for i in ids if startswith(i, pre))
    from_S = hubs("response/specific_heat_from_entropy/")
    from_U = hubs("response/specific_heat_fdt/")
    @test !isempty(from_S) && !isempty(from_U)
    @test !isempty(intersect(from_S, from_U))   # hubs covered by both routes
end

@testset "model-axis derivatives are opt-in and finite-difference only" begin
    # The allow-list is the defence the cross-check cannot provide: for a
    # transverse-field model both backends agree on -∂F/∂h = ⟨σˣ⟩, which is not
    # M_z, so nothing numerical would flag it.  TFIM must therefore be absent.
    ids = [c.id for c in generated_checks(; kinds=(:response,))]
    mag = filter(startswith("response/magnetization_response/"), ids)
    @test !isempty(mag)
    @test any(occursin("/CurieWeissIsing/"), mag)
    @test !any(occursin("/TFIM/"), mag)          # transverse field — must not be checked
    # An edge with no allow-list keeps generating everywhere.
    @test any(occursin("/TFIM/"), filter(startswith("response/entropy_response/"), ids))

    # A model axis is pinned to finite differences: rebuilding a struct whose
    # fields are ::Float64 with an AD dual destroys the derivative silently.
    fd_axis = QAtlas.∂(QAtlas.FreeEnergy, :h)
    st_axis = QAtlas.∂(QAtlas.FreeEnergy, :T)
    @test QAtlas._axis_backend(fd_axis, QAtlas.ForwardDiffBackend()) isa
        QAtlas.FiniteDifference
    @test QAtlas._axis_backend(st_axis, QAtlas.ForwardDiffBackend()) isa
        QAtlas.ForwardDiffBackend
end
