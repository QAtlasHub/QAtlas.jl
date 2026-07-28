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
    # A MODEL axis is now supported (`_diff_target` rebuilds the model with
    # `_with_param`), so `:h` is no longer rejected — it is pinned to finite
    # differences instead.  What must still be rejected is a field the model
    # does not have; that throws at evaluation, from `_with_param`.
    @test QAtlas.∂(QAtlas.FreeEnergy, :h) isa QAtlas.DerivedInput
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

@testset "model-axis derivatives are opt-in, placed, and finite-difference only" begin
    # The allow-list is the defence the cross-check cannot provide: for a
    # transverse-field model both backends agree on -∂F/∂h = ⟨σˣ⟩, which is not
    # M_z, so nothing numerical would flag it.
    ids = [c.id for c in generated_checks(; kinds=(:response,))]
    mag = filter(startswith("response/magnetization_response/"), ids)
    sus = filter(startswith("response/susceptibility_response/"), ids)
    # Log what was generated: when one of these assertions fails the useful
    # information is WHICH hubs appeared, and `@test any(...)` does not print it.
    @info "model-axis edges" mag sus
    @test !isempty(mag)
    @test !isempty(sus)
    # The allow-list is opt-IN, so only the two longitudinal-field models appear.
    @test all(i -> occursin("/CurieWeissIsing/", i) || occursin("/IsingChain1D/", i), mag)
    @test all(i -> occursin("/CurieWeissIsing/", i), sus)
    # A transverse-field model must never be checked against M_z = -∂F/∂h.
    @test !any(occursin("/TFIM/"), mag)
    # IsingChain1D's SusceptibilityZZ is h = 0 only, and this edge sits at h ≠ 0.
    @test !any(occursin("/IsingChain1D/"), sus)
    # An edge with no `models` allow-list reaches strictly more hubs than one
    # with it.  (`entropy_response` happens to exclude TFIM for an unrelated
    # reason — see _THERMO_DERIVATIVE_EXCLUSIONS — so assert the structural
    # property, not a particular model.)
    open_edge = filter(startswith("response/entropy_response/"), ids)
    hubs_of(v) = Set(join(split(i, "/")[3:4], "/") for i in v)
    @test length(hubs_of(open_edge)) > length(hubs_of(mag))

    # A model axis is pinned to finite differences: rebuilding a struct whose
    # fields are ::Float64 with an AD dual destroys the derivative silently.
    fd_axis = QAtlas.∂(QAtlas.FreeEnergy, :h)
    st_axis = QAtlas.∂(QAtlas.FreeEnergy, :T)
    @test QAtlas._axis_backend(fd_axis, QAtlas.ForwardDiffBackend()) isa
        QAtlas.FiniteDifference
    @test QAtlas._axis_backend(st_axis, QAtlas.ForwardDiffBackend()) isa
        QAtlas.ForwardDiffBackend
end

@testset "a placed edge is evaluated where the identity has content" begin
    # `at` exists because the DEFAULT model is the useless point for a field
    # response: at h = 0 both M_z and -∂F/∂h vanish by symmetry, so the check
    # passes while testing nothing — and below T_c it is worse than useless,
    # because F has a kink there and a central difference straddles the jump.
    # Measured on CurieWeissIsing at βJ = 2: m(0⁺) = 0.9575, central diff = 0.
    m0 = QAtlas.CurieWeissIsing(; J=1.0, h=0.0)
    placed = QAtlas._at_params(m0, (h=0.35,))
    @test placed.h == 0.35
    @test placed.J == m0.J
    @test QAtlas._at_params(m0, NamedTuple()) === m0
    # A field the model does not have is a declaration bug, not a silent no-op.
    @test_throws ArgumentError QAtlas._at_params(m0, (nosuchfield=1.0,))

    # The vacuity this guards against, stated as an assertion.
    @test QAtlas.fetch(m0, QAtlas.Magnetization{:z}(), QAtlas.Infinite(); beta=0.5) == 0.0
    @test QAtlas.fetch(placed, QAtlas.Magnetization{:z}(), QAtlas.Infinite(); beta=0.5) >
        0.1
end
