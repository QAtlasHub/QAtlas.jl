# Constraint-edge layer (#697): kernel + @symmetry/@identity/@dual/@limits_to.
#
# Unit tests of the SHARED kernel (EDGE_STORES registration, GeneratedCheck
# protocol, hub enumeration) and of each edge type's store, macro, queries and
# coherence checks — including synthetic violations pushed and popped in the
# test_reduces.jl probe-row style.  The PHYSICS of the generated checks runs
# in test/generated/; here we test the machinery.

using QAtlas, Test
using QAtlas:
    TFIM,
    XXZ1D,
    Heisenberg1D,
    S1Heisenberg1D,
    Kitaev1D,
    Infinite,
    OBC,
    PBC,
    Energy,
    FreeEnergy,
    ThermalEntropy,
    SpecificHeat,
    MassGap,
    ChargeGap,
    SusceptibilityXX,
    SusceptibilityYY,
    SusceptibilityZZ,
    MagnetizationX,
    MagnetizationZLocal,
    XXStructureFactor,
    VonNeumannEntropy,
    LuttingerVelocity,
    component,
    symmetry_profile,
    models_with_symmetry,
    identities_for,
    participants,
    dualities,
    limits_from,
    limits_into,
    generated_checks,
    run_generated_check

# Probe model types for synthetic store entries (structs must be top-level).
struct _LSMProbe <: QAtlas.AbstractQAtlasModel end
struct _DegProbe <: QAtlas.AbstractQAtlasModel end

# ── quantity taxonomy (#690 slice) ───────────────────────────────────
@testset "quantity family supertypes are purely additive" begin
    @test FreeEnergy <: AbstractThermalPotential <: QAtlas.AbstractQuantity
    @test Energy{:per_site} <: AbstractThermalPotential
    @test SusceptibilityXX <: AbstractSusceptibility
    @test MagnetizationX <: AbstractMagnetization
    @test MagnetizationZLocal <: AbstractMagnetization
    @test XXStructureFactor <: AbstractStructureFactor
    @test MassGap <: AbstractGap && ChargeGap <: AbstractGap
    @test LuttingerVelocity <: AbstractVelocity
    @test VonNeumannEntropy{:equilibrium} <: AbstractEntanglementMeasure
end

@testset "component trait recovers the name-encoded index" begin
    @test component(SusceptibilityXX) === :xx
    @test component(SusceptibilityYY()) === :yy
    @test component(MagnetizationX) === :x
    @test component(MassGap) === :mass
    @test component(ChargeGap) === :charge
    # …Local variants deliberately carry NO component (different fetch shape)
    @test component(MagnetizationZLocal) === nothing
    # default: quantities without a family index
    @test component(FreeEnergy) === nothing
end

# ── kernel: store registration ───────────────────────────────────────
@testset "EDGE_STORES registration covers every declarative store" begin
    names = Set(s.name for s in QAtlas.EDGE_STORES)
    for expected in
        (:registry, :realizes, :reduces, :about, :symmetry, :identity, :dual, :limits_to)
        @test expected in names
    end
    # accessors work on a real row of every non-empty store
    for spec in QAtlas.EDGE_STORES
        isempty(spec.store) && continue
        row = first(spec.store)
        @test spec.references_of(row) isa Vector{String}
        @test spec.location_of(row) isa String && !isempty(spec.location_of(row))
    end
    # duplicate registration is rejected
    @test_throws ArgumentError QAtlas.register_edge_store!(:registry, QAtlas.REGISTRY)
end

# ── kernel: hub enumeration ──────────────────────────────────────────
@testset "_implemented_hubs finds models covering a quantity set" begin
    hubs = QAtlas._implemented_hubs((Energy{:per_site}, FreeEnergy, ThermalEntropy))
    @test (model=TFIM, bc=Infinite) in hubs
    # deterministic: sorted by (model, bc) names
    @test hubs == sort(hubs; by=h -> (string(h.model), string(h.bc)))
    # Universality / Bound namespaces never appear as hubs
    @test all(!(h.model <: QAtlas.Universality) && !(h.model <: QAtlas.Bound) for h in hubs)
    @test isempty(QAtlas._implemented_hubs(()))
end

# ── @symmetry ────────────────────────────────────────────────────────
@testset "@symmetry profiles + queries" begin
    p = symmetry_profile(Heisenberg1D)
    @test p !== nothing
    @test p.internal === :SU2 && p.translation && p.site_spin == 1//2
    @test p.gapped === false
    @test symmetry_profile(Heisenberg1D()) === p          # instance form
    @test Heisenberg1D in models_with_symmetry(:SU2)
    @test XXZ1D in models_with_symmetry(:U1)              # generic-parameter group
    @test !(XXZ1D in models_with_symmetry(:SU2))          # Δ=1 enhancement is NOT generic
    @test symmetry_profile(Kitaev1D).internal === :Z2     # fermion parity
    @test symmetry_profile(Kitaev1D).site_spin === nothing # spinless fermions
    @test symmetry_profile(_DegProbe) === nothing          # absence is legal, not an error
    # validation: one profile per family; gs_degeneracy needs gapped=true
    @test_throws ArgumentError QAtlas.symmetry!(Heisenberg1D; internal=:SU2)
    @test_throws ArgumentError QAtlas.symmetry!(_DegProbe; internal=:U1, gs_degeneracy=2)
end

@testset "C10 LSM consistency fires on a violating declaration" begin
    # shipped profiles must be LSM-consistent
    @test isempty(filter(f -> f.severity === :error, QAtlas.check_lsm_consistency()))

    # gapped + unique GS + SU(2) + translation + spin-1/2 ⟹ contradicts LSM.
    # try/finally so a failing assertion cannot leak the probe into the global
    # store and corrupt every later test file.
    try
        QAtlas.symmetry!(
            _LSMProbe;
            internal=:SU2,
            translation=true,
            site_spin=1//2,
            gapped=true,
            gs_degeneracy=1,
        )
        fs = QAtlas.check_lsm_consistency()
        @test any(f -> f.severity === :error && occursin("_LSMProbe", f.message), fs)
    finally
        pop!(QAtlas.SYMMETRY_PROFILES)
    end

    # gapped but degeneracy undeclared ⟹ self-reported :gap, not :error
    try
        QAtlas.symmetry!(
            _LSMProbe; internal=:U1, translation=true, site_spin=3//2, gapped=true
        )
        fs = QAtlas.check_lsm_consistency()
        @test any(f -> f.severity === :gap && occursin("_LSMProbe", f.message), fs)
        @test !any(f -> f.severity === :error && occursin("_LSMProbe", f.message), fs)
    finally
        pop!(QAtlas.SYMMETRY_PROFILES)
    end

    # integer spin (Haldane): gapped + unique is fine — S1Heisenberg1D ships so
    @test !any(f -> occursin("S1Heisenberg1D", f.message), QAtlas.check_lsm_consistency())
end

@testset "C10b symmetry corroboration runner detects a contradiction" begin
    # every shipped corroboration check passes (gapped facts match MassGap)
    shipped = generated_checks(; kinds=(:symmetry,))
    @test !isempty(shipped)
    for c in shipped
        @test run_generated_check(c).status === :pass
    end

    # Flip Heisenberg1D's profile to a CONTRADICTING gapped=true while its
    # registered MassGap@Infinite is exactly 0 (gapless) — the runner must
    # return :fail.  Swap-and-restore under try/finally so the global store is
    # untouched after the test, regardless of assertion outcome.
    idx = findfirst(p -> p.model === Heisenberg1D, QAtlas.SYMMETRY_PROFILES)
    @test idx !== nothing
    saved = QAtlas.SYMMETRY_PROFILES[idx]
    try
        deleteat!(QAtlas.SYMMETRY_PROFILES, idx)
        QAtlas.symmetry!(
            Heisenberg1D;
            internal=:SU2,
            translation=true,
            site_spin=1//2,
            gapped=true,
            gs_degeneracy=2,   # keeps C10 LSM satisfied (degenerate, not unique)
        )
        chk = only(
            filter(
                c -> startswith(c.id, "symmetry/gapped/Heisenberg1D/"),
                generated_checks(; kinds=(:symmetry,)),
            ),
        )
        out = run_generated_check(chk)
        @test out.status === :fail   # declares gapped=true but the fetched gap is 0
    finally
        i2 = findfirst(p -> p.model === Heisenberg1D, QAtlas.SYMMETRY_PROFILES)
        i2 === nothing || deleteat!(QAtlas.SYMMETRY_PROFILES, i2)
        insert!(QAtlas.SYMMETRY_PROFILES, idx, saved)
    end
end

# ── @identity ────────────────────────────────────────────────────────
@testset "@identity store + queries" begin
    edges = identities_for(FreeEnergy)
    @test any(e -> e.name === :gibbs, edges)
    gibbs = only(filter(e -> e.name === :gibbs, QAtlas.IDENTITIES))
    @test gibbs isa QAtlas.TupleIdentityEdge
    @test Set(participants(gibbs)) == Set([Energy{:per_site}, FreeEnergy, ThermalEntropy])
    # family identities pick up members through the taxonomy + component
    iso = only(filter(e -> e.name === :su2_susceptibility_isotropy, QAtlas.IDENTITIES))
    @test iso isa QAtlas.IsotropyIdentityEdge
    @test SusceptibilityXX in participants(iso)
    @test iso in identities_for(SusceptibilityZZ)
    # the two concrete edge types share AbstractIdentityEdge
    @test gibbs isa QAtlas.AbstractIdentityEdge && iso isa QAtlas.AbstractIdentityEdge
    # validation
    @test_throws ArgumentError QAtlas.identity!(:gibbs; family=AbstractSusceptibility)
    @test_throws ArgumentError QAtlas.identity!(:_neither)
    @test_throws ArgumentError QAtlas.identity!(
        :_both;
        quantities=(a=FreeEnergy, b=ThermalEntropy),
        check=(v, p) -> (v.a, v.b),
        family=AbstractSusceptibility,
    )
    @test_throws ArgumentError QAtlas.identity!(:_concrete_family; family=FreeEnergy)
    # tolerance footgun: rtol ≥ 1 passes everything → rejected
    @test_throws ArgumentError QAtlas.identity!(
        :_bad_rtol;
        quantities=(a=FreeEnergy, b=ThermalEntropy),
        check=(v, p) -> (v.a, v.b),
        rtol=1.0,
    )
end

# ── @dual ────────────────────────────────────────────────────────────
@testset "@dual store + queries + C12" begin
    ds = dualities(TFIM)
    @test any(d -> d.kind === :kramers_wannier, ds)
    @test any(d -> d.kind === :jordan_wigner && d.target === Kitaev1D, ds)
    @test any(d -> d.kind === :jordan_wigner, dualities(Kitaev1D))  # backlink
    # shipped edges: param_map sanity holds (zero :error)
    @test isempty(filter(f -> f.severity === :error, QAtlas.check_duality_maps()))

    # synthetic: param_map landing on the wrong type is an :error, and the
    # malformed image must NOT crash the check (the C12 `continue` guard) — so
    # the whole report still returns rather than throwing.
    try
        QAtlas.dual!(
            :_broken_probe,
            TFIM,
            Kitaev1D;
            param_map=m -> m,                   # returns a TFIM, not a Kitaev1D
            kind=:probe,
            quantities=[(quantity=MassGap, bc=Infinite)],
            examples=[TFIM(; J=1.0, h=0.5)],
        )
        fs = QAtlas.check_duality_maps()        # must not throw
        @test any(f -> f.severity === :error && occursin("_broken_probe", f.message), fs)
    finally
        pop!(QAtlas.DUALITIES)
    end

    # validation: empty examples / empty quantities rejected
    @test_throws ArgumentError QAtlas.dual!(
        :_no_examples,
        TFIM,
        TFIM;
        param_map=identity,
        kind=:probe,
        quantities=[(quantity=MassGap, bc=Infinite)],
        examples=[],
    )
    # an involution maps a manifold to itself: distinct endpoints rejected
    @test_throws ArgumentError QAtlas.dual!(
        :_inv_mismatch,
        TFIM,
        Kitaev1D;
        param_map=identity,
        kind=:probe,
        involution=true,
        quantities=[(quantity=MassGap, bc=Infinite)],
        examples=[TFIM(; J=1.0, h=0.5)],
    )
end

# ── @limits_to ───────────────────────────────────────────────────────
@testset "@limits_to store + queries + C13" begin
    ls = limits_from(XXZ1D)
    @test any(l -> l.target === Heisenberg1D && l.param === :Δ, ls)
    @test any(l -> l.source === XXZ1D, limits_into(Heisenberg1D))   # backlink
    @test isempty(filter(f -> f.severity === :error, QAtlas.check_limit_edges()))
    # validation: a non-monotone approach sequence is rejected
    @test_throws ArgumentError QAtlas.limits_to!(
        :_bad_approach,
        XXZ1D,
        Heisenberg1D;
        param=:Δ,
        approach=[1.1, 1.3, 1.01],
        regime="probe",
        quantities=[(quantity=Energy{:per_site}, bc=Infinite, final_atol=1e-3)],
    )
    @test_throws ArgumentError QAtlas.limits_to!(
        :_short_approach,
        XXZ1D,
        Heisenberg1D;
        param=:Δ,
        approach=[1.1],
        regime="probe",
        quantities=[(quantity=Energy{:per_site}, bc=Infinite, final_atol=1e-3)],
    )
end

# ── kernel: generated-check protocol ─────────────────────────────────
@testset "generated_checks aggregates deterministically" begin
    checks = generated_checks()
    @test !isempty(checks)
    @test issorted([c.id for c in checks])
    @test allunique(c.id for c in checks)
    @test all(c.kind in (:identity, :dual, :limit, :symmetry) for c in checks)
    # kind selection
    only_ids = generated_checks(; kinds=(:identity,))
    @test all(c -> c.kind === :identity, only_ids)
    @test !isempty(only_ids)
    # the same registry state generates the same plan (determinism)
    @test [c.id for c in generated_checks()] == [c.id for c in checks]
end

@testset "a generated check runs to a CheckOutcome" begin
    checks = generated_checks(; kinds=(:identity,))
    i = findfirst(c -> startswith(c.id, "identity/gibbs/TFIM/Infinite"), checks)
    @test i !== nothing
    out = run_generated_check(checks[i])
    @test out isa QAtlas.CheckOutcome
    @test out.status === :pass
    @test out.abs_err < 1e-8
    # a throwing runner is reported as :error (NOT :fail — a thrown exception is
    # a config/dispatch bug, not a numerical disagreement), never raises
    boom = QAtlas.GeneratedCheck(:identity, "probe/boom", "probe", () -> error("boom"))
    out = run_generated_check(boom)
    @test out.status === :error && occursin("boom", out.detail)
end

# ── the hard invariant stays: zero coherence errors with C10–C13 wired ──
@testset "coherence_report stays error-free with the constraint layer" begin
    fs = QAtlas.coherence_report()
    errs = QAtlas.coherence_errors(fs)
    isempty(errs) || foreach(println, errs)
    @test isempty(errs)
end
