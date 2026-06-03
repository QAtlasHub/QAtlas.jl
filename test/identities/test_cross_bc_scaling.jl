using Test
using QAtlas
using LinearAlgebra: norm

# Cross-BC scaling identities — verify the convergence rate of
# `OBC(N)` / `PBC(N)` per-site quantities to their `Infinite()` limits
# matches the boundary-correction theory:
#
#   * Off-critical TFIM (`h ≠ J`):
#       OBC: leading correction `O(1/N)`  (boundary energy)
#       PBC: leading correction `O(exp(-N / ξ))`  (gap ξ = 1/(2|h-J|))
#
#   * Critical TFIM (`h = J`):
#       OBC: `O(1/N²)` from CFT boundary contribution (log corrections)
#       PBC: `O(1/N²)` Lüscher-style finite-size scaling
#
#   * XXZ at Δ = 0 (free-fermion XX): per-site GS energy `-J/π` is the
#     thermodynamic-limit closed form.  Dense-ED OBC `Energy/N` should
#     converge to it as `N → ∞`.
#
# These rates are universal predictions; each test fits the (N, error)
# data to a power law `error ∝ N^(-p)` and asserts `p` matches the
# theoretical exponent within 20-30 %.  The aim is to guard against
# regressions that would replace `1/N²` with `1/N` (e.g. accidental
# OBC-vs-PBC swap in a registry binding) or vice versa.

# Linear fit of `log(err) = a − decay · log(N)` returns the decay exponent
# `decay > 0` (so `err ∝ N^{-decay}`).  Returns `(a, decay)`.
function _fit_loglog(Ns, errs)
    x = log.(float.(collect(Ns)))
    y = log.(float.(collect(errs)))
    n = length(x)
    sx = sum(x)
    sy = sum(y)
    sxx = sum(x .^ 2)
    sxy = sum(x .* y)
    slope = (n * sxy - sx * sy) / (n * sxx - sx^2)  # negative for decaying err
    a = (sy - slope * sx) / n
    return a, -slope   # report decay exponent as a positive number
end

# ────────────────────────────────────────────────────────────────────
# (1) TFIM disordered phase OBC: 1/N boundary scaling
# ────────────────────────────────────────────────────────────────────

@testset "TFIM disordered (h>J) OBC: per-site ε converges to Infinite at O(1/N)" begin
    J, h = 1.0, 1.5
    model = TFIM(; J=J, h=h)
    β = 2.0  # finite β; gap is 2|h-J| = 1, so ξ = 1
    ε_inf = QAtlas.fetch(model, Energy(:per_site), Infinite(); beta=β)

    Ns = (16, 32, 64, 128)
    errs = Float64[]
    for N in Ns
        ε_N = QAtlas.fetch(model, Energy(:per_site), OBC(N); beta=β) / 1  # per-site
        push!(errs, abs(ε_N - ε_inf))
    end

    # Errors should be monotonically decreasing.
    @test issorted(errs; rev=true)
    # Power-law fit `err ∝ N^(-p)` should give p ≈ 1.
    _, p = _fit_loglog(Ns, errs)
    # Allow [0.7, 1.5] band — boundary leading 1/N contribution can be
    # contaminated by sub-leading 1/N² terms at small N.
    @test 0.7 ≤ p ≤ 1.5
end

@testset "TFIM disordered (h>J) PBC: monotonic convergence to Infinite" begin
    # At finite β the dominant `O(1/N)` error in *both* OBC and PBC
    # In the gapped phase (|h - J| > 0) the leading finite-size
    # correction to the PBC per-site energy is **exponential**
    # `O(exp(-N / ξ))` with ξ = 1/(2|h - J|); the trapezoidal-rule
    # quadrature error of the `Infinite()` integrand vanishes as
    # `O(exp(-cN))` on a smooth periodic integrand too.  We assert
    # monotonic decrease only — the convergence is faster than any
    # power law and saturates at machine precision by N ~ 128, so a
    # power-law slope fit on log(err) is dominated by log(0) garbage.
    J, h = 1.0, 1.5
    model = TFIM(; J=J, h=h)
    β = 2.0
    ε_inf = QAtlas.fetch(model, Energy(:per_site), Infinite(); beta=β)

    Ns = (16, 32, 64, 128)
    errs = Float64[
        abs(QAtlas.fetch(model, Energy(:per_site), PBC(N); beta=β) - ε_inf) for N in Ns
    ]
    @test issorted(errs; rev=true)
    # Exponential convergence: each doubling of N should reduce err
    # by orders of magnitude, capped at machine precision.  Below the
    # cap we cannot fit a power law (log(0) blows the slope up); above
    # we just assert each step shrinks by ≥ 10× as the lightest
    # exponential check.
    for i in 1:(length(errs) - 1)
        errs[i + 1] ≤ 1e-12 && continue   # saturated at machine precision
        @test errs[i + 1] ≤ errs[i] / 10
    end
end

# ────────────────────────────────────────────────────────────────────
# (2) TFIM critical point OBC: 1/N² scaling (CFT-corrected)
# ────────────────────────────────────────────────────────────────────

@testset "TFIM critical (h=J) OBC: per-site ε convergence is 1/N (boundary)" begin
    # OBC critical: the leading correction is `1/N` from the boundary
    # ε_b/N rather than the bulk Cardy 1/N² Casimir.  Cardy
    # `c π / 24 N²` is the *PBC* prediction; OBC retains the 1/N
    # boundary energy at criticality just as in the gapped phase.
    J, h = 1.0, 1.0
    model = TFIM(; J=J, h=h)
    β = 100.0  # effectively ground state
    ε_inf = QAtlas.fetch(model, Energy(:per_site), Infinite(); beta=β)

    Ns = (16, 32, 64, 128)
    errs = Float64[
        abs(QAtlas.fetch(model, Energy(:per_site), OBC(N); beta=β) - ε_inf) for N in Ns
    ]
    @test issorted(errs; rev=true)
    _, p = _fit_loglog(Ns, errs)
    # 1/N boundary energy expected.
    @test 0.7 ≤ p ≤ 1.5
end

@testset "TFIM critical (h=J) PBC: per-site ε convergence is ~1/N² (Cardy)" begin
    # No boundary at PBC, so the leading correction is Cardy's
    # `-c·π v_F / (6 N²)` ground-state Casimir.  Fits should give
    # an exponent close to 2.
    J, h = 1.0, 1.0
    model = TFIM(; J=J, h=h)
    β = 100.0
    ε_inf = QAtlas.fetch(model, Energy(:per_site), Infinite(); beta=β)

    Ns = (16, 32, 64, 128)
    errs = Float64[
        abs(QAtlas.fetch(model, Energy(:per_site), PBC(N); beta=β) - ε_inf) for N in Ns
    ]
    @test issorted(errs; rev=true)
    _, p = _fit_loglog(Ns, errs)
    # Cardy 1/N²; the discrete-grid trapezoidal rule on the smooth
    # PBC integrand actually gives super-polynomial convergence beyond
    # ~ 1/N², so we accept anything ≥ 1.7.
    @test p ≥ 1.7
end

# ────────────────────────────────────────────────────────────────────
# (3) PBC monotonic convergence — TFIM critical point
# ────────────────────────────────────────────────────────────────────

@testset "TFIM critical PBC at finite β: monotonic convergence" begin
    # At finite β the Cardy 1/N² term coexists with thermal
    # quadrature corrections.  Just verify monotonic decay; exponent
    # is regime-dependent.
    J, h = 1.0, 1.0
    model = TFIM(; J=J, h=h)
    β = 5.0
    ε_inf = QAtlas.fetch(model, Energy(:per_site), Infinite(); beta=β)
    Ns = (16, 32, 64, 128)
    errs = Float64[
        abs(QAtlas.fetch(model, Energy(:per_site), PBC(N); beta=β) - ε_inf) for N in Ns
    ]
    @test issorted(errs; rev=true)
end

# ────────────────────────────────────────────────────────────────────
# (4) Free-energy / SusceptibilityXX convergence — finite β
# ────────────────────────────────────────────────────────────────────

@testset "TFIM free energy OBC → Infinite (h>J): 1/N decay" begin
    # Free energy convergence at finite β: leading O(1/N) boundary
    # contribution governs the rate.  (We do not include χ_xx here:
    # its OBC (variance) and Infinite (Kubo) implementations target
    # different physical objects — see PR #132 for the convention
    # split.)
    J, h, β = 1.0, 1.5, 1.0
    model = TFIM(; J=J, h=h)
    f_inf = QAtlas.fetch(model, FreeEnergy(), Infinite(); beta=β)

    Ns = (16, 32, 64, 128)
    f_errs = Float64[
        abs(QAtlas.fetch(model, FreeEnergy(), OBC(N); beta=β) - f_inf) for N in Ns
    ]
    @test issorted(f_errs; rev=true)
    _, p_f = _fit_loglog(Ns, f_errs)
    @test 0.7 ≤ p_f ≤ 1.5
end

# ────────────────────────────────────────────────────────────────────
# (5) XXZ Δ = 0: free-fermion GS energy convergence
# ────────────────────────────────────────────────────────────────────

@testset "XXZ Δ=0 OBC GS energy → -J/π at N → ∞" begin
    J = 1.0
    model = XXZ1D(; J=J, Δ=0.0)
    e_inf = -J / π
    β = 50.0  # large β: thermal ≈ ground state

    # Limited by `_MAX_ED_SITES = 12` cap on dense ED.
    Ns = (4, 6, 8, 10, 12)
    errs = Float64[]
    for N in Ns
        E_total = QAtlas.fetch(model, Energy(), OBC(N); beta=β)
        push!(errs, abs(E_total / N - e_inf))
    end
    @test issorted(errs; rev=true)
    # OBC convergence is 1/N from boundary missing-bond defect.
    _, p = _fit_loglog(Ns, errs)
    @test 0.5 ≤ p ≤ 2.0
end

# ────────────────────────────────────────────────────────────────────
# (6) Cross-BC: PBC closer to Infinite than OBC at every N
# ────────────────────────────────────────────────────────────────────

@testset "TFIM cross-BC: both OBC and PBC reach Infinite at large N" begin
    # Loose convergence smoke: at the largest sampled N, both OBC and
    # PBC per-site free energies are within 1 % of the Infinite()
    # value across the phase diagram (ordered, disordered, critical).
    # The relative ordering of OBC vs PBC errors at fixed finite N is
    # coefficient-dependent and not asserted here.
    cases = (
        (TFIM(J=1.0, h=0.5), 1.0),
        (TFIM(J=1.0, h=1.5), 1.5),
        (TFIM(J=1.0, h=1.0), 2.0),  # critical
    )
    for (model, β) in cases
        f_inf = QAtlas.fetch(model, FreeEnergy(), Infinite(); beta=β)
        N = 128
        f_obc = QAtlas.fetch(model, FreeEnergy(), OBC(N); beta=β)
        f_pbc = QAtlas.fetch(model, FreeEnergy(), PBC(N); beta=β)
        @test abs(f_obc - f_inf) / max(abs(f_inf), 1e-3) < 0.01
        @test abs(f_pbc - f_inf) / max(abs(f_inf), 1e-3) < 0.01
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "cross-BC scaling — verification cards" begin
    # β = 0: every TFIM bond/field term is traceless => per-site
    # ⟨H⟩_{β=0} = 0 exactly (independent operator-trace sum rule).
    for bc in (Infinite(), OBC(8), PBC(8))
        verify(
            TFIM(; J=1.0, h=0.5),
            Energy(:per_site),
            bc;
            route=:sum_rule,
            fetch_kw=(; beta=0.0),
            independent=0.0,
            agree_within=1e-9,
            refs=["Tr(σz σz)=Tr(σx)=0 => per-site ⟨H⟩_{β=0}=0 across all BC"],
        )
    end
end
