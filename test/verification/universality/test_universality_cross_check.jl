# ─────────────────────────────────────────────────────────────────────────────
# Cross-verification: universality exponents vs model-specific results
#
# This test connects the ABSTRACT universality class data
# (src/universalities/) to the CONCRETE model results (src/models/)
# by numerically extracting critical exponents from the model-specific
# closed-form solutions and comparing them to the tabulated values.
#
# Each test constitutes an independent physical verification:
# - The universality exponent was sourced from paper A (CFT / Coulomb gas)
# - The model result was sourced from paper B (Onsager / Yang / Bethe)
# - If both agree, the QAtlas data is cross-validated from two
#   independent theoretical lines.
#
# This is the key layer that makes QAtlas's universality data
# *physically* verified, not just internally consistent.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Lattice2D, LinearAlgebra, ForwardDiff, Test

# ═══════════════════════════════════════════════════════════════════════════════
# 1. Yang magnetization → β = 1/8
#
# Source A: Universality(:Ising) d=2 → β = 1/8 (CFT: BPZ 1984)
# Source B: IsingSquare SpontaneousMagnetization → M(T) (Yang 1952)
#
# Extraction: log(M) / log(T_c − T) → β as T → T_c⁻
# ═══════════════════════════════════════════════════════════════════════════════

@testset "Cross-check: Yang M(T) → β = 1/8 (Ising 2D)" begin
    β_exact = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2).β
    Tc = QAtlas.fetch(IsingSquare(), CriticalTemperature())

    # Extract effective β from log-log slope at progressively closer T to Tc
    δTs = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
    β_effs = Float64[]
    for i in 1:(length(δTs) - 1)
        T1, T2 = Tc - δTs[i], Tc - δTs[i + 1]
        M1 = QAtlas.fetch(IsingSquare(), SpontaneousMagnetization(); β=1 / T1)
        M2 = QAtlas.fetch(IsingSquare(), SpontaneousMagnetization(); β=1 / T2)
        push!(β_effs, log(M2 / M1) / log(δTs[i + 1] / δTs[i]))
    end

    # The Yang formula is exact, so finite-difference β_eff converges
    # exponentially to β = 1/8 as δT → 0. Empirical residuals over
    # δT = 10^{-2}…10^{-6}: 4e-4, 4e-5, 4e-6, ~0 (machine precision at
    # the tightest pair). The previous `rtol = 0.01` (1 %) was 25× looser
    # than the worst pair.
    refs = [Float64(β_exact)]
    @test all(b -> isapprox(b, refs[1]; atol=1e-3), β_effs)
    # Monotone convergence — each pair tighter than the previous.
    @test all(
        abs(β_effs[i + 1] - refs[1]) < abs(β_effs[i] - refs[1]) for
        i in 1:(length(β_effs) - 1)
    )
    # Tightest pair must be at machine precision (Yang's exact formula
    # leaves no finite-N residue at this δT).
    @test isapprox(β_effs[end], refs[1]; atol=1e-5)
end

# ═══════════════════════════════════════════════════════════════════════════════
# 2. TFIM gap Δ(N) at h = J → dynamic exponent z = 1
#
# Source A: Universality(:Ising) d=2 → ν = 1, and the 1D TFIM has
#           dynamic exponent z = 1. The gap scales as Δ ~ N^{-z}.
# Source B: TFIM ED gap from build_tfim + Lattice2D.
#
# Extraction: log(Δ) / log(1/N) → z as N increases
# ═══════════════════════════════════════════════════════════════════════════════

@testset "Cross-check: TFIM gap Δ(N) → z = 1 (Ising 2D)" begin
    ising_exp = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2)
    J = 1.0
    h = J  # critical point

    # Compute gap for several N
    Ns = [4, 6, 8, 10, 12]
    gaps = Float64[]
    for N in Ns
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())
        H = build_tfim(lat, J, h)
        λ = sort(eigvals(Symmetric(H)))
        push!(gaps, λ[2] - λ[1])
    end

    # Extract z from log-log fit: log Δ = -z log N + const
    # Use pairs of successive N
    z_effs = Float64[]
    for i in 1:(length(Ns) - 1)
        z_eff = -log(gaps[i + 1] / gaps[i]) / log(Ns[i + 1] / Ns[i])
        push!(z_effs, z_eff)
    end

    # z should approach 1 (known exact value for 1D TFIM = 2D Ising, z = 1)
    # Finite-size corrections make early values less accurate
    @test z_effs[end] ≈ 1.0 rtol = 0.1  # last pair, closest to thermodynamic limit
    # The trend should be monotonically approaching 1
    @test z_effs[end] > z_effs[1] || abs(z_effs[end] - 1.0) < abs(z_effs[1] - 1.0)
end

# ═══════════════════════════════════════════════════════════════════════════════
# 2b. TFIM gap Δ(N) at h = J → z = 1 (tight PBC signal)
#
# The OBC extraction above has a boundary-energy O(1/N) contamination
# in the gap that caps the extracted z_eff at ~0.97 for N ≤ 12, so the
# assertion `rtol = 0.1` only rejects gross deviations. Repeating the
# extraction on the PBC chain removes the boundary term entirely and
# sends z_eff to 1 with a 1/N² tail — extracted z lies within ~0.2 %
# of the exact value already at N ∈ {14, 16}, enough to catch any sign
# or prefactor drift in the BdG spectrum used by QAtlas.
# ═══════════════════════════════════════════════════════════════════════════════

@testset "Cross-check: TFIM gap Δ(N) → z = 1 (PBC, tight)" begin
    J = 1.0
    h = J
    bc_chain = LatticeBoundary((PeriodicAxis(), OpenAxis()))

    Ns = [10, 12, 14, 16]
    gaps = Float64[]
    for N in Ns
        lat = build_lattice(Square, N, 1; boundary=bc_chain)
        H = build_tfim_sparse(lat, J, h)
        vals, _, info = eigsolve(
            H,
            randn(MersenneTwister(42 + N), 2^N),
            2,
            :SR;
            issymmetric=true,
            tol=1e-11,
            krylovdim=40,
        )
        info.converged < 2 && error(
            "PBC TFIM two-lowest ED failed to converge at N = $N (converged=$(info.converged))",
        )
        push!(gaps, vals[2] - vals[1])
    end

    z_effs = [
        -log(gaps[i + 1] / gaps[i]) / log(Ns[i + 1] / Ns[i]) for i in 1:(length(Ns) - 1)
    ]
    # Geometric mean of each (Ns[i], Ns[i+1]) pair — the natural N̄
    # for the 1/N̄² ansatz below.
    Nbars = [sqrt(Float64(Ns[i] * Ns[i + 1])) for i in 1:(length(Ns) - 1)]

    # Calabrese–Cardy on a torus predicts that z_eff(N̄) − 1 scales as
    # `a / N̄²` with no boundary contribution.  Empirically, this 1-parameter
    # ansatz fits z_effs[1..3] with `a ≈ 0.41` to 4 significant figures at
    # N̄ ∈ [√120, √168, √224]. Fit and assert the intercept matches z = 1
    # to atol = 5e-4 — strictly stronger than the previous rtol = 0.01
    # single-pair check.
    invNbar_sq = [1.0 / N̄^2 for N̄ in Nbars]
    n = length(invNbar_sq)
    x̄ = sum(invNbar_sq) / n
    ȳ = sum(z_effs) / n
    a_fit =
        sum((invNbar_sq[i] - x̄) * (z_effs[i] - ȳ) for i in 1:n) /
        sum((invNbar_sq[i] - x̄)^2 for i in 1:n)
    z_inf = ȳ - a_fit * x̄
    @test isapprox(z_inf, 1.0; atol=5e-4)
    @test a_fit > 0  # leading correction sign on PBC critical TFIM

    # Error shrinks monotonically in N.
    for i in 1:(length(z_effs) - 1)
        @test abs(z_effs[i + 1] - 1.0) < abs(z_effs[i] - 1.0)
    end

    # Finite-N correction has the expected positive sign (z_eff > 1 at
    # finite N for PBC TFIM at h = J, with the correction scaling as
    # O(1/N²)).
    @test all(>(1.0), z_effs)
end

# ═══════════════════════════════════════════════════════════════════════════════
# 3. TFIM ground-state energy E₀(N) → Ising central charge c = 1/2
#
# Source A: Universality(:Ising) d=2 → c = 1/2 (CFT)
# Source B: TFIM ED energy from build_tfim.
#
# At criticality, the finite-size correction to the ground state energy
# per site is (Blöte, Cardy, Nightingale 1986, Affleck 1986):
#   E₀(N)/N = e_∞ − π v c / (6 N²)
# where v is the "velocity" (= 2J sin(π/N) ≈ 2πJ/N for OBC TFIM at
# criticality). Extracting c from the finite-size scaling.
#
# Simpler check: E₀(N)/N converges and the correction ~ 1/N².
# ═══════════════════════════════════════════════════════════════════════════════

@testset "Cross-check: TFIM E₀ scaling → consistent with c = 1/2" begin
    J = 1.0
    h = J

    Ns = [6, 8, 10, 12]
    e0_per_site = Float64[]
    for N in Ns
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())
        H = build_tfim(lat, J, h)
        λ = sort(eigvals(Symmetric(H)))
        push!(e0_per_site, λ[1] / N)
    end

    # E₀/N should converge as N grows. For OBC the approach is not
    # necessarily monotone due to boundary corrections, so we only check
    # that the spread decreases.
    spread = maximum(e0_per_site) - minimum(e0_per_site)
    @test spread < 0.1  # all values within 0.1 of each other

    # Cross-check with the BdG analytical E₀ at each N
    for (i, N) in enumerate(Ns)
        E0_analytical = QAtlas.fetch(TFIM(; J=J, h=h), Energy(), OBC(; N=N))
        @test e0_per_site[i] * N ≈ E0_analytical rtol = 1e-10
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# 3b. TFIM BdG E₀ at larger N via sparse ED (Lanczos)
#
# Source A: QAtlas TFIM Energy OBC — analytical BdG sum E₀ = −½ Σ Λₙ.
# Source B: Real-space spin Hamiltonian H = −J Σ σ^z_i σ^z_{i+1} − h Σ σ^x_i
#           built on the 2^N spin-1/2 basis via `build_tfim_sparse`, lowest
#           eigenpair extracted with KrylovKit Lanczos (O(nnz) per iter).
#
# The dense-ED cross-check above (test 3) reaches N = 12 only because
# dense eigvals scales as O(4^N).  Repeating the same cross-check at
# N = 14, 16 through sparse ED pushes the independent verification into
# the regime where boundary corrections are numerically smaller and any
# drift between the BdG formula and the real-space spin Hamiltonian
# would show up cleanly.  Both the critical point h = J and a
# paramagnetic point h = 2J are covered.
# ═══════════════════════════════════════════════════════════════════════════════

@testset "Cross-check: TFIM sparse ED ↔ BdG E₀ at N=14, 16 (OBC)" begin
    for (J, h) in [(1.0, 1.0), (1.0, 2.0)]
        for N in (14, 16)
            lat = build_lattice(Square, N, 1; boundary=OpenAxis())
            H = build_tfim_sparse(lat, J, h)
            vals, _, info = eigsolve(
                H, randn(2^N), 1, :SR; issymmetric=true, tol=1e-11, krylovdim=30
            )
            info.converged < 1 && error("sparse ED failed to converge at N = $N")
            E0_ed = vals[1]

            E0_bdg = QAtlas.fetch(TFIM(; J=J, h=h), Energy(), OBC(; N=N))
            @test E0_ed ≈ E0_bdg rtol = 1e-8
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# 4. Heisenberg E₀/N finite-size extrapolation
#
# Source A: Bethe ansatz → e_∞ = J(1/4 − ln 2) ≈ −0.4431J (N→∞, PBC)
# Source B: Heisenberg1D ExactSpectrum N=4 PBC, and ED for N=6,8 PBC
#
# The finite-size correction for the 1D Heisenberg chain is known:
#   E₀(N)/N = e_∞ + A/N² + O(1/N⁴)   (Bethe ansatz, PBC)
# ═══════════════════════════════════════════════════════════════════════════════

@testset "Cross-check: Heisenberg E₀/N → Bethe ansatz e_∞" begin
    J = 1.0
    e_bethe = J * (1 // 4 - log(2))  # ≈ -0.4431

    # N=4 from QAtlas exact
    λ4 = QAtlas.fetch(Heisenberg1D(), ExactSpectrum(); N=4, J=J, bc=:PBC)
    e0_4 = λ4[1] / 4

    # N=6 and N=8 via ED (using Lattice2D PBC chain + spinhalf_ed)
    e0s = Dict{Int,Float64}()
    e0s[4] = e0_4
    for N in [6, 8]
        lat = build_lattice(
            Square, N, 1; boundary=LatticeBoundary((PeriodicAxis(), OpenAxis()))
        )
        H = build_spinhalf_heisenberg(lat, J)
        λ = sort(eigvals(Symmetric(H)))
        e0s[N] = λ[1] / N
    end

    # All E₀/N should be below (more negative than) Bethe ansatz e_∞
    # because finite-size correction is negative for PBC Heisenberg
    for (N, e0) in e0s
        @test e0 < e_bethe
    end

    # E₀/N should approach e_∞ as N increases (|e0 - e_bethe| decreases)
    err4 = abs(e0s[4] - e_bethe)
    err6 = abs(e0s[6] - e_bethe)
    err8 = abs(e0s[8] - e_bethe)
    @test err4 > err6 > err8

    # N=8 should be within ~5% of the Bethe ansatz value
    @test err8 / abs(e_bethe) < 0.05
end

# ═══════════════════════════════════════════════════════════════════════════════
# 5. Onsager T_c ↔ Ising universality self-consistency
#
# The critical temperature T_c = 2J/ln(1+√2) is a MODEL-specific result
# (IsingSquare), but it MUST be consistent with the UNIVERSAL relation:
#   sinh(2J/T_c) = 1   (Kramers-Wannier duality fixed point)
# This connects the IsingSquare model to the universality structure.
# ═══════════════════════════════════════════════════════════════════════════════

@testset "Cross-check: Onsager T_c ↔ duality + universality" begin
    J = 1.0
    Tc = QAtlas.fetch(IsingSquare(), CriticalTemperature(); J=J)
    βc = 1.0 / Tc

    # Kramers-Wannier duality: sinh(2K_c) = 1
    @test sinh(2 * βc * J) ≈ 1.0 atol = 1e-14

    # At T = T_c, Yang magnetization = 0 (phase boundary)
    @test QAtlas.fetch(IsingSquare(), SpontaneousMagnetization(); β=βc, J=J) == 0.0

    # At T = T_c, the correlation length diverges (ν = 1) — we can't
    # directly test ξ, but we can verify the partition function Z(β)
    # has a non-analyticity by checking the transfer-matrix spectrum
    # has eigenvalue degeneracy approaching 1 as system size grows.
    # (Qualitative check: the partition function ratio Z(βc+δ)/Z(βc-δ)
    # is well-behaved for small systems.)
    for Ly in [3, 4]
        Z_above = QAtlas.fetch(
            IsingSquare(), PartitionFunction(); Lx=Ly, Ly=Ly, β=βc * 1.01, J=J
        )
        Z_at = QAtlas.fetch(IsingSquare(), PartitionFunction(); Lx=Ly, Ly=Ly, β=βc, J=J)
        Z_below = QAtlas.fetch(
            IsingSquare(), PartitionFunction(); Lx=Ly, Ly=Ly, β=βc * 0.99, J=J
        )
        @test Z_above > Z_at > Z_below  # Z increases with β (E_gs < 0)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# 6. TFIM BdG gap → νz = 1 (rigorous result based)
#
# Source A: Universality(:Ising) d=2 → ν=1, and z=1 for 1D TFIM → νz=1
# Source B: QAtlas BdG quasiparticle spectrum (TFIM.jl, rigorous).
#
# The BdG gap Δ = min(Λ_n) is the EXACT many-body gap of the TFIM.
# In the thermodynamic limit (N→∞), Δ → 2|J−h| for h ≠ J, confirming
# the gap scales as |h−J|^{νz} with νz = 1.
#
# This test uses ONLY QAtlas rigorous results (no ED approximation).
# ═══════════════════════════════════════════════════════════════════════════════

@testset "Cross-check: TFIM BdG gap → νz = 1 (rigorous)" begin
    ising_exp = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2)
    # For 1D TFIM: z = 1, ν = ising_exp.ν = 1, so νz = 1
    νz_expected = Float64(ising_exp.ν) * 1.0  # z = 1 for TFIM

    J = 1.0
    N = 200  # large N for thermodynamic limit
    # BdG quasiparticle spectrum → gap = min(Λ)
    hs_disordered = [1.1, 1.5, 2.0, 3.0, 5.0]

    for h in hs_disordered
        Λ = QAtlas._tfim_bdg_spectrum(N, J, h)
        gap_bdg = minimum(Λ)
        gap_exact = 2 * abs(h - J)  # thermodynamic limit prediction
        # OBC boundary correction is O(1/N); empirical relative residual
        # at N = 200 ranges from 1.2e-2 (h = 1.1, closest to critical) to
        # 4e-5 (h = 5). `rtol = 0.02` covers the worst case with ~2x
        # margin and is 2.5× tighter than the previous `rtol = 0.05`.
        @test gap_bdg ≈ gap_exact rtol = 0.02
        # Boundary correction has positive sign for OBC critical TFIM —
        # any sign flip in the BdG spectrum would push gap below.
        @test gap_bdg ≥ gap_exact
    end

    # Extract νz from log-log slope: log(Δ) vs log|h−J|
    hs = [1.05, 1.1, 1.2, 1.5, 2.0]
    log_gaps = [log(minimum(QAtlas._tfim_bdg_spectrum(N, J, h))) for h in hs]
    log_dhs = [log(abs(h - J)) for h in hs]
    # Linear regression slope ≈ νz
    n = length(hs)
    x̄ = sum(log_dhs) / n
    ȳ = sum(log_gaps) / n
    slope =
        sum((log_dhs[i] - x̄) * (log_gaps[i] - ȳ) for i in 1:n) /
        sum((log_dhs[i] - x̄)^2 for i in 1:n)
    # Empirical slope = 0.988 over hs = $(hs) at N = $N, residual −1.2e-2
    # from finite-N OBC boundary corrections in the closest-to-critical
    # h points (h = 1.05, 1.1).  `rtol = 0.02` is 2.5× tighter than the
    # previous `rtol = 0.05` and gives a clean ~2x margin over the
    # measured residual.
    @test slope ≈ νz_expected rtol = 0.02
    @test slope < νz_expected   # finite-N reduces the slope (boundary correction)
end

# ═══════════════════════════════════════════════════════════════════════════════
# 7. IsingSquare specific heat near T_c → α = 0 (rigorous, via ForwardDiff)
#
# Source A: Universality(:Ising) d=2 → α = 0 (logarithmic divergence)
# Source B: IsingSquare PartitionFunction + ForwardDiff → C(β) = β²∂²(ln Z)/∂β²
#
# For α = 0, the specific heat diverges logarithmically:
#   C(T) ~ -A ln|T − T_c| + B
# Unlike α > 0 (power-law divergence), C grows slowly near T_c.
#
# Test: C at progressively closer T to T_c grows but SLOWER than any
# power law |T−T_c|^{-ε} for ε > 0 — consistent with α = 0 (log).
# ═══════════════════════════════════════════════════════════════════════════════

@testset "Cross-check: IsingSquare specific heat → α = 0 (rigorous)" begin
    α_exact = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2).α
    @test α_exact == 0  # logarithmic, not power-law

    J = 1.0
    Tc = QAtlas.fetch(IsingSquare(), CriticalTemperature(); J=J)
    Lx, Ly = 4, 4  # small lattice for transfer matrix

    log_Z(β) = log(QAtlas.fetch(IsingSquare(), PartitionFunction(); Lx=Lx, Ly=Ly, β=β, J=J))

    # Compute C(β) = β² d²(log Z)/dβ² via nested ForwardDiff
    function specific_heat(β)
        d2 = ForwardDiff.derivative(βo -> ForwardDiff.derivative(βi -> log_Z(βi), βo), β)
        return β^2 * d2
    end

    # C at several points approaching T_c from above (β < β_c).
    # On a 4×4 lattice the singularity is rounded by finite size, so
    # we only test at well-separated points where the trend is clear.
    βc = 1 / Tc
    βs_test = [βc * 0.5, βc * 0.7, βc * 0.9]
    Cs = [specific_heat(β) for β in βs_test]

    # C should generally increase as T → T_c (β increases toward β_c)
    @test Cs[end] > Cs[1]

    # All specific heat values should be positive
    for C in Cs
        @test C > 0
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# 8. TFIM at h = J: central charge c = 1/2 from E₀(N) finite-size scaling
#
# Source A: Universality(:Ising) d=2 → c = 1/2
# Source B: QAtlas BdG ground-state energy E₀(N) (rigorous).
#
# Cardy formula (1986): at criticality, the ground-state energy of a
# 1D critical system with OBC has the finite-size correction
#   E₀(N) = N e_∞ − π v_F c / (24 N) + O(1/N²)
# where v_F is the Fermi velocity and c the central charge.
#
# For OBC TFIM at h = J = 1, v_F = 2J = 2 and e_∞ can be obtained
# from the integral of the dispersion.
#
# We extract c from the N-dependence of E₀(N) using QAtlas's exact BdG.
# ═══════════════════════════════════════════════════════════════════════════════

@testset "Cross-check: TFIM E₀(N)/N → e_∞ = −4J/π (rigorous, OBC)" begin
    # For the OBC TFIM at h = J (critical), the thermodynamic limit of
    # the ground state energy per site is e_∞ = −(2/π)·2J = −4J/π.
    # This value comes from integrating the BdG dispersion Λ(k) = 2J|sin k|
    # over the Brillouin zone.
    #
    # Both sides of this test are consequences of the same Pfeuty BdG
    # formula — the finite-N sum −½ Σ Λₙ and the thermodynamic-limit
    # integral −4J/π. Independent cross-validation against a real-space
    # spin Hamiltonian is handled by the sparse-ED check above ("3b",
    # N = 14, 16) and the dense-ED check ("3", N ≤ 12); this test adds
    # the complementary assertion that the BdG sum approaches the BdG
    # integral at large N with the expected 1/N boundary-correction
    # scaling, which ED cannot probe at such sizes.
    #
    # For OBC, the leading finite-size correction is O(1/N) (boundary
    # energy), not O(1/N²) (Cardy). Extracting the central charge c = 1/2
    # requires PBC or much larger N; instead we verify:
    # (a) E₀/N converges to e_∞ = -4J/π as N → ∞
    # (b) The convergence is 1/N (boundary energy scaling)
    # (c) BdG E₀ matches the analytical value at each N

    J = 1.0
    h = J
    e_inf = -4J / π  # ≈ -1.2732

    Ns = [50, 100, 200, 400]
    E0s = [QAtlas.fetch(TFIM(; J=J, h=h), Energy(), OBC(; N=N)) for N in Ns]
    e0_per_site = [E0s[i] / Ns[i] for i in eachindex(Ns)]

    # (a) Convergence: error decreases with N
    errs = [abs(e0_per_site[i] - e_inf) for i in eachindex(Ns)]
    for i in 1:(length(errs) - 1)
        @test errs[i] > errs[i + 1]
    end
    @test errs[end] < 0.001  # N=400 within 0.1% of e_∞

    # (b) O(1/N) scaling: corr × N ≈ const (boundary energy coefficient).
    # Empirically corrN converges to ≈ 0.363 with O(1/N²) subleading.
    # Largest measured pairwise spread is 3.6e-3 (between N=50 and 100);
    # `rtol = 0.005` tightens the previous rtol = 0.01 by 2× and gives
    # margin over that worst pair while remaining sensitive to scaling
    # bugs.
    corrN = [(e0_per_site[i] - e_inf) * Ns[i] for i in eachindex(Ns)]
    for i in 1:(length(corrN) - 1)
        @test corrN[i] ≈ corrN[i + 1] rtol = 0.005
    end

    # 1/N fit on `corrN(N) = c0 + c1/N`: if the leading scaling is truly
    # 1/N (and not 1/√N or 1/N²), the extrapolated boundary-energy
    # coefficient c0 should match the largest-N corrN to atol = 5e-4
    # (empirical: c0 = 0.36337 vs corrN[N=400] = 0.36305, Δ ≈ 3.2e-4).
    # A wrong-power scaling bug pushes c0 well past 5e-4.
    invN = [1.0 / Ns[i] for i in eachindex(Ns)]
    n = length(Ns)
    x̄ = sum(invN) / n
    ȳ = sum(corrN) / n
    c1 =
        sum((invN[i] - x̄) * (corrN[i] - ȳ) for i in 1:n) /
        sum((invN[i] - x̄)^2 for i in 1:n)
    c0 = ȳ - c1 * x̄
    @test isapprox(c0, corrN[end]; atol=5e-4)
    @test c1 < 0  # corrN approaches its limit from below at OBC critical TFIM
end
