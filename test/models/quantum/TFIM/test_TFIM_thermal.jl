# =============================================================================
# Tests for TFIM finite-temperature observables (TFIM_thermal.jl + the
# thermal extensions of TFIM_dynamics.jl).
#
# Layers:
#   1. T → 0 limit: every thermal observable should match the T = 0
#      result already covered by the energy / sz_sz_correlation entries.
#   2. T → ∞ (high-temperature) analytic limits:
#         f(β→0)  → -(1/β) log 2          (single spin entropy log 2 per site)
#         s(β→0)  → log 2
#         c_v(β→0)→ 0
#         m_x(β→0)→ 0
#         χ_xx(β→0)→ ?  (free-spin Curie law: limit of integral)
#   3. Thermodynamic identities:
#         ε(β) = f(β) + T s(β),                           T = 1/β
#         c_v(β) = -β² ∂²(βf)/∂β²  (numerical second derivative)
#         χ_xx(β) = ∂m_x/∂h        (numerical first derivative)
#   4. ED comparison at small N:
#         all observables match a brute-force diagonalisation of the dense
#         2^N × 2^N TFIM Hamiltonian.
#   5. OBC ↔ Infinite consistency: per-site values for OBC at large N
#      converge to the Infinite values, and the deviation shrinks with N.
# =============================================================================

# ---- ED helpers (small N only) ---------------------------------------------

# Per-site thermodynamic averages from the dense spectrum at temperature β.
function _ed_thermo(N, J, h, β)
    H = _build_tfim_dense(N, J, h)
    E = eigvals(H)
    Z_shift = sum(exp.(-β .* (E .- E[1])))   # = e^{βE_min} · Z_true (stable)
    weights = exp.(-β .* (E .- E[1])) ./ Z_shift
    ε = sum(weights .* E) / N
    # log(Z_true) = -β E_min + log(Z_shift) ⇒ f = -T log(Z_true)/N = (E_min - T log(Z_shift)) / N
    f = (E[1] - log(Z_shift) / β) / N
    s = β * (ε - f)
    H2 = sum(weights .* E .^ 2)
    Havg = sum(weights .* E)
    c_v = β^2 * (H2 - Havg^2) / N
    return (; ε, f, s, c_v)
end

function _ed_transverse(N, J, h, β)
    H = _build_tfim_dense(N, J, h)
    E, V = eigen(H)
    Z = sum(exp.(-β .* (E .- E[1])))
    Mx = sum(_op_site(_SX, i, N) for i in 1:N)
    weights = exp.(-β .* (E .- E[1])) ./ Z
    diagMx = real.(diag(V' * Mx * V))
    mx_avg = sum(weights .* diagMx) / N
    return mx_avg
end

# ============================================================================

@testset "TFIM thermal observables" begin

    # ───────────────────────────────────────────────────────────────────────
    # Layer 1: T = 0 consistency with the existing entries.
    # Use the disordered phase (h > J) so the BdG gap is finite (= 2(h - J))
    # and there is no near-zero edge mode of the OBC ordered phase to spoil
    # the convergence at moderate β.
    # ───────────────────────────────────────────────────────────────────────
    @testset "T → 0 consistency" begin
        N, J, h = 10, 1.0, 1.5
        β_low = 80.0   # well below the gap (Δ = 2|h-J| = 1)

        # Per-site GS energy from BdG
        ε_thermal = QAtlas.fetch(TFIM(; J=J, h=h), Energy(), OBC(; N=N); beta=β_low) / N
        ε_gs = QAtlas.fetch(TFIM(; J=J, h=h), Energy(), OBC(; N=N)) / N
        @test ε_thermal ≈ ε_gs atol=1e-12

        # Per-site free energy → ε at T = 0
        f_thermal = QAtlas.fetch(TFIM(; J=J, h=h), FreeEnergy(), OBC(; N=N); beta=β_low)
        @test f_thermal ≈ ε_gs rtol=1e-6

        # Entropy → 0 at T = 0
        s_low = QAtlas.fetch(TFIM(; J=J, h=h), ThermalEntropy(), OBC(; N=N); beta=β_low)
        @test abs(s_low) < 1e-6

        # Specific heat → 0 at T = 0
        c_low = QAtlas.fetch(TFIM(; J=J, h=h), SpecificHeat(), OBC(; N=N); beta=β_low)
        @test abs(c_low) < 1e-6

        # Static σᶻ σᶻ thermal at β = ∞ matches the existing T = 0 entry —
        # this works in any phase since it does not rely on β being finite.
        for h_phase in (0.5, 1.0, 1.5)
            for r in 1:4
                v_th = QAtlas.fetch(
                    TFIM(; J=1.0, h=h_phase),
                    ZZCorrelation{:static}(),
                    OBC(; N=10);
                    beta=Inf,
                    i=2,
                    j=2 + r,
                )
                v_gs = real(
                    QAtlas.fetch(
                        TFIM(; J=1.0, h=h_phase),
                        ZZCorrelation{:dynamic}(),
                        OBC(; N=10);
                        i=2,
                        j=2 + r,
                        t=0.0,
                    ),
                )
                @test v_th ≈ v_gs atol=1e-12
            end
        end
    end

    # ───────────────────────────────────────────────────────────────────────
    # Layer 2: T → ∞ analytic limits (single-spin paramagnet).
    # ───────────────────────────────────────────────────────────────────────
    @testset "T → ∞ analytic limits" begin
        N, J, h = 8, 1.0, 0.5
        β_high = 1e-4

        # Free energy per site at β → 0:
        #   f → -T log 2 + O(β · h, β · J)
        f_hi = QAtlas.fetch(TFIM(; J=J, h=h), FreeEnergy(), OBC(; N=N); beta=β_high)
        @test f_hi ≈ -log(2) / β_high atol=1e-3 / β_high

        # Entropy per site at β → 0: s → log 2
        s_hi = QAtlas.fetch(TFIM(; J=J, h=h), ThermalEntropy(), OBC(; N=N); beta=β_high)
        @test s_hi ≈ log(2) atol=1e-6

        # Specific heat at β → 0: c_v → 0 quadratically in β
        c_hi = QAtlas.fetch(TFIM(; J=J, h=h), SpecificHeat(), OBC(; N=N); beta=β_high)
        @test abs(c_hi) < 1e-6

        # Transverse magnetisation at β → 0: m_x → 0 (Curie regime)
        mx_hi = QAtlas.fetch(TFIM(; J=J, h=h), MagnetizationX(), OBC(; N=N); beta=β_high)
        @test abs(mx_hi) < 1e-3

        # Same in the thermodynamic limit.
        f_inf = QAtlas.fetch(TFIM(; J=J, h=h), FreeEnergy(), Infinite(); beta=β_high)
        s_inf = QAtlas.fetch(TFIM(; J=J, h=h), ThermalEntropy(), Infinite(); beta=β_high)
        @test f_inf ≈ -log(2) / β_high atol=1e-3 / β_high
        @test s_inf ≈ log(2) atol=1e-6
    end

    # ───────────────────────────────────────────────────────────────────────
    # Layer 3: Thermodynamic identities.
    # ───────────────────────────────────────────────────────────────────────
    @testset "thermodynamic identities" begin
        N, J, h = 10, 1.0, 0.5
        for β in (0.3, 0.7, 1.5, 3.0)
            ε = QAtlas.fetch(TFIM(; J=J, h=h), Energy(), OBC(; N=N); beta=β) / N
            f = QAtlas.fetch(TFIM(; J=J, h=h), FreeEnergy(), OBC(; N=N); beta=β)
            s = QAtlas.fetch(TFIM(; J=J, h=h), ThermalEntropy(), OBC(; N=N); beta=β)

            # ε = f + T s = f + s/β
            @test ε ≈ f + s / β atol=1e-10

            # c_v from numerical second derivative of (β f).
            δ = 1e-3 * β
            f_p = QAtlas.fetch(TFIM(; J=J, h=h), FreeEnergy(), OBC(; N=N); beta=β + δ)
            f_m = QAtlas.fetch(TFIM(; J=J, h=h), FreeEnergy(), OBC(; N=N); beta=β - δ)
            βf = β * f
            βpfp = (β + δ) * f_p
            βmfm = (β - δ) * f_m
            c_num = -β^2 * (βpfp - 2 * βf + βmfm) / δ^2
            c_an = QAtlas.fetch(TFIM(; J=J, h=h), SpecificHeat(), OBC(; N=N); beta=β)
            @test c_an ≈ c_num rtol=1e-3
        end
    end

    @testset "thermodynamic identities (Infinite)" begin
        J, h = 1.0, 0.7
        for β in (0.3, 0.7, 1.5, 3.0)
            ε = QAtlas.fetch(TFIM(; J=J, h=h), Energy(), Infinite(); beta=β)
            f = QAtlas.fetch(TFIM(; J=J, h=h), FreeEnergy(), Infinite(); beta=β)
            s = QAtlas.fetch(TFIM(; J=J, h=h), ThermalEntropy(), Infinite(); beta=β)
            @test ε ≈ f + s / β atol=1e-10

            δ = 1e-3 * β
            f_p = QAtlas.fetch(TFIM(; J=J, h=h), FreeEnergy(), Infinite(); beta=β + δ)
            f_m = QAtlas.fetch(TFIM(; J=J, h=h), FreeEnergy(), Infinite(); beta=β - δ)
            c_num = -β^2 * ((β + δ) * f_p - 2 * β * f + (β - δ) * f_m) / δ^2
            c_an = QAtlas.fetch(TFIM(; J=J, h=h), SpecificHeat(), Infinite(); beta=β)
            @test c_an ≈ c_num rtol=1e-3

            δh = 1e-4
            mp = QAtlas.fetch(TFIM(; J=J, h=h + δh), MagnetizationX(), Infinite(); beta=β)
            mm = QAtlas.fetch(TFIM(; J=J, h=h - δh), MagnetizationX(), Infinite(); beta=β)
            χ_num = (mp - mm) / (2 * δh)
            χ_an = QAtlas.fetch(TFIM(; J=J, h=h), SusceptibilityXX(), Infinite(); beta=β)
            @test χ_an ≈ χ_num rtol=5e-3
        end
    end

    # ───────────────────────────────────────────────────────────────────────
    # Layer 4: ED comparison at small N (per site).
    # ───────────────────────────────────────────────────────────────────────
    @testset "ED comparison N=4" begin
        N, J, h = 4, 1.0, 0.7
        for β in (0.5, 1.0, 2.5)
            ed = _ed_thermo(N, J, h, β)
            mx_ed = _ed_transverse(N, J, h, β)

            ε = QAtlas.fetch(TFIM(; J=J, h=h), Energy(), OBC(; N=N); beta=β) / N
            f = QAtlas.fetch(TFIM(; J=J, h=h), FreeEnergy(), OBC(; N=N); beta=β)
            s = QAtlas.fetch(TFIM(; J=J, h=h), ThermalEntropy(), OBC(; N=N); beta=β)
            c = QAtlas.fetch(TFIM(; J=J, h=h), SpecificHeat(), OBC(; N=N); beta=β)
            mx = QAtlas.fetch(TFIM(; J=J, h=h), MagnetizationX(), OBC(; N=N); beta=β)

            @test ε ≈ ed.ε atol=1e-10
            @test f ≈ ed.f atol=1e-10
            @test s ≈ ed.s atol=1e-10
            @test c ≈ ed.c_v atol=1e-10
            @test mx ≈ mx_ed atol=1e-10
        end
    end

    @testset "ED comparison: static σᶻ σᶻ thermal" begin
        N, J, h = 4, 1.0, 0.7
        H = _build_tfim_dense(N, J, h)
        E, V = eigen(H)
        for β in (0.5, 1.0, 2.5)
            Z = sum(exp.(-β .* (E .- E[1])))
            ρ = V * (Diagonal(exp.(-β .* (E .- E[1]))) ./ Z) * V'
            for i in 1:N, j in i:N
                op = _op_site(_SZ, i, N) * _op_site(_SZ, j, N)
                ed_val = real(tr(ρ * op))
                qa_val = QAtlas.fetch(
                    TFIM(; J=J, h=h), ZZCorrelation{:static}(), OBC(; N=N); beta=β, i=i, j=j
                )
                @test qa_val ≈ ed_val atol=1e-10
            end
        end
    end

    @testset "ED comparison: longitudinal susceptibility" begin
        N, J, h = 4, 1.0, 0.7
        for β in (0.5, 1.0, 2.5)
            # χ_zz = β · ⟨M²⟩_c / N from the ED density matrix
            H = _build_tfim_dense(N, J, h)
            E, V = eigen(H)
            Z = sum(exp.(-β .* (E .- E[1])))
            ρ = V * (Diagonal(exp.(-β .* (E .- E[1]))) ./ Z) * V'
            Mz = sum(_op_site(_SZ, i, N) for i in 1:N)
            M2 = real(tr(ρ * (Mz * Mz)))
            M1 = real(tr(ρ * Mz))
            χ_ed = β * (M2 - M1^2) / N

            χ_qa = QAtlas.fetch(TFIM(; J=J, h=h), SusceptibilityZZ(), OBC(; N=N); beta=β)
            @test χ_qa ≈ χ_ed atol=1e-10
        end
    end

    @testset "ED comparison: transverse susceptibility" begin
        # SusceptibilityXX(OBC) returns β·Var(Mx)/N (equal-time fluctuation),
        # which for a quantum system differs from the Kubo response ∂m_x/∂h.
        # Validate against the same ED density-matrix variance formula.
        N, J, h = 4, 1.0, 0.7
        for β in (0.5, 1.0, 2.5)
            H = _build_tfim_dense(N, J, h)
            E, V = eigen(H)
            Z = sum(exp.(-β .* (E .- E[1])))
            ρ = V * (Diagonal(exp.(-β .* (E .- E[1]))) ./ Z) * V'
            Mx = sum(_op_site(_SX, i, N) for i in 1:N)
            M2 = real(tr(ρ * (Mx * Mx)))
            M1 = real(tr(ρ * Mx))
            χ_ed = β * (M2 - M1^2) / N

            χ_qa = QAtlas.fetch(TFIM(; J=J, h=h), SusceptibilityXX(), OBC(; N=N); beta=β)
            @test χ_qa ≈ χ_ed atol=1e-10
        end
    end

    # ───────────────────────────────────────────────────────────────────────
    # Layer 5: OBC ↔ Infinite consistency at large N.
    # ───────────────────────────────────────────────────────────────────────
    @testset "OBC → Infinite at large N" begin
        J, h, β = 1.0, 0.5, 1.0
        Ns = (20, 60, 120)
        errs_f = Float64[]
        errs_s = Float64[]
        errs_c = Float64[]
        errs_m = Float64[]
        for N in Ns
            f_obc = QAtlas.fetch(TFIM(; J=J, h=h), FreeEnergy(), OBC(; N=N); beta=β)
            s_obc = QAtlas.fetch(TFIM(; J=J, h=h), ThermalEntropy(), OBC(; N=N); beta=β)
            c_obc = QAtlas.fetch(TFIM(; J=J, h=h), SpecificHeat(), OBC(; N=N); beta=β)
            m_obc = QAtlas.fetch(TFIM(; J=J, h=h), MagnetizationX(), OBC(; N=N); beta=β)
            f_inf = QAtlas.fetch(TFIM(; J=J, h=h), FreeEnergy(), Infinite(); beta=β)
            s_inf = QAtlas.fetch(TFIM(; J=J, h=h), ThermalEntropy(), Infinite(); beta=β)
            c_inf = QAtlas.fetch(TFIM(; J=J, h=h), SpecificHeat(), Infinite(); beta=β)
            m_inf = QAtlas.fetch(TFIM(; J=J, h=h), MagnetizationX(), Infinite(); beta=β)
            push!(errs_f, abs(f_obc - f_inf))
            push!(errs_s, abs(s_obc - s_inf))
            push!(errs_c, abs(c_obc - c_inf))
            push!(errs_m, abs(m_obc - m_inf))
        end
        # All errors should be small at the largest N
        @test errs_f[end] < 0.05
        @test errs_s[end] < 0.05
        @test errs_c[end] < 0.05
        @test errs_m[end] < 0.05
        # Convergence: largest-N error should be the smallest
        @test errs_f[end] ≤ errs_f[1]
        @test errs_s[end] ≤ errs_s[1]
        @test errs_c[end] ≤ errs_c[1]
        @test errs_m[end] ≤ errs_m[1]
    end
end
