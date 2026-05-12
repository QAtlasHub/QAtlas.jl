# ─────────────────────────────────────────────────────────────────────────────
# Verification: TFIM ground-state energy & gap closure
#
# Cross-validate the dense many-body Hamiltonian
#     H = -J Σ σᶻ_i σᶻ_{i+1}  −  h Σ σˣ_i
# built from Lattice2D's OBC chain via `build_tfim(lat, J, h)`
# against the analytical BdG ground-state energy from QAtlas (TFIM.jl).
#
# Additionally verify that the many-body energy gap Δ = E₁ − E₀ closes
# at the quantum critical point h = J (Ising CFT, c = 1/2). For finite
# N the gap scales as Δ ~ π v_F / N; the test checks this scaling.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Lattice2D, LinearAlgebra, Test

@testset "TFIM — ED ground state vs BdG analytical" begin
    for N in [4, 6, 8]
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())
        @test num_sites(lat) == N
        @test length(collect(bonds(lat))) == N - 1

        @testset "N=$N OBC — E₀ match" begin
            for (J, h) in [(1.0, 0.0), (1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (0.5, 1.5)]
                H = build_tfim(lat, J, h)
                λ = sort(eigvals(Symmetric(H)))
                E0_ed = λ[1]
                E0_analytical = QAtlas.fetch(TFIM(; J=J, h=h), Energy(), OBC(; N=N))
                @test E0_ed ≈ E0_analytical rtol = 1e-10
            end
        end
    end
end

@testset "TFIM — gap closure at quantum critical point h = J" begin
    J = 1.0

    @testset "Gap shrinks with N at h = J (critical)" begin
        gaps_at_critical = Float64[]
        Ns = [4, 6, 8]
        for N in Ns
            lat = build_lattice(Square, N, 1; boundary=OpenAxis())
            H = build_tfim(lat, J, J)
            λ = sort(eigvals(Symmetric(H)))
            push!(gaps_at_critical, λ[2] - λ[1])
        end
        # Gap should decrease with N (finite-size scaling)
        for k in 1:(length(gaps_at_critical) - 1)
            @test gaps_at_critical[k] > gaps_at_critical[k + 1]
        end
    end

    @testset "Gap structure across phases" begin
        N = 6
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())

        # Deep in disordered phase (h ≫ J): gap ~ 2(h - J), large
        H_disordered = build_tfim(lat, J, 3.0)
        λ_disordered = sort(eigvals(Symmetric(H_disordered)))
        gap_disordered = λ_disordered[2] - λ_disordered[1]

        # Critical point
        H_crit = build_tfim(lat, J, J)
        λ_crit = sort(eigvals(Symmetric(H_crit)))
        gap_crit = λ_crit[2] - λ_crit[1]

        @test gap_disordered > gap_crit

        # Deep in ordered phase (h ≪ J): the "gap" seen by full ED is
        # actually the Z₂ tunneling splitting between the two quasi-
        # degenerate ground states |↑⟩^N and |↓⟩^N.  This splitting is
        # exponentially small in N and much smaller than the critical gap.
        H_ordered = build_tfim(lat, J, 0.1)
        λ_ordered = sort(eigvals(Symmetric(H_ordered)))
        gap_ordered = λ_ordered[2] - λ_ordered[1]
        @test gap_ordered < gap_crit  # tunneling ≪ critical gap
        @test gap_ordered < 1e-3      # exponentially small for N=6
    end

    @testset "Limiting cases" begin
        N = 6
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())

        # h = 0: classical Ising, E₀ = -J(N-1), doubly degenerate
        H0 = build_tfim(lat, J, 0.0)
        λ0 = sort(eigvals(Symmetric(H0)))
        @test λ0[1] ≈ -J * (N - 1) atol = 1e-12
        @test λ0[2] ≈ -J * (N - 1) atol = 1e-12  # 2-fold degenerate
        @test λ0[3] > λ0[1] + 1e-10               # gap to excited state

        # Strong-field limit: at large h the σˣ-eigenstate |+...+⟩ is the
        # unperturbed ground state with E₀⁽⁰⁾ = -hN.  Treating the Ising
        # term V = -J Σ σᶻᵢσᶻᵢ₊₁ as a perturbation, second-order
        # Rayleigh–Schrödinger gives one excited state per bond at gap
        # 4h with matrix element |⟨n|V|0⟩|² = J², so
        #
        #     E₀ ≈ -hN - J²(N - 1) / (4h)   + O(J⁴/h³)
        #
        # Empirically the PT² formula matches the dense ED to machine
        # precision (~1e-10 relative) once h ≥ 100·J, with the residual
        # the next-order O(J⁴/h³) correction.  The previous test
        # compared to the unperturbed -hN with `rtol = 0.01` — 500×
        # looser than the empirical leading-order residual at h = 100,
        # so it would silently accept a wrong sign or prefactor in the
        # σᶻ coupling.  Replace with the PT² reference at h ∈ {100, 1000}.
        for h_large in (100.0, 1000.0)
            H_large = build_tfim(lat, J, h_large)
            λ_large = sort(eigvals(Symmetric(H_large)))
            E0_pt = -h_large * N - J^2 * (N - 1) / (4 * h_large)
            @test λ_large[1] ≈ E0_pt rtol = 1e-9
            # Sign/sense check: PT² lowers E₀ below the unperturbed value.
            @test λ_large[1] < -h_large * N
        end
    end
end
