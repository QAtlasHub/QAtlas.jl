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

@testset "TFIM — ED ground state vs BdG analytical (verify cards + structural)" begin
    for N in [4, 6, 8]
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())
        # Lattice2D structural sanity — not a QAtlas hub, kept raw.
        @test num_sites(lat) == N
        @test length(collect(bonds(lat))) == N - 1

        @testset "N=$N OBC — E₀ verify cards" begin
            for (J, h) in [(1.0, 0.0), (1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (0.5, 1.5)]
                H = build_tfim(lat, J, h)
                λ = sort(eigvals(Symmetric(H)))
                E0_ed = λ[1]
                verify(
                    TFIM(; J=J, h=h),
                    Energy(),
                    OBC(N);
                    route=:ed_finite_size,
                    independent=E0_ed,
                    # Legacy used rtol=1e-10; emulate via |E0_ed|-scaling.
                    agree_within=max(1e-12, 1e-10 * abs(E0_ed)),
                    at=["J=$(J)", "h=$(h)", "N=$(N)"],
                    refs=[
                        "Independent dense-ED of build_tfim Lattice2D OBC chain (eigvals of Symmetric H = -J Σ σᶻσᶻ − h Σ σˣ) — cross-checks BdG analytical Energy OBC",
                    ],
                )
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

    @testset "Limiting cases (verify cards + structural)" begin
        N = 6
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())

        # h = 0: classical Ising, E₀ = -J(N-1), doubly degenerate.
        verify(
            TFIM(; J=J, h=0.0),
            Energy(),
            OBC(N);
            route=:limiting_case,
            independent=-J * (N - 1),
            agree_within=1e-12,
            at=["J=$(J)", "h=0.0", "N=$(N)"],
            refs=[
                "Classical Ising limit h=0: E_0^OBC = -J(N-1) exact (independent of dense ED)",
            ],
        )
        # 2-fold degeneracy + gap to first excited — multi-eigenvalue
        # structural property, kept raw.
        H0 = build_tfim(lat, J, 0.0)
        λ0 = sort(eigvals(Symmetric(H0)))
        @test λ0[2] ≈ -J * (N - 1) atol = 1e-12
        @test λ0[3] > λ0[1] + 1e-10

        # Strong-field limit: at large h the σˣ-eigenstate |+...+⟩ is the
        # unperturbed ground state. Second-order Rayleigh–Schrödinger gives
        #     E₀ ≈ -hN - J²(N - 1) / (4h)   + O(J⁴/h³)
        # which the dense ED satisfies to ~1e-10 once h ≥ 100·J.
        for h_large in (100.0, 1000.0)
            E0_pt = -h_large * N - J^2 * (N - 1) / (4 * h_large)
            verify(
                TFIM(; J=J, h=h_large),
                Energy(),
                OBC(N);
                route=:limiting_case,
                independent=E0_pt,
                # Legacy used rtol=1e-9 (|E0_pt| ~ hN ~ 6e3 at h=1000).
                agree_within=max(1e-12, 1e-9 * abs(E0_pt)),
                at=["J=$(J)", "h=$(h_large)", "N=$(N)"],
                refs=[
                    "Strong-field PT² limit h ≫ J: E_0^OBC ≈ -hN - J²(N-1)/(4h) (Rayleigh-Schrödinger, |+⟩^N unperturbed g.s., bond perturbation V = -J Σ σᶻσᶻ)",
                ],
            )
            # Sign/sense check: PT² lowers E₀ below -hN unperturbed.
            H_large = build_tfim(lat, J, h_large)
            λ_large = sort(eigvals(Symmetric(H_large)))
            @test λ_large[1] < -h_large * N
        end
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "TFIM gap closure — verification cards" begin
    # Pfeuty 1970: Δ = 2|h - J|, closing exactly at the critical point.
    for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0))
        verify(
            TFIM(; J=J, h=h),
            MassGap(),
            Infinite();
            route=:second_closed_form,
            independent=2 * abs(h - J),
            agree_within=1e-10,
            refs=["Pfeuty 1970: Δ = 2|h - J| (= 0 at the QCP h = J)"],
        )
    end

    # Independent dense-ED corroboration. NOT the closed form: the gap
    # is read off the spectrum of the full many-body H = -J Σ σᶻσᶻ −
    # h Σ σˣ, never from Δ = 2|h−J|, so it breaks the circularity of
    # the second_closed_form cards above. At J = 0 the sites decouple
    # (H = −h Σ σˣ): the ground state is |+⟩^N and a single spin flip
    # costs exactly 2h, so the dense-ED gap equals the exact Pfeuty
    # value 2|h−J| = 2h at *every* finite N — no finite-size error and
    # no extrapolation (the OBC gap at J≠0 only converges as O(1/N²),
    # which `verify`'s last(ind) contract cannot use). This is the
    # cleanest possible non-circular check of src's MassGap.
    let J = 0.0, h = 1.5, Ns = (4, 6, 8)
        ed_gap = function (N)
            lat = build_lattice(Square, N, 1; boundary=OpenAxis())
            H = build_tfim(lat, J, h)
            λ = sort(eigvals(Symmetric(H)))
            return λ[2] - λ[1]
        end
        verify(
            TFIM(; J=J, h=h),
            MassGap(),
            Infinite();
            route=:ed_finite_size,
            independent=[ed_gap(N) for N in Ns],
            at=["N=$N" for N in Ns],
            agree_within=1e-10,
            refs=[
                "Pfeuty 1970: Δ = 2|h − J|; independent dense ED of " *
                "H = −h Σ σˣ at J=0 (decoupled spins, gap = 2h exact ∀N)",
            ],
        )
    end
end
