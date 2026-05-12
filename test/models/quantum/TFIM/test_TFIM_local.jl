# =============================================================================
# Tests for TFIM site-local thermal observables (TFIM_local.jl).
#
# Cross-checks:
#   1. Consistency with existing scalar quantities:
#        mean(magnetization_x_local)  ≈ transverse_magnetization
#        sum(energy_local)            ≈ energy        (exact per-site sum)
#        magnetization_z_local        ≡ 0
#   2. Consistency of the nearest-neighbour bond extraction against the full
#      Pfaffian correlator:
#        energy_local[i] recomputed from :zz_static_thermal bonds agrees.
#   3. Direct comparison against dense ED at small N for per-site values.
#   4. High-temperature limit: every local entry → 0 as β → 0.
# =============================================================================

@testset "TFIM site-local observables (OBC)" begin

    # ───────────────────────────────────────────────────────────────────────
    # Layer 1: consistency with existing scalar entries.
    # Use disordered phase (h > J) at moderate β so OBC edge modes do not
    # complicate things.
    # ───────────────────────────────────────────────────────────────────────
    @testset "consistency with scalar entries" begin
        for (N, J, h, β) in [
            (10, 1.0, 1.5, 0.8),
            (12, 1.0, 1.0, 1.5),   # critical point
            (8, 1.0, 0.5, 2.0),    # ordered phase, β below edge-mode scale
        ]
            mx_loc = QAtlas.fetch(
                TFIM(; J=J, h=h), MagnetizationXLocal(), OBC(; N=N); beta=β
            )
            @test mx_loc isa Vector{Float64}
            @test length(mx_loc) == N

            mx_avg_scalar = QAtlas.fetch(
                TFIM(; J=J, h=h), MagnetizationX(), OBC(; N=N); beta=β
            )
            @test sum(mx_loc) / N ≈ mx_avg_scalar atol=1e-12

            mz_loc = QAtlas.fetch(
                TFIM(; J=J, h=h), MagnetizationZLocal(), OBC(; N=N); beta=β
            )
            @test mz_loc == zeros(Float64, N)

            ε_loc = QAtlas.fetch(TFIM(; J=J, h=h), EnergyLocal(), OBC(; N=N); beta=β)
            @test ε_loc isa Vector{Float64}
            @test length(ε_loc) == N

            E_total = QAtlas.fetch(TFIM(; J=J, h=h), Energy(), OBC(; N=N); beta=β)
            @test sum(ε_loc) ≈ E_total atol=1e-10
        end
    end

    # ───────────────────────────────────────────────────────────────────────
    # Layer 2: the fast covariance-based bond extraction agrees with the
    # Pfaffian-based full correlator matrix.
    # ───────────────────────────────────────────────────────────────────────
    @testset "bond extraction vs Pfaffian correlator" begin
        N, J, h, β = 8, 1.0, 1.2, 1.0
        ε_loc = QAtlas.fetch(TFIM(; J=J, h=h), EnergyLocal(), OBC(; N=N); beta=β)
        mx_loc = QAtlas.fetch(TFIM(; J=J, h=h), MagnetizationXLocal(), OBC(; N=N); beta=β)
        C = QAtlas.fetch(TFIM(; J=J, h=h), ZZCorrelation{:static}(), OBC(; N=N); beta=β)

        # ε_i reconstructed from the Pfaffian nearest-neighbour correlator.
        ε_ref = Vector{Float64}(undef, N)
        for i in 1:N
            left = i > 1 ? C[i - 1, i] : 0.0
            right = i < N ? C[i, i + 1] : 0.0
            ε_ref[i] = -(J / 2) * (left + right) - h * mx_loc[i]
        end
        @test ε_loc ≈ ε_ref atol=1e-10
    end

    # ───────────────────────────────────────────────────────────────────────
    # Layer 3: dense ED comparison at small N for each site.
    # ───────────────────────────────────────────────────────────────────────
    @testset "per-site ED comparison (N = 4)" begin
        N, J, h, β = 4, 1.0, 0.8, 1.3
        mx_loc = QAtlas.fetch(TFIM(; J=J, h=h), MagnetizationXLocal(), OBC(; N=N); beta=β)
        ε_loc = QAtlas.fetch(TFIM(; J=J, h=h), EnergyLocal(), OBC(; N=N); beta=β)

        H = _build_tfim_dense(N, J, h)
        E, V = eigen(Matrix(H))
        w = exp.(-β .* (E .- E[1]))
        Z = sum(w)
        w ./= Z

        # ⟨σˣ_i⟩ per site
        for i in 1:N
            Sxi = _op_site(_SX, i, N)
            mx_ed = real(sum(w .* diag(V' * Sxi * V)))
            @test mx_loc[i] ≈ mx_ed atol=1e-10
        end

        # Local energy: bonds split symmetrically (rebuild operator by hand).
        for i in 1:N
            op = -h * _op_site(_SX, i, N)
            if i > 1
                op -= (J / 2) * _op_site(_SZ, i - 1, N) * _op_site(_SZ, i, N)
            end
            if i < N
                op -= (J / 2) * _op_site(_SZ, i, N) * _op_site(_SZ, i + 1, N)
            end
            ε_ed = real(sum(w .* diag(V' * op * V)))
            @test ε_loc[i] ≈ ε_ed atol=1e-10
        end
    end

    # ───────────────────────────────────────────────────────────────────────
    # Layer 4: β → 0 high-temperature limit — every local quantity vanishes.
    # ───────────────────────────────────────────────────────────────────────
    @testset "β → 0 limit" begin
        N, J, h, β = 10, 1.0, 1.0, 1e-8
        mx_loc = QAtlas.fetch(TFIM(; J=J, h=h), MagnetizationXLocal(), OBC(; N=N); beta=β)
        ε_loc = QAtlas.fetch(TFIM(; J=J, h=h), EnergyLocal(), OBC(; N=N); beta=β)
        @test maximum(abs, mx_loc) < 1e-6
        @test maximum(abs, ε_loc) < 1e-6
    end

    # ───────────────────────────────────────────────────────────────────────
    # Layer 5: boundary vs bulk — OBC breaks translation invariance, so the
    # local profile is non-uniform while the bulk value converges to the
    # Infinite result.  We just sanity-check non-uniformity and symmetry.
    # ───────────────────────────────────────────────────────────────────────
    @testset "profile symmetry" begin
        N, J, h, β = 16, 1.0, 1.5, 1.5
        mx_loc = QAtlas.fetch(TFIM(; J=J, h=h), MagnetizationXLocal(), OBC(; N=N); beta=β)
        ε_loc = QAtlas.fetch(TFIM(; J=J, h=h), EnergyLocal(), OBC(; N=N); beta=β)

        # Reflection symmetry about the centre (OBC chain is invariant under
        # i ↔ N + 1 - i).
        @test mx_loc ≈ reverse(mx_loc) atol=1e-10
        @test ε_loc ≈ reverse(ε_loc) atol=1e-10

        # Non-uniformity: bulk ≠ boundary (profile should actually depend on
        # site — if this test passes trivially we've lost the information
        # we built the local quantities for).
        @test maximum(mx_loc) - minimum(mx_loc) > 1e-4
    end
end
