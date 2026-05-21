# =============================================================================
# KitaevHoneycomb - sector-enumerated formula vs sparse ED (testsets 8-10, heaviest)
#
# Split out of test/models/quantum/KitaevHoneycomb/test_KitaevHoneycomb.jl
# (2.5 min on s06) so each top-level group runs on its own shard.
# =============================================================================

using QAtlas: KitaevHoneycomb, Energy, MassGap, Infinite, PBC, OBC
using Test
using LinearAlgebra: eigvals, Hermitian, I as I_mat
using SparseArrays: sparse, spzeros
using KrylovKit: eigsolve
using Lattice2D: build_lattice, Honeycomb, OpenAxis, PeriodicAxis
using LatticeCore: num_sites, bonds

# Reference value taken from Baskaran, Mandal & Shankar (PRL 98, 247201)
# and reproduced in Feng et al. DMRG benchmarks of the Kitaev model:
# at the isotropic point (Kx = Ky = Kz = 1) the ground state energy
# per site for H = -Σ Kγ σγ σγ is
#
#     ε₀ ≈ -0.7872986...
#
# i.e. `-⟨|f(k)|⟩ / 2` in the conventions of the module.

const ε_isotropic_TL = -0.7872986216706852

"""
    _kitaev_ed_gs_per_site(lat, Kx, Ky, Kz) -> Float64

Dense ED of the Kitaev Hamiltonian on the given Lattice2D lattice.
Iterates `bonds(lat)` and emits σˣσˣ / σʸσʸ / σᶻσᶻ on bond.type ∈
(`:type_2`, `:type_3`, `:type_1`) respectively, matching the QAtlas
convention.
"""
function _kitaev_ed_gs_per_site(lat, Kx::Real, Ky::Real, Kz::Real)
    N = num_sites(lat)
    dim = 2^N
    σx = ComplexF64[0 1; 1 0]
    σy = ComplexF64[0 -im; im 0]
    σz = ComplexF64[1 0; 0 -1]
    function embed2(A, B, i, j)
        @assert i != j
        i1, j1, A1, B1 = i < j ? (i, j, A, B) : (j, i, B, A)
        L = ComplexF64.(I_mat(2^(i1 - 1)))
        M = ComplexF64.(I_mat(2^(j1 - i1 - 1)))
        R = ComplexF64.(I_mat(2^(N - j1)))
        return kron(kron(kron(kron(L, A1), M), B1), R)
    end
    H = zeros(ComplexF64, dim, dim)
    for b in bonds(lat)
        if b.type === :type_1
            H .-= Kz * embed2(σz, σz, b.i, b.j)
        elseif b.type === :type_2
            H .-= Kx * embed2(σx, σx, b.i, b.j)
        elseif b.type === :type_3
            H .-= Ky * embed2(σy, σy, b.i, b.j)
        else
            error("unexpected Kitaev bond type: $(b.type)")
        end
    end
    E_gs = real(eigvals(Hermitian(H))[1])
    return E_gs / N
end

"""
    _kitaev_ed_gs_per_site_sparse(lat, Kx, Ky, Kz) -> Float64

Sparse-matrix Kitaev ED via Lanczos (`eigsolve` with
`ishermitian=true`). Needed once the Hilbert space exceeds the
~2^16 dense limit; the Kitaev Hamiltonian is complex Hermitian
(σʸ makes it so) so the existing `ground_state_krylov` — which
assumes real symmetric — doesn't apply directly.
"""
function _kitaev_ed_gs_per_site_sparse(lat, Kx::Real, Ky::Real, Kz::Real; seed::Int=0)
    σx = sparse(ComplexF64[0 1; 1 0])
    σy = sparse(ComplexF64[0 -im; im 0])
    σz = sparse(ComplexF64[1 0; 0 -1])
    function embed2(A, B, i, j, N)
        i1, j1, A1, B1 = i < j ? (i, j, A, B) : (j, i, B, A)
        L = sparse(ComplexF64.(I_mat(2^(i1 - 1))))
        M = sparse(ComplexF64.(I_mat(2^(j1 - i1 - 1))))
        R = sparse(ComplexF64.(I_mat(2^(N - j1))))
        return kron(kron(kron(kron(L, A1), M), B1), R)
    end
    N = num_sites(lat)
    D = 2^N
    H = spzeros(ComplexF64, D, D)
    for b in bonds(lat)
        if b.type === :type_1
            H -= Kz * embed2(σz, σz, b.i, b.j, N)
        elseif b.type === :type_2
            H -= Kx * embed2(σx, σx, b.i, b.j, N)
        elseif b.type === :type_3
            H -= Ky * embed2(σy, σy, b.i, b.j, N)
        end
    end
    # Deterministic start vector so tests are reproducible.
    Random.seed!(seed)
    x0 = randn(ComplexF64, D)
    vals, _, info = eigsolve(H, x0, 1, :SR; ishermitian=true, tol=1e-10, krylovdim=30)
    info.converged < 1 && @warn "Kitaev sparse ED failed to converge" info
    return real(vals[1]) / N
end

@testset "KitaevHoneycomb PBC: sector-enumerated formula matches sparse ED at 3×3" begin
    # Lx = Ly = 3 is the smallest symmetric torus where the flux-free
    # sector formula reproduces the true spin-Hamiltonian ground
    # state exactly (the 2×2 isotropic case sits in a vortex-bearing
    # sector and is skipped — see the separate smaller-torus testset).
    # 18 sites = 2^18 = 262144-dim Hilbert space, sparse Lanczos
    # runs in a few seconds.
    lat = build_lattice(Honeycomb, 3, 3; boundary=PeriodicAxis())
    @test num_sites(lat) == 18
    # Isotropic & generic anisotropic B-phase points: exact agreement
    # (residual at machine precision).
    for (Kx, Ky, Kz) in [(1.0, 1.0, 1.0), (0.3, 0.7, 1.0), (0.8, 1.0, 1.1)]
        m = KitaevHoneycomb(; Kx=Kx, Ky=Ky, Kz=Kz)
        E_ed = _kitaev_ed_gs_per_site_sparse(lat, Kx, Ky, Kz)
        E_fl = QAtlas.fetch(m, Energy(), PBC(0); Lx=3, Ly=3)
        @info "PBC 3×3 formula vs sparse ED (B-phase)" Kx Ky Kz E_ed E_fl Δ = E_ed - E_fl
        @test abs(E_fl - E_ed) < 1e-8
    end
    # Gapped A-phase points: flux-free ansatz does not always sit in
    # the true vortex sector on small tori, so agreement relaxes to
    # O(10⁻³). Tight agreement requires either explicit vortex-sector
    # enumeration or larger L.
    for (Kx, Ky, Kz) in [(2.0, 0.5, 0.5), (0.5, 0.5, 2.0)]
        m = KitaevHoneycomb(; Kx=Kx, Ky=Ky, Kz=Kz)
        E_ed = _kitaev_ed_gs_per_site_sparse(lat, Kx, Ky, Kz)
        E_fl = QAtlas.fetch(m, Energy(), PBC(0); Lx=3, Ly=3)
        @info "PBC 3×3 formula vs sparse ED (A-phase)" Kx Ky Kz E_ed E_fl Δ = E_ed - E_fl
        @test abs(E_fl - E_ed) < 5e-3
    end
end

@testset "KitaevHoneycomb PBC: sector-enumerated formula matches ED on small tori" begin
    # The PBC fetch enumerates all four topological flux sectors
    # (νx, νy) ∈ {0, 1/2}² and takes the minimum. On finite tori this
    # is essential; without sector enumeration the formula undershoots
    # ED by up to ~15% at Lx=Ly=2 isotropic.
    #
    # Cross-check: asymmetric (3,2) / (2,3) tori where all four Wilson
    # loops can be distinguished hit machine precision agreement.
    # Anisotropic K's at (2,2) also agree to < 0.5% (residual is the
    # known vortex-sector physics pathology of the smallest isotropic
    # torus).
    cases = [
        ((3, 2), (1.0, 1.0, 1.0), 1e-10),      # exact match
        ((2, 3), (1.0, 1.0, 1.0), 1e-10),      # exact match
        ((3, 2), (0.3, 0.7, 1.0), 1e-10),      # exact match
        ((2, 3), (2.0, 0.5, 0.5), 1e-10),      # exact match (Ax-phase)
        ((2, 2), (2.0, 0.5, 0.5), 5e-3),       # Ax-phase: within 0.5%
    ]
    for ((Lx, Ly), (Kx, Ky, Kz), tol) in cases
        m = KitaevHoneycomb(; Kx=Kx, Ky=Ky, Kz=Kz)
        lat = build_lattice(Honeycomb, Lx, Ly; boundary=PeriodicAxis())
        E_ed = _kitaev_ed_gs_per_site(lat, Kx, Ky, Kz)
        E_fl = QAtlas.fetch(m, Energy(), PBC(0); Lx=Lx, Ly=Ly)
        @info "PBC formula vs ED" Lx Ly Kx Ky Kz E_ed E_fl Δ = E_ed - E_fl
        @test abs(E_fl - E_ed) < tol
    end
end

@testset "KitaevHoneycomb OBC: flux-free formula matches ED on small clusters" begin
    # (a) Single z-bond (Lx = Ly = 1 OBC on Honeycomb reduces to just the
    # same-cell z-bond). H = -Kz σᶻσᶻ on 2 sites ⇒ Egs/site = -Kz/2.
    # Both the formula and ED must land on -0.5 at Kz = 1.
    m = KitaevHoneycomb(; Kx=1.0, Ky=1.0, Kz=1.0)
    lat1 = build_lattice(Honeycomb, 1, 1; boundary=OpenAxis())
    E_ed1 = _kitaev_ed_gs_per_site(lat1, m.Kx, m.Ky, m.Kz)
    E_fl1 = QAtlas.fetch(m, Energy(), OBC(0); Lx=1, Ly=1)
    @test E_ed1 ≈ -0.5 rtol = 1e-12
    @test E_fl1 ≈ E_ed1 rtol = 1e-12

    # (b) Full 2 × 2 OBC: 8 sites, 7 bonds (4 z + 1 x + 2 y per the
    # Lattice2D Honeycomb topology). Non-trivial cross-check that
    # flux-free ansatz matches ED on the real spin Hilbert space.
    lat2 = build_lattice(Honeycomb, 2, 2; boundary=OpenAxis())
    @test num_sites(lat2) == 8
    for (Kx, Ky, Kz) in [
        (1.0, 1.0, 1.0),           # isotropic B-phase
        (2.0, 0.5, 0.5),           # Ax-phase (gapped)
        (0.3, 0.7, 1.0),
    ]           # anisotropic
        m_ed = KitaevHoneycomb(; Kx=Kx, Ky=Ky, Kz=Kz)
        E_ed = _kitaev_ed_gs_per_site(lat2, Kx, Ky, Kz)
        E_fl = QAtlas.fetch(m_ed, Energy(), OBC(0); Lx=2, Ly=2)
        @test E_fl ≈ E_ed rtol = 1e-10
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
