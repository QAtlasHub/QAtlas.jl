# ─────────────────────────────────────────────────────────────────────────────
# test/util/hubbard_ed.jl
#
# Dense exact-diagonalization free energy for the 1D Hubbard chain, used as the
# INDEPENDENT oracle for the JKS QTM/NLIE row.  It shares no machinery with the
# NLIE — different representation, different algorithm — which is the whole
# point: a second expression of the same integral equations could carry the same
# convention error and agree with it (see the Kaufman/transfer-matrix lesson in
# the IsingSquare stream).
#
# Fock space is the full 4^N (2N spin-orbitals as a 2N-bit integer, Jordan-Wigner
# ordering: orbital 2(i-1)+sigma, sigma = 0 up / 1 down), so this is grand
# canonical at chemical potential `mu` with no particle-number projection —
# exactly what `f = -T log Z / N` needs.
#
#   N = 4 -> 256 states, N = 5 -> 1024, N = 6 -> 4096 (a few seconds).
#
# Dependencies (expected to be `using`'d by the including test file):
#   LinearAlgebra  — `Hermitian`, `eigvals`
#   SparseArrays   — `sparse`
# ─────────────────────────────────────────────────────────────────────────────

"""
    _ed_hubbard_free_energy(N, t, U, mu, beta; pbc=false) -> Float64

Per-site grand-canonical free energy `f = -log(Z)/(beta N)` of the `N`-site 1D
Hubbard chain by dense ED of the full `4^N` Fock space.

`pbc=true` closes the ring, which reaches the thermodynamic limit faster at a
given `N`; comparing the two is how a finite-size effect is told apart from a
disagreement with the thermodynamic limit itself.

The ground-state energy is factored out before exponentiating, so the sum stays
in range at large `beta`.
"""
function _ed_hubbard_free_energy(
    N::Int, t::Real, U::Real, mu::Real, beta::Real; pbc::Bool=false
)
    N >= 2 || throw(ArgumentError("need at least 2 sites; got N = $N"))
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    n_modes = 2 * N
    dim = 2^n_modes

    # Diagonal: on-site repulsion and chemical potential.
    H_diag = zeros(Float64, dim)
    for s in 0:(dim - 1)
        e = 0.0
        for i in 1:N
            su, sd = 2 * (i - 1), 2 * (i - 1) + 1
            nu = (s >> su) & 1
            nd = (s >> sd) & 1
            e += U * nu * nd - mu * (nu + nd)
        end
        H_diag[s + 1] = e
    end

    # Hopping, with the Jordan-Wigner sign from the occupations passed over.
    rows, cols, vals = Int[], Int[], Float64[]
    bonds = pbc ? collect(1:N) : collect(1:(N - 1))
    for i in bonds
        ip = pbc ? (i % N) + 1 : i + 1
        for sigma in 0:1
            si, sip = 2 * (i - 1) + sigma, 2 * (ip - 1) + sigma
            for s in 0:(dim - 1)
                ((s >> sip) & 1 == 1 && (s >> si) & 1 == 0) || continue
                c1 = sum((s >> k) & 1 for k in 0:(sip - 1); init=0)
                s2 = s & ~(1 << sip)
                c2 = sum((s2 >> k) & 1 for k in 0:(si - 1); init=0)
                sg = (iseven(c1) ? 1 : -1) * (iseven(c2) ? 1 : -1)
                s3 = s2 | (1 << si)
                push!(rows, s3 + 1)
                push!(cols, s + 1)
                push!(vals, -t * sg)
                push!(rows, s + 1)
                push!(cols, s3 + 1)
                push!(vals, -t * sg)
            end
        end
    end

    H = Matrix(sparse(rows, cols, vals, dim, dim))
    for i in 1:dim
        H[i, i] += H_diag[i]
    end
    eigs = eigvals(Hermitian(H))
    emin = minimum(eigs)
    return (-log(sum(exp.(-beta .* (eigs .- emin)))) / beta + emin) / N
end
