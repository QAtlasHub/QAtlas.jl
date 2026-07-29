# ─────────────────────────────────────────────────────────────────────────────
# test/util/aklt_dense_ed.jl
#
# Many-body dense ED for the OBC AKLT chain — a TEST INSTRUMENT, not an atlas
# implementation.
#
# WHY IT LIVES HERE.  QAtlas serves reference values to other packages, so what
# it registers has to be something those packages can afford to ask for.  Dense
# ED on a `3^N` Hilbert space is not: it is exponential and capped at N ≤ 8,
# which makes it a poor oracle however correct it is.  The atlas therefore
# advertises only routes that scale — closed forms, Bethe ansatz / NLIE, and
# free-fermion (BdG, single-particle) diagonalisation.
#
# ED keeps its real job, which is the opposite one: it is how the ANALYTIC claims
# get checked.  `test_aklt_structural.jl` uses this to confirm that
# E₀(N) = -(2/3)(N-1) exactly, that the OBC ground state is 4-fold degenerate
# (the AKLT edge-mode theorem), and that the spectrum has the right shape —
# statements the atlas asserts and ED independently confirms at small N.
#
# Moved verbatim from `src/models/quantum/AKLT/AKLT1D.jl`.
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra: Hermitian, I, eigvals, kron
using QAtlas: AKLT1D, OBC, _MAX_ED_SITES_S1, _S1_x, _S1_y, _S1_z

"""
    _aklt_dense_hamiltonian(model::AKLT1D, N::Int) -> Matrix{ComplexF64}

Dense `3^N × 3^N` OBC AKLT Hamiltonian

    H = J Σᵢ [ Sᵢ · Sᵢ₊₁ + (1/3) (Sᵢ · Sᵢ₊₁)² ]

built from the spin-1 primitives in `HeisenbergS1.jl`.  Capped by
`_MAX_ED_SITES_S1`.
"""
function _aklt_dense_hamiltonian(model::AKLT1D, N::Int)
    N ≥ 2 || throw(ArgumentError("AKLT1D OBC chain needs N ≥ 2 (got N = $N)"))
    N ≤ _MAX_ED_SITES_S1 || throw(
        ArgumentError("spin-1 dense ED is capped at N ≤ $(_MAX_ED_SITES_S1) (got N = $N)"),
    )
    J = model.J
    D = 3^N
    # 9×9 single-bond block: S₁·S₂ + (1/3)(S₁·S₂)²
    SdotS = kron(_S1_x, _S1_x) + kron(_S1_y, _S1_y) + kron(_S1_z, _S1_z)
    bond = J * (SdotS + (1.0 / 3.0) * (SdotS * SdotS))
    H = zeros(ComplexF64, D, D)
    for i in 1:(N - 1)
        d_left = 3^(i - 1)
        d_right = 3^(N - i - 1)
        H .+= kron(
            Matrix{ComplexF64}(I, d_left, d_left),
            bond,
            Matrix{ComplexF64}(I, d_right, d_right),
        )
    end
    return H
end

"""
    _aklt_dense_spectrum(model::AKLT1D, N::Int) -> Vector{Float64}

Sorted full eigenvalue spectrum of the OBC AKLT chain on `N` sites.
"""
function _aklt_dense_spectrum(model::AKLT1D, N::Int)
    return sort(real.(eigvals(Hermitian(_aklt_dense_hamiltonian(model, N)))))
end
