# test/util/generic_ed.jl — model-agnostic spin-S dense ED.
#
# Independent-route building blocks for `verify(...)` cards.  These
# reconstruct the physics from a *Hamiltonian spec the test supplies
# from the model definition* — they deliberately do NOT call any
# `QAtlas._internal` builder, so a bug in a src matrix-builder cannot
# also corrupt the independent cross-check (true black-box independence).

using LinearAlgebra: Hermitian, Diagonal, eigen, eigvals, I, kron

_eyed(n) = Matrix{ComplexF64}(I, n, n)

"""
    spin_ops(S) -> (Sx, Sy, Sz)

Dense `(2S+1)×(2S+1)` spin-`S` operators (ħ = 1).  S ∈ {1/2, 1, 3/2, …}.
"""
function spin_ops(S::Real)
    d = Int(round(2S + 1))
    m = [S - (k - 1) for k in 1:d]                     # S, S-1, …, -S
    Sz = Matrix{ComplexF64}(Diagonal(m))
    Sp = zeros(ComplexF64, d, d)
    for k in 1:(d - 1)
        mk = m[k + 1]                                  # lower state's m
        Sp[k, k + 1] = sqrt(S * (S + 1) - mk * (mk + 1))
    end
    Sm = collect(Sp')
    Sx = (Sp + Sm) / 2
    Sy = (Sp - Sm) / (2im)
    return Sx, Sy, Sz
end

"""
    chain_hamiltonian(d, N, bond; onsite=nothing) -> Hermitian

Dense `dᴺ×dᴺ` OBC Hamiltonian from a `d²×d²` nearest-neighbour `bond`
operator (+ optional `d×d` `onsite`).  The caller supplies `bond` from
the model's *written Hamiltonian*, not from src.
"""
function chain_hamiltonian(d::Int, N::Int, bond::AbstractMatrix; onsite=nothing)
    D = d^N
    H = zeros(ComplexF64, D, D)
    for i in 1:(N - 1)
        H .+= kron(_eyed(d^(i - 1)), bond, _eyed(d^(N - i - 1)))
    end
    if onsite !== nothing
        for i in 1:N
            H .+= kron(_eyed(d^(i - 1)), onsite, _eyed(d^(N - i)))
        end
    end
    return Hermitian(H)
end

"""
    chain_hamiltonian_pbc(d, N, terms) -> Hermitian

PBC chain where the bond is `Σ A⊗B` over `terms::Vector{Tuple}` (e.g.
`[(Sx,Sx),(Sy,Sy),(Sz,Sz)]` for S·S).  Adds NN bonds i=1..N-1 plus the
wrap (site N, site 1).  Exact for sum-of-products bonds (the
AKLT/Heisenberg family).
"""
function chain_hamiltonian_pbc(d::Int, N::Int, terms::Vector{<:Tuple})
    D = d^N
    H = zeros(ComplexF64, D, D)
    bond = sum(kron(A, B) for (A, B) in terms)
    for i in 1:(N - 1)
        H .+= kron(_eyed(d^(i - 1)), bond, _eyed(d^(N - i - 1)))
    end
    for (A, B) in terms                                # wrap (N,1)
        opN = kron(_eyed(d^(N - 1)), A)
        op1 = kron(B, _eyed(d^(N - 1)))
        H .+= opN * op1
    end
    return Hermitian(H)
end

dense_spectrum(H) = sort(real.(eigvals(Hermitian(Matrix(H)))))

function ground_state(H)
    F = eigen(Hermitian(Matrix(H)))
    return F.values[1], F.vectors[:, 1]
end

"""
    site_op(o, d, N, site) -> Matrix

Embed a single-site `d×d` operator at `site` in an `N`-site chain.
"""
function site_op(o::AbstractMatrix, d::Int, N::Int, site::Int)
    return kron(_eyed(d^(site - 1)), o, _eyed(d^(N - site)))
end

"""
    two_point(ψ, d, N, o, i, j) -> Float64

`⟨ψ| oᵢ oⱼ |ψ⟩` (real part) for a single-site operator `o`.
"""
function two_point(ψ, d::Int, N::Int, o::AbstractMatrix, i::Int, j::Int)
    Oi = site_op(o, d, N, i)
    Oj = site_op(o, d, N, j)
    return real(ψ' * (Oi * (Oj * ψ)))
end

"""
    thermo_from_spectrum(evals, β) -> (E, F, S, C)  [totals]

Exact canonical thermodynamics from an eigenvalue list (log-sum-exp
shifted).  Per-site = divide by N at the call site.
"""
function thermo_from_spectrum(evals::AbstractVector, β::Real)
    emin = minimum(evals)
    w = exp.(-β .* (evals .- emin))
    Z = sum(w)
    E = sum(evals .* w) / Z
    E2 = sum((evals .^ 2) .* w) / Z
    F = -(log(Z) - β * emin) / β
    Sent = β * (E - F)
    Cv = β^2 * (E2 - E^2)
    return E, F, Sent, Cv
end
