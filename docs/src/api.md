# API Primer

This page covers the **universal interface** of QAtlas.jl — the parts that
apply regardless of which model or quantity you are working with.
For model-specific documentation (Hamiltonian, parameters, stored quantities,
verification cards) see the individual model pages under
[Models](models/index.md).

---

## `fetch` — the primary entry point

```julia
QAtlas.fetch(model, quantity, bc; kwargs...) -> Number | NamedTuple
```

| Argument | Type | Description |
|----------|------|-------------|
| `model` | `<: AbstractModel` | The physical model (e.g. `TFIM(J=1.0, h=0.5)`) |
| `quantity` | `<: AbstractQuantity` | What to compute (e.g. `Energy(:per_site)`) |
| `bc` | `<: BoundaryCondition` | System size / topology |
| `kwargs` | varies | Model-specific parameters (e.g. `β`, `N`) |

The return type depends on the quantity:
- Scalar quantities (`Energy`, `MassGap`, …) → `Float64` or `Rational`
- Multi-component quantities (`CriticalExponents`, `AnyonStatistics`, …) → `NamedTuple`

```julia
using QAtlas

# Scalar
E = QAtlas.fetch(TFIM(J=1.0, h=0.5), Energy(:per_site), OBC(16))

# NamedTuple
e = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2)
e.β    # 1//8
e.ν    # 1//1
```

---

## Boundary Conditions

```julia
OBC(N)        # open chain / slab of N sites
PBC(N)        # periodic ring of N sites
Infinite()    # thermodynamic limit (k-space / Bethe ansatz)
```

All three share the supertype `BoundaryCondition`.  Not all models support
all boundary conditions — see the Quantity × BC matrix on each model page.

---

## Quantity types

Quantities are parameterised types whose type parameter encodes the
*variant* of the quantity:

```julia
Energy(:total)       # total ground-state energy E₀
Energy(:per_site)    # E₀ / N

FreeEnergy()         # F = -T ln Z (or -β⁻¹ ln Z)
ThermalEntropy()     # S = -∂F/∂T
SpecificHeat()       # C = -T ∂²F/∂T²

MassGap()            # many-body gap Δ = E₁ - E₀
CorrelationLength()  # ξ from exponential decay
CentralCharge()      # CFT central charge c

VonNeumannEntropy()  # S_vN = -tr(ρ ln ρ)
RenyiEntropy(α)      # S_α = (1-α)⁻¹ ln tr(ρ^α)

CriticalExponents()  # returns NamedTuple (β, ν, γ, η, δ, α, c)
```

All quantity types are subtypes of `AbstractQuantity`.

---

## Model constructors

Every model is a Julia struct with keyword-argument fields for its
physical parameters.  Default values are always provided.

```julia
# Quantum spin models
TFIM(; J=1.0, h=1.0)
Heisenberg1D(; J=1.0)
XXZ1D(; J=1.0, Δ=1.0)
Hubbard1D(; t=1.0, U=4.0)

# Classical lattice models
IsingSquare(; J=1.0)
IsingTriangular(; J=1.0)

# Universality classes
Universality(:Ising)        # equivalent to Universality{:Ising}()
MeanField()
MinimalModel(3, 4)          # M(p, p') minimal model
WZWSU2(k)                   # SU(2)_k WZW model
```

---

## Checking what is implemented

```julia
QAtlas.implementation_status(TFIM())
# returns a Markdown table of all registered (quantity, bc) pairs

QAtlas.implementation_status_markdown(TFIM())
# same, as a String
```

---

## Operator conventions

QAtlas standardises on the **spin operator** convention (`S^α` with
eigenvalues in `{-S, …, +S}`) for all spin-observable return values.
For fermionic and topological models, see the [Conventions](conventions.md)
page.

---

## See also

- [Models](models/index.md) — per-model Hamiltonian, parameters, quantity matrix
- [Universality Classes](universalities/index.md) — `Universality{C}` API
- [Conventions](conventions.md) — operator normalisation rules
- [Verification](verification/index.md) — how correctness is ensured

---

## Full API Index

```@autodocs
Modules = [QAtlas, QAtlas.XXZKlumperNLIE]
```
