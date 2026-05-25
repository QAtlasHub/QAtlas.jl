# API Reference

Complete Julia docstring index for QAtlas.jl, organized by category.
For narrative documentation see the [Models](models/index.md),
[Universality Classes](universalities/index.md), and [Methods](methods/index.md) pages.

---

## Core API

```@docs
QAtlas.fetch
```

### Boundary Conditions

```@docs
BoundaryCondition
OBC
PBC
Infinite
```

### Types

```@docs
AbstractModel
AbstractQuantity
```

---

## Quantities

```@autodocs
Modules = [QAtlas]
Filter = t -> t isa Type && t <: QAtlas.AbstractQuantity
```

---

## Models — Classical

```@autodocs
Modules = [QAtlas]
Filter = t -> begin
    !(t isa Type) && return false
    t <: QAtlas.AbstractModel || return false
    m = parentmodule(t)
    occursin("classical", string(m)) || occursin("Classical", string(m))
end
```

---

## Models — Quantum

```@autodocs
Modules = [QAtlas]
Filter = t -> begin
    !(t isa Type) && return false
    t <: QAtlas.AbstractModel || return false
    m = parentmodule(t)
    occursin("quantum", string(m)) || occursin("Quantum", string(m))
end
```

---

## Universality Classes

```@docs
Universality
CriticalExponents
GrowthExponents
MeanField
MinimalModel
WZWSU2
```

---

## Full Index

All exported symbols not covered above.

```@autodocs
Modules = [QAtlas]
```
