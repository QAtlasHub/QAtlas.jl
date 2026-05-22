# Triangular Lattice

## Overview

The nearest-neighbor tight-binding model on the triangular lattice
is the simplest frustrated lattice hopping problem. Each site has 6
nearest neighbours, the lattice has a single sublattice per unit cell,
and the resulting single band exhibits a characteristic asymmetric
dispersion that is the hallmark of geometric frustration.

$``H = -t \sum_{\langle i,j \rangle} \bigl(c^\dagger_i c_j + c^\dagger_j c_i\bigr)``$

**Lattice properties**: 1 sublattice per unit cell, 6 nearest
neighbours, **not bipartite** (frustrated), no flat band.

**Key physics**: Because the triangular lattice is non-bipartite, there
is no chiral (particle-hole) symmetry, and the spectrum is **not**
symmetric about ``E = 0``. The band ranges from ``-6t`` (at ``\Gamma``) to
``+3t`` (at the ``K``-points), an asymmetric window that is a direct
consequence of frustration. The density of states exhibits a Van Hove
singularity within the band.

---

## Bloch Spectrum

### Statement

Since the triangular lattice has one sublattice, the Bloch
"Hamiltonian" is a scalar at each ``\mathbf{k}``-point:

$``E_{mn} = -2t\left[\cos\theta_1 + \cos\theta_2 + \cos(\theta_2 - \theta_1)\right]``$

where ``\theta_1 = 2\pi m/L_x`` and ``\theta_2 = 2\pi n/L_y``, with
``\mathbf{a}_1 = (1, 0)``, ``\mathbf{a}_2 = (1/2, \sqrt{3}/2)`` the
primitive vectors of `Lattice2D`'s triangular topology.

The full spectrum has ``L_x L_y`` eigenvalues, one per allowed momentum.

### Band edges

| Point | ``(m/L_x, n/L_y)`` | Energy |
|-------|-------------------|--------|
| ``\Gamma`` | ``(0, 0)`` | ``-6t`` (global minimum, unique) |
| ``K`` | ``(1/3, 2/3)`` | ``+3t`` (global maximum) |
| ``K'`` | ``(2/3, 1/3)`` | ``+3t`` (global maximum) |

The band range ``[-6t, +3t]`` has total width ``9t``, compared to
``[-4t, +4t]`` (width ``8t``) for the square lattice. The asymmetry
``|E_{\min}| = 6t \neq |E_{\max}| = 3t`` is a direct manifestation
of the absence of bipartite symmetry.

### Van Hove singularity

At certain energies within the band, the density of states
``g(E) = \sum_{\mathbf{k}} \delta(E - E(\mathbf{k}))`` diverges
logarithmically due to saddle points in the dispersion ``E(\mathbf{k})``.
These Van Hove singularities are responsible for electronic
instabilities (magnetism, superconductivity) in triangular-lattice
materials.

### K-point commensurability

The ``K``-point eigenvalue ``E = +3t`` appears in the finite-size spectrum
only when both ``L_x`` and ``L_y`` are divisible by 3, since the
``K``-point momenta ``(1/3, 2/3)`` and ``(2/3, 1/3)`` must be commensurate
with the discrete Brillouin zone grid. When this condition is met,
exactly 2 eigenvalues sit at ``+3t``.

### References

- G. H. Wannier, "Antiferromagnetism. The Triangular Ising Net",
  Phys. Rev. **79**, 357 (1950) -- triangular lattice frustration.
- T. Koretsune, M. Ogata, "Electronic structures of triangular lattice
  models", J. Phys. Soc. Jpn. **76**, 074706 (2007) -- NN
  tight-binding spectrum and Van Hove singularity.

### QAtlas API

```julia
# Sorted single-particle spectrum, 6×6 triangular PBC
λ = QAtlas.fetch(QAtlas.Triangular(), TightBindingSpectrum(); Lx=6, Ly=6, t=1.0)
# → 36 eigenvalues, ranging from -6.0 to +3.0
```

### Verification

| Test file | Method | What is checked |
|-----------|--------|-----------------|
| `test_triangular_tight_binding.jl` | Real-space ED via Lattice2D | ``\lambda_{\text{real}} = \lambda_{\text{Bloch}}`` for ``3 \times 3`` through ``6 \times 6`` |
| `test_triangular_tight_binding.jl` | Band edges | ``E_{\min} = -6t`` (unique), ``E_{\max} = +3t`` (when ``3 | L``) |
| `test_triangular_tight_binding.jl` | Frustration | Spectrum is NOT symmetric about zero |
| `test_triangular_tight_binding.jl` | ``K``-point degeneracy | ``\geq 2`` eigenvalues at ``+3t`` when both ``L_x, L_y \equiv 0 \pmod{3}`` |
| `test_triangular_tight_binding.jl` | Structural | ``\text{tr}\,H = 0`` |
| `test_triangular_tight_binding.jl` | ``t`` scaling | ``\lambda(t) = t \cdot \lambda(1)`` |
| `test_bloch_generic.jl` | Generic Bloch builder | Scalar formula ``=`` `bloch_tb_spectrum` |

---

## Connections

- **Lattice family**: Part of the [tight-binding model catalogue](index.md).
  The triangular lattice is frustrated (non-bipartite) like
  [Kagome](kagome.md), but has only one sublattice and no flat band.
- **Bipartite contrast**: The square lattice (also 1 sublattice) has a
  symmetric band ``[-4t, +4t]`` because it is bipartite. The triangular
  band ``[-6t, +3t]`` breaks this symmetry.
- **Magnetic frustration**: The non-bipartite structure underlies the
  classical Ising antiferromagnet frustration problem on the triangular
  lattice (Wannier 1950).
- **Methods**: Computed via the
  [Bloch Hamiltonian](../../../methods/bloch-hamiltonian/index.md) method.
