# Model description cards (the "what is this model" one-liner + Hamiltonian).
#
# Included after all model types are defined.  Each `@about` row gives a model
# its Wikipedia-first-sentence `summary` and its `hamiltonian` in LaTeX, shown
# at the top of the model page.  Authored incrementally — a model without a
# card falls back to its struct docstring in the docs generator, so this set
# can grow over time without leaving any model blank.
#
# Use `raw"…"` for both fields: it preserves LaTeX backslashes and lets the
# `summary` carry inline `$…$` math without Julia string interpolation.

@about TFIM summary = raw"The 1D transverse-field Ising model — the canonical exactly-solvable quantum phase transition, Jordan-Wigner dual to free fermions and critical at $h = J$." hamiltonian = raw"H = -J\sum_i \sigma^z_i \sigma^z_{i+1} - h\sum_i \sigma^x_i"

@about Heisenberg1D summary = raw"The spin-$\tfrac12$ antiferromagnetic Heisenberg chain — the Bethe-ansatz-integrable model of quantum magnetism, gapless with an $SU(2)_1$ WZW critical point ($c = 1$)." hamiltonian = raw"H = J\sum_{\langle i,j\rangle} \mathbf{S}_i \cdot \mathbf{S}_j, \qquad J > 0"

@about XXZ1D summary = raw"The spin-$\tfrac12$ XXZ chain — the uniaxially anisotropic Heisenberg model, Bethe-ansatz integrable, with a critical Luttinger-liquid line for $-1 < \Delta \le 1$." hamiltonian = raw"H = \sum_i \left( S^x_i S^x_{i+1} + S^y_i S^y_{i+1} + \Delta\, S^z_i S^z_{i+1} \right)"

@about Hubbard1D summary = raw"The 1D Hubbard model — the minimal lattice model of interacting electrons, Lieb-Wu Bethe-ansatz solvable; a Mott insulator with spin-charge separation at half filling." hamiltonian = raw"H = -t\sum_{i,\sigma} \left( c^\dagger_{i\sigma} c_{i+1\sigma} + \text{h.c.} \right) + U\sum_i n_{i\uparrow} n_{i\downarrow}"

@about IsingSquare summary = raw"The 2D classical Ising model on the square lattice — Onsager's exactly-solved ferromagnet with a finite-temperature order-disorder transition." hamiltonian = raw"H = -J\sum_{\langle i,j\rangle} \sigma_i \sigma_j, \qquad \sigma_i = \pm 1"

@about DimerLattice summary = raw"The close-packed dimer model — perfect matchings (dominoes) tiling the square lattice, solved exactly by the Kasteleyn-Temperley-Fisher Pfaffian, with residual entropy $G/\pi$ per site." hamiltonian = raw"Z = \#\{\text{perfect matchings of the } L_x \times L_y \text{ grid}\}"

@about IsingTriangular summary = raw"The 2D classical Ising model on the triangular lattice — exactly solved and in the 2D Ising universality class; its antiferromagnet is geometrically frustrated." hamiltonian = raw"H = -J\sum_{\langle i,j\rangle} \sigma_i \sigma_j, \qquad \sigma_i = \pm 1"

@about Kitaev1D summary = raw"The Kitaev chain — a 1D spinless $p$-wave superconductor whose topological phase binds unpaired Majorana zero modes at its ends." hamiltonian = raw"H = \sum_i \left( -t\, c^\dagger_i c_{i+1} + \Delta\, c_i c_{i+1} + \text{h.c.} \right) - \mu\sum_i c^\dagger_i c_i"

@about SSH summary = raw"The Su-Schrieffer-Heeger chain — a 1D dimerised tight-binding model (chiral class BDI) whose topological phase ($|w|>|v|$) carries winding number $W=1$ and protected end zero modes." hamiltonian = raw"H = \sum_i \left( v\, c^\dagger_{i,A} c_{i,B} + w\, c^\dagger_{i,B} c_{i+1,A} + \text{h.c.} \right)"

@about KitaevHoneycomb summary = raw"The Kitaev honeycomb model — an exactly-solvable bond-dependent spin model realizing a $\mathbb{Z}_2$ quantum spin liquid with emergent Majorana fermions." hamiltonian = raw"H = -\sum_{\langle i,j\rangle_\gamma} J_\gamma\, \sigma^\gamma_i \sigma^\gamma_j, \qquad \gamma \in \{x,y,z\}"

@about ToricCode summary = raw"The toric code — Kitaev's exactly-solvable stabilizer model, the paradigmatic $\mathbb{Z}_2$ topological order with deconfined anyonic excitations." hamiltonian = raw"H = -\sum_v A_v - \sum_p B_p, \qquad A_v = \prod_{i\in v}\sigma^x_i,\ \ B_p = \prod_{i\in p}\sigma^z_i"

@about MajumdarGhosh summary = raw"The Majumdar-Ghosh chain — the $J_2 = J_1/2$ frustrated Heisenberg chain with an exact, doubly-degenerate nearest-neighbour dimerized ground state." hamiltonian = raw"H = J_1\sum_i \mathbf{S}_i \cdot \mathbf{S}_{i+1} + J_2\sum_i \mathbf{S}_i \cdot \mathbf{S}_{i+2}, \qquad J_2 = J_1/2"
