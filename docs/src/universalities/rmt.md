# Random Matrix Theory (RMT) and Poisson Universality

## Overview

A spectrum of energy levels can be classified by the statistics of its
*nearest-neighbour level spacings* and *largest-eigenvalue extremes*.
For a wide class of complex quantum systems --- chaotic many-body
Hamiltonians, heavy nuclei, disordered conductors --- the local
statistics are universal and depend only on the underlying global
symmetry through the Dyson index $\beta \in \{1, 2, 4\}$
(Wigner-Dyson three-fold way):

| Symmetry class               | $\beta$ | Ensemble |
|------------------------------|---------|----------|
| Time-reversal invariant      | $1$     | GOE      |
| No time-reversal             | $2$     | GUE      |
| TRI + half-integer spin only | $4$     | GSE      |

An *integrable* (or many-body localised) system, in contrast, has a
spectrum that is asymptotically uncorrelated and produces **Poisson**
level statistics.

QAtlas exposes both classes through the parametric `Universality{C}`
API (`C` = `:RMT` or `:Poisson`).

---

## Wigner Surmise --- Closed Form

The Wigner surmise is the exact $N = 2$ Gaussian-ensemble spacing
distribution; it is also a celebrated, accurate approximation to the
$N \to \infty$ bulk spacing distribution.

$$P_1(s) = \frac{\pi s}{2} \exp\!\left(-\frac{\pi s^2}{4}\right) \qquad \text{(GOE)}$$

$$P_2(s) = \frac{32 s^2}{\pi^2} \exp\!\left(-\frac{4 s^2}{\pi}\right) \qquad \text{(GUE)}$$

$$P_4(s) = \frac{2^{18} s^4}{3^6 \pi^3} \exp\!\left(-\frac{64 s^2}{9\pi}\right) \qquad \text{(GSE)}$$

All three are normalised so that $\int_0^\infty P_\beta(s)\,ds = 1$
and $\int_0^\infty s\,P_\beta(s)\,ds = 1$ (mean spacing $= 1$).
Small-$s$ behaviour $P_\beta(s) \sim s^\beta$ is the **level
repulsion** characteristic of each ensemble.

The Poisson counterpart is

$$P_{\text{Poisson}}(s) = e^{-s},$$

with no level repulsion and the same mean spacing.

---

## Tracy-Widom $F_\beta$ --- Largest-Eigenvalue Distribution

In the limit $N \to \infty$, the largest eigenvalue of a Gaussian
$\beta$-ensemble matrix obeys

$$\mathbb{P}\!\left[\lambda_{\max} \le \lambda_c + \sigma N^{-2/3} x\right] \to F_\beta(x),$$

where $F_\beta$ is the **Tracy-Widom** distribution. Closed forms
require the Painleve II transcendent $q(x)$ via

$$F_2(x) = \exp\!\left[-\int_x^\infty (s - x)\, q(s)^2\,ds\right],$$

with $F_1$ and $F_4$ given by analogous Pfaffian formulas
(Tracy-Widom 1996).

### Phase 1 implementation

QAtlas Phase 1 evaluates $F_\beta$ from an embedded high-precision
table compiled from
[Bornemann (Math. Comp. **79**, 871, 2010), Table 1](https://www.ams.org/journals/mcom/2010-79-270/S0025-5718-09-02280-7/),
covering $x \in [-4, 4]$ for all three $\beta$. Inside the table
support the interpolant is piecewise-linear and monotone
non-decreasing; outside the support the function returns the
Tracy-Widom 1994/1996 tail asymptotics

$$F_\beta(x) \sim \tau_\beta \exp\!\left[-\tfrac{\beta}{24} |x|^3\right] \qquad (x \to -\infty),$$

$$1 - F_\beta(x) \sim \exp\!\left[-\tfrac{2\beta}{3} x^{3/2}\right] \qquad (x \to +\infty),$$

continuously matched at the table boundary. Reference checkpoints
pinned by the standalone test:

| $x = 0$              | $F_\beta(0)$         |
|----------------------|----------------------|
| $\beta = 1$ (GOE)    | $\approx 0.8319$     |
| $\beta = 2$ (GUE)    | $\approx 0.9694$     |
| $\beta = 4$ (GSE)    | $\approx 0.99966$    |

A direct Painleve-II ODE integrator (DifferentialEquations.jl-based)
is deferred to **Phase 2** of issue #151.

---

## Mean Ratio $\langle r \rangle$

For a level sequence $\{E_n\}$ with spacings $s_n = E_{n+1} - E_n$,
the consecutive-spacing ratio

$$r_n = \frac{\min(s_n, s_{n+1})}{\max(s_n, s_{n+1})}$$

is dimensionless (no spectral unfolding needed). Its mean takes
universal values:

| Ensemble                   | $\langle r \rangle$              |
|----------------------------|----------------------------------|
| Poisson (integrable / MBL) | $2 \log 2 - 1 \approx 0.3863$    |
| GOE ($\beta = 1$)          | $\approx 0.5307$                 |
| GUE ($\beta = 2$)          | $\approx 0.5996$                 |
| GSE ($\beta = 4$)          | $\approx 0.6744$                 |

QAtlas returns these literature values directly from
[Atas-Bogomolny-Giraud-Roux, Phys. Rev. Lett. **110**, 084101 (2013)](https://doi.org/10.1103/PhysRevLett.110.084101).

---

## QAtlas API

```julia
using QAtlas

# Wigner-Dyson surmise
P_GUE_at_1 = QAtlas.fetch(Universality(:RMT), WignerSurmise(); beta=2, s=1.0)

# Tracy-Widom CDF
F_GUE_at_0 = QAtlas.fetch(Universality(:RMT), TracyWidom(); beta=2, x=0.0)
# ~ 0.9694 (Bornemann 2010 Table 1)

# Mean ratio
r_GUE = QAtlas.fetch(Universality(:RMT), MeanRatio(); beta=2)   # 0.5996

# Poisson (integrable baseline)
P_int = QAtlas.fetch(Universality(:Poisson), WignerSurmise(); s=1.0)  # 1/e
r_int = QAtlas.fetch(Universality(:Poisson), MeanRatio())             # 2 log 2 - 1
```

(In the actual code use the Greek letter `beta` keyword as `β`.)

---

## References

- M. L. Mehta, *Random Matrices*, 3rd ed., Elsevier (2004).
- E. P. Wigner, *Conference on Neutron Physics by Time-of-Flight*,
  Oak Ridge Natl. Lab. Rep. ORNL-2309, 59 (1957) --- surmise.
- F. J. Dyson, J. Math. Phys. **3**, 140 (1962) --- three-fold way.
- C. A. Tracy, H. Widom, *Level-spacing distributions and the Airy
  kernel*, Commun. Math. Phys. **159**, 151 (1994) --- $F_2$.
- C. A. Tracy, H. Widom, *On orthogonal and symplectic matrix
  ensembles*, Commun. Math. Phys. **177**, 727 (1996) --- $F_1, F_4$.
- F. Bornemann, *On the numerical evaluation of Fredholm determinants*,
  Math. Comp. **79**, 871 (2010) --- high-precision $F_\beta$ table.
- Y. Y. Atas, E. Bogomolny, O. Giraud, G. Roux, *Distribution of the
  ratio of consecutive level spacings in random matrix ensembles*,
  Phys. Rev. Lett. **110**, 084101 (2013) --- $\langle r \rangle$.
- V. Oganesyan, D. A. Huse, Phys. Rev. B **75**, 155111 (2007) ---
  ratio statistic in many-body localisation.
