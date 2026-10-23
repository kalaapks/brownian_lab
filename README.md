# brownian-lab

Three constructions of Brownian motion on [0, 1] in Python, with plotting and
quadratic-variation visualisation. Companion code to my master's thesis on
stochastic calculus and the Black-Scholes PDE (IISER Kolkata) where the rigorous proofs of the results we are coding here can be found.

## The three generators

| Function | Construction | Grid |
|---|---|---|
| `donsker_walk(n, seed)` | Scaled Gaussian random walk (Donsker's theorem) | n steps |
| `levy(level, seed)` | Lévy midpoint-displacement construction | 2^level steps |
| `cholesky(n, seed)` | Exact simulation via Cholesky factor of min(s, t) | n steps |

Each generator returns `(t, W)`, the time grid and the path values. Passing a
list returns a dictionary with one path per entry.

Each generator also has `plot_*`, `qv_*` and `plotqv_*` functions, e.g.
`plot_levy`, `qv_levy`, `plotqv_levy`.

## Installation

```
pip install git+https://github.com/YOUR_USERNAME/brownian-lab.git
```

## Usage

```python
from brownian_lab import donsker_walk, plot_donsker, qv_donsker, plotqv_donsker

t, W = donsker_walk(100, 43)          # prints coordinates, returns (t, W)
plot_donsker([10, 100, 1000], 43)     # several paths on one plot
qv_donsker([10, 20, 50, 60], 43)      # quadratic variation for each n
plotqv_donsker(20, 43)                # path with increments e_i drawn in red
```

The same pattern works for `levy` and `cholesky`:

```python
from brownian_lab import levy, plot_levy, plotqv_levy

plot_levy([2, 4, 6, 10], 43)          # same seed: each level refines the last
plotqv_levy(4, 43)
```

## Plots

![Donsker walks](images/donsker.png)
![Lévy refinement](images/levy.png)
![Quadratic variation](images/plotqv.png)

## Background

- **Donsker:** the rescaled walk W_n(t) = S_{⌊nt⌋}/√n converges in law to
  Brownian motion.
- **Lévy:** the path is built coarse to fine, setting each midpoint to the
  average of its endpoints plus an independent Gaussian correction.
- **Cholesky:** B = LZ where LLᵀ = Σ, Σ_ij = min(t_i, t_j). It is exact, but costs
  O(n³) to factorise.
- **Quadratic variation:** Σ(ΔW)² tends to t as the grid is refined.

## Status

Generators and plotting are complete. Numerical validation (increment and
covariance tests, convergence rates) is in progress.

## Requirements

Python ≥ 3.9, NumPy, SciPy, Matplotlib.
