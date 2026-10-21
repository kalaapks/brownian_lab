import numpy as np
np.set_printoptions(precision=4)

sigma = 1

def donsker_walk(n, seed=None):
    rng = np.random.default_rng(seed)
    xi = rng.normal(0, sigma, n)
    S = np.concatenate([[0.0], np.cumsum(xi)])
    t = np.arange(n + 1) / n
    return t, S / np.sqrt(n)
    