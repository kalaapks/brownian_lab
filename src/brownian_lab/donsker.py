import numpy as np
import matplotlib.pyplot as plt

sigma = 1


def donsker_walk(n, seed=None, show=True):
    if np.ndim(n) > 0:                       # a list / array of n's
        return {int(k): donsker_walk(int(k), seed, show) for k in n}

    rng = np.random.default_rng(seed)
    xi = rng.normal(0, sigma, n)
    S = np.concatenate([[0.0], np.cumsum(xi)])
    t = np.arange(n + 1) / n
    W = S / (sigma*np.sqrt(n))
    if show:
        print(f"\nn = {n}: (Time, Walk) coordinates\n")
        for time_pt, walk_val in zip(t, W):
            print(f"({time_pt:.2f}, {walk_val:.4f})")
    return t, W


def plot_donsker(ns, seed=None, ax=None, linewidth=1.2, alpha=0.75):
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6), dpi=120)
    else:
        fig = ax.figure

    paths = donsker_walk(ns, seed=seed, show=False)
    for n, (t, W) in paths.items():
        ax.plot(t, W, label=f"$n={n}$", linewidth=linewidth, alpha=alpha)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.minorticks_on()
    ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.grid(True, which="minor", linestyle=":", linewidth=0.5, alpha=0.3)

    ax.set_title("Donsker Invariance Principle: Scaled Random Walks to Brownian Motion",
                 fontsize=13, pad=12)
    ax.set_xlabel("$t$", fontsize=11)
    ax.set_ylabel(r"$W_n(t) = \frac{S_{\lfloor nt \rfloor}}{\sqrt{n}}$", fontsize=11)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    return fig, ax
    
def qv_donsker(ns, seed=None, limit=50):

    results = {}
    for n in np.atleast_1d(ns):
        n = int(n)
        t, W = donsker_walk(n, seed=seed, show=False)
        qv = np.sum(np.diff(W) ** 2)
        results[n] = qv
        print(f"n={n}: The Quadratic variance is {qv:.6f}.")
        if n <= limit:
            print("      Plot this using plotqv_donsker(n, seed).")
        else:
            print("      n is too large to plot cleanly.")
    return results

def plotqv_donsker(n, seed=None, ax=None, show_labels=None):
    
    t, W = donsker_walk(int(n), seed=seed, show=False)
    e = np.diff(W)                       # e_1, ..., e_n
    qv = np.sum(e ** 2)

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6), dpi=120)
    else:
        fig = ax.figure

    # the path itself
    ax.plot(t, W, color="steelblue", linewidth=1.2, alpha=0.8,
            marker="o", markersize=3, zorder=2)

    # one red vertical segment (with caps) at each t_i, spanning W[i-1] to W[i]
    mid = (W[:-1] + W[1:]) / 2
    ax.errorbar(t[1:], mid, yerr=np.abs(e) / 2, fmt="none",
                ecolor="red", elinewidth=1.2, capsize=3, zorder=3)

    # tiny e_i labels, switched off automatically when they would clutter
    if show_labels is None:
        show_labels = n <= 30
    if show_labels:
        for i in range(1, int(n) + 1):
            ax.annotate(f"$e_{{{i}}}$", (t[i], mid[i - 1]), xytext=(3, 0),
                        textcoords="offset points", fontsize=6,
                        color="red", va="center")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.minorticks_on()
    ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.grid(True, which="minor", linestyle=":", linewidth=0.5, alpha=0.3)

    ax.set_title(f"Quadratic Variation of the Donsker Walk ($n={int(n)}$)",
                 fontsize=13, pad=28)
    ax.set_title(rf"$\sum_{{i=1}}^{{{int(n)}}} e_i^2 = {qv:.4f}$",
                 fontsize=13, color="red", pad=10)
    ax.set_xlabel("$t$", fontsize=11)
    ax.set_ylabel(r"$W_n(t)$", fontsize=11)
    fig.tight_layout()
    return fig, ax