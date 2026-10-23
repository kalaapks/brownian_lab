import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import TextArea, HPacker, AnchoredOffsetbox

sigma = 1

def cholesky(n, seed=None, show=True):
    if np.ndim(n) > 0:                        # a list / array of n's
        return {int(k): cholesky(int(k), seed, show) for k in n}

    n = int(n)
    rng = np.random.default_rng(seed)

    t = np.arange(1, n + 1) / n               #1/n, ..., 1 (no t=0 like donkser and levy had)
    Sigma = np.minimum.outer(t, t)            #as in outer produt created by all ts from previous line and then np.minimum applied
    L = np.linalg.cholesky(Sigma)             # Sigma = L @ L.T
    Z = rng.normal(0, 1, size=n)              # a n-D Gaussian vector Z

    B = L @ Z                                 # now correlated like BM
    B = np.concatenate(([0.0], B))            # B(0) = 0
    t = np.concatenate(([0.0], t))
    B = sigma * B

    if show:
        print(f"\nn = {n}: (Time, Path) coordinates\n")
        for time_pt, path_val in zip(t, B):
            print(f"({time_pt:.4f}, {path_val:.4f})")
    return t, B


def plot_cholesky(ns, seed=None, ax=None, linewidth=1.2, alpha=0.75):
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6), dpi=120)
    else:
        fig = ax.figure

    paths = cholesky(ns, seed=seed, show=False)
    for n, (t, B) in paths.items():
        ax.plot(t, B, label=f"$n={n}$", linewidth=linewidth, alpha=alpha)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.minorticks_on()
    ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.grid(True, which="minor", linestyle=":", linewidth=0.5, alpha=0.3)

    ax.set_title("Cholesky Simulation: Exact Brownian Motion on a Grid",
                 fontsize=13, pad=12)
    ax.set_xlabel("$t$", fontsize=11)
    ax.set_ylabel(r"$B(t)$", fontsize=11)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    return fig, ax


def qv_cholesky(ns, seed=None, limit=50):
    results = {}
    for n in np.atleast_1d(ns):
        n = int(n)
        t, B = cholesky(n, seed=seed, show=False)
        qv = np.sum(np.diff(B) ** 2)
        results[n] = qv
        print(f"n={n}: The Quadratic variance is {qv:.6f}.")
        if n <= limit:
            print(f"      Plot this using plotqv_cholesky({n}).")
        else:
            print("      n is too large to plot cleanly.")
    return results


def plotqv_cholesky(n, seed=None, ax=None, show_labels=None):
    n = int(n)
    t, B = cholesky(n, seed=seed, show=False)
    e = np.diff(B)                       # e_1, ..., e_n
    qv = np.sum(e ** 2)

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6), dpi=120)
    else:
        fig = ax.figure

    ax.plot(t, B, color="steelblue", linewidth=1.2, alpha=0.8,
            marker="o", markersize=3, zorder=2)

    mid = (B[:-1] + B[1:]) / 2
    ax.errorbar(t[1:], mid, yerr=np.abs(e) / 2, fmt="none",
                ecolor="red", elinewidth=1.2, capsize=3, zorder=3)

    if show_labels is None:
        show_labels = n <= 30
    if show_labels:
        for i in range(1, n + 1):
            ax.annotate(f"$e_{{{i}}}$", (t[i], mid[i - 1]), xytext=(3, 0),
                        textcoords="offset points", fontsize=6,
                        color="red", va="center")

    props = dict(fontsize=13)
    box = HPacker(
        children=[
            TextArea("Quadratic Variation of the Cholesky Path (",
                     textprops=dict(color="black", **props)),
            TextArea(rf"$\sum_{{i=1}}^{{{n}}} e_i^2 = {qv:.4f}$",
                     textprops=dict(color="red", **props)),
            TextArea(")", textprops=dict(color="black", **props)),
        ],
        align="baseline", pad=0, sep=0)
    ax.add_artist(AnchoredOffsetbox(
        loc="lower center", child=box, frameon=False, borderpad=0,
        bbox_to_anchor=(0.5, 1.02), bbox_transform=ax.transAxes))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.minorticks_on()
    ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.grid(True, which="minor", linestyle=":", linewidth=0.5, alpha=0.3)

    ax.set_xlabel("$t$", fontsize=11)
    ax.set_ylabel(r"$B(t)$", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig, ax