import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import TextArea, HPacker, AnchoredOffsetbox

sigma = 1


def levy(level, seed=None, show=True):
    if np.ndim(level) > 0:                    # a list / array of levels
        if seed is None:
            seed = int(np.random.default_rng().integers(0, 2**32))
        return {int(k): levy(int(k), seed, show) for k in level}

    level = int(level)
    N = 2 ** level
    B = np.zeros(N + 1)                       # boxes for times 0, 1/N, ..., 1
    rng = np.random.default_rng(seed)         
    B[0] = 0.0
    B[N] = rng.normal(0, 1)                   

    #r goes to D_0,D_1,... one by one and for each r, k goes through each point and tries to assign values to the midpoints of the points in D_r, hence creating the fucntion on D_r+1
    for r in range(level):
        step = N // 2 ** r                    # interval length, in boxes
        eps = rng.normal(0, 1, size=2 ** r)   # one fresh epsilon per interval
        for k in range(2 ** r):
            left = k * step
            right = (k + 1) * step
            mid = left + step // 2
            B[mid] = (B[left] + B[right]) / 2 + eps[k] / 2 ** (r / 2 + 1)

    t = np.arange(N + 1) / N
    B = sigma * B

    if show:
        print(f"\nlevel = {level} ({N} steps): (Time, Path) coordinates\n")
        for time_pt, path_val in zip(t, B):
            print(f"({time_pt:.4f}, {path_val:.4f})")
    return t, B


def plot_levy(levels, seed=None, ax=None, linewidth=1.2, alpha=0.75):
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6), dpi=120)
    else:
        fig = ax.figure

    paths = levy(levels, seed=seed, show=False)
    for r, (t, B) in paths.items():
        ax.plot(t, B, label=rf"level ${r}$ ($N={2 ** r}$)",
                linewidth=linewidth, alpha=alpha)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.minorticks_on()
    ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.grid(True, which="minor", linestyle=":", linewidth=0.5, alpha=0.3)

    ax.set_title("Lévy Construction: Midpoint Refinement to Brownian Motion",
                 fontsize=13, pad=12)
    ax.set_xlabel("$t$", fontsize=11)
    ax.set_ylabel(r"$B(t)$", fontsize=11)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    return fig, ax


def qv_levy(levels, seed=None, limit=6):
    levels = np.atleast_1d(levels)
    if seed is None:
        seed = int(np.random.default_rng().integers(0, 2**32))

    results = {}
    for r in levels:
        r = int(r)
        t, B = levy(r, seed=seed, show=False)
        qv = np.sum(np.diff(B) ** 2)
        results[r] = qv
        print(f"level={r} (N={2 ** r}): The Quadratic variance is {qv:.6f}.")
        if r <= limit:
            print(f"      Plot this using plotqv_levy({r}).")
        else:
            print("      level is too large to plot cleanly.")
    return results


def plotqv_levy(level, seed=None, ax=None, show_labels=None):
    level = int(level)
    N = 2 ** level
    t, B = levy(level, seed=seed, show=False)
    e = np.diff(B)                       # e_1, ..., e_N
    qv = np.sum(e ** 2)

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6), dpi=120)
    else:
        fig = ax.figure

    #the path
    ax.plot(t, B, color="steelblue", linewidth=1.2, alpha=0.8,
            marker="o", markersize=3, zorder=2)

    # one red vertical segment (with bars) at each t_i, spanning B[i-1] to B[i]
    mid = (B[:-1] + B[1:]) / 2
    ax.errorbar(t[1:], mid, yerr=np.abs(e) / 2, fmt="none",
                ecolor="red", elinewidth=1.2, capsize=3, zorder=3)

    # tiny e_i labels, switched off automatically when they clutter
    if show_labels is None:
        show_labels = N <= 30
    if show_labels:
        for i in range(1, N + 1):
            ax.annotate(f"$e_{{{i}}}$", (t[i], mid[i - 1]), xytext=(3, 0),
                        textcoords="offset points", fontsize=6,
                        color="red", va="center")

    props = dict(fontsize=13)
    box = HPacker(
        children=[
            TextArea("Quadratic Variation of the Lévy Path (",
                     textprops=dict(color="black", **props)),
            TextArea(rf"$\sum_{{i=1}}^{{{N}}} e_i^2 = {qv:.4f}$",
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