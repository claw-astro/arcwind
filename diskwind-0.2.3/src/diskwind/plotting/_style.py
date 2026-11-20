"""Shared figure helpers."""
from __future__ import annotations

import numpy as np

RC = {"font.size": 9, "axes.linewidth": 0.8, "xtick.direction": "in",
      "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
      "savefig.bbox": "tight"}


def plt():
    """matplotlib.pyplot with the diskwind defaults applied."""
    import matplotlib.pyplot as _plt
    _plt.rcParams.update(RC)
    return _plt


def finish(fig, save=None, dpi=200):
    """Save the figure if ``save`` is a path; return it."""
    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight", pad_inches=0.05)
    return fig


def cut_label(y0, zc=None):
    return f"y = {y0:+.2f}\"" + (f"  (z ≈ {zc:.0f} au)" if zc is not None and y0 > 0 else "")


def draw_curve(ax, curve, **kw):
    X, v = curve.masked()
    return ax.plot(X, v, **kw)


def grid(n, ncols, panel=(3.0, 2.5), sharex=True, sharey=True):
    p = plt()
    nc = min(ncols, n)
    nr = int(np.ceil(n / nc))
    fig, axs = p.subplots(nr, nc, figsize=(panel[0] * nc, panel[1] * nr), sharex=sharex,
                          sharey=sharey, squeeze=False)
    for ax in axs.flat[n:]:
        ax.set_visible(False)
    return fig, axs


def tidy_grid_ticks(axs, nbins=5):
    """Shared-axis panel grids: drop the end ticks so labels of neighboring
    panels do not run into each other."""
    from matplotlib.ticker import MaxNLocator
    for ax in np.ravel(axs):
        ax.xaxis.set_major_locator(MaxNLocator(nbins, steps=[1, 2, 5, 10], prune="both"))


def label_grid(fig, axs, xlabel, ylabel, title=None, pad_in=0.45):
    """Shared x/y labels and a title for a panel grid, placed a fixed distance
    (in inches) outside the panels, so they never touch tick labels, also for
    single-row grids."""
    boxes = [a.get_position() for a in np.ravel(axs) if a.get_visible()]
    x0 = min(b.x0 for b in boxes); x1 = max(b.x1 for b in boxes)
    y0 = min(b.y0 for b in boxes); y1 = max(b.y1 for b in boxes)
    w, h = fig.get_size_inches()
    fig.text(0.5 * (x0 + x1), y0 - pad_in / h, xlabel, ha="center", va="top")
    fig.text(x0 - (pad_in + 0.12) / w, 0.5 * (y0 + y1), ylabel, ha="right", va="center",
             rotation=90)
    if title:
        fig.text(0.5 * (x0 + x1), y1 + 0.12 / h, title, ha="center", va="bottom", fontsize=11)
