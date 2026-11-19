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
        fig.savefig(save, dpi=dpi)
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
