"""Shared figure helpers."""
from __future__ import annotations

import numpy as np

#: Plot style of every diskwind figure: matplotlib's 'classic' style with these
#: settings on top.
STYLE = "classic"
RC = {"font.size": 14, "axes.linewidth": 2.0, "font.family": "sans-serif",
      "font.weight": "500", "mathtext.bf": "serif:normal",
      "savefig.bbox": "tight",
      # readability on top of 'classic': no offset labels, scientific notation for
      # very small/large numbers, a little room above and below curves
      "axes.formatter.useoffset": False, "axes.formatter.limits": (-4, 5),
      "axes.ymargin": 0.05}


def set_style():
    """Apply the diskwind plot style (see ``STYLE`` and ``RC``)."""
    import matplotlib.pyplot as _plt
    _plt.style.use(STYLE)
    _plt.rcParams.update(RC)


def plt():
    """matplotlib.pyplot with the diskwind style applied."""
    import matplotlib.pyplot as _plt
    set_style()
    return _plt


def thin_ticks(fig, nbins=5):
    """At most ~``nbins`` labeled ticks on the default (linear) axes, so the large
    tick labels of the style do not run into each other."""
    from matplotlib.ticker import AutoLocator, MaxNLocator
    for ax in fig.axes:
        for axis, scale in ((ax.xaxis, ax.get_xscale()), (ax.yaxis, ax.get_yscale())):
            if scale == "linear" and type(axis.get_major_locator()) is AutoLocator:
                axis.set_major_locator(MaxNLocator(nbins, steps=[1, 2, 2.5, 5, 10]))


def finish(fig, save=None, dpi=200):
    """Save the figure if ``save`` is a path; return it."""
    thin_ticks(fig)
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
        fig.text(0.5 * (x0 + x1), y1 + 0.12 / h, title, ha="center", va="bottom", fontsize="large")
