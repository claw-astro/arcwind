"""PV-diagram figures: one diagram, a grid of cuts, a high-contrast gallery."""
from __future__ import annotations

import numpy as np

from ..pv import contour_levels
from ._style import cut_label, draw_curve, finish, grid, label_grid, plt, tidy_grid_ticks

__all__ = ["plot_pv", "plot_pv_grid", "plot_gallery"]


def plot_pv(pv, ax=None, contours=True, show_mask=False, cmap="magma", gamma=0.5,
            start_sigma=5.0, step_sigma=5.0, model=None, trace=None, xlim=None, vlim=None,
            save=None):
    """Show one :class:`~diskwind.pv.PVDiagram`.

    Parameters
    ----------
    pv : PVDiagram
    ax : matplotlib Axes, optional
    contours : bool
        Contours from ``start_sigma`` in steps of ``step_sigma`` x rms.
    show_mask : bool
        Outline the pixels excluded from the trace (cyan).
    cmap, gamma : str, float
        Color map and power-law stretch.
    model : ModelCurve, FitResult or list, optional
        Curve(s) to overlay (a FitResult is evaluated at ``pv.y0``).
    trace : Trace, optional
        Trace points to overlay.
    xlim, vlim : tuple, optional
    save : str, optional
        File to write.

    Returns
    -------
    matplotlib Axes

    Examples
    --------
    >>> ax = cube.pv(1.0).plot(show_mask=True, model=result)     # doctest: +SKIP
    """
    from matplotlib.colors import PowerNorm
    p = plt()
    if ax is None:
        fig, ax = p.subplots(figsize=(5, 4))
    d = pv.data
    ax.imshow(np.clip(d, 0, None), origin="lower", aspect="auto", extent=pv.extent, cmap=cmap,
              interpolation="nearest",
              norm=PowerNorm(gamma, vmin=0, vmax=max(np.nanmax(d), 5 * pv.rms)))
    if contours:
        lev = contour_levels(pv.rms, np.nanmax(d), start_sigma, step_sigma)
        ax.contour(pv.x, pv.v, d, levels=lev, colors="w", linewidths=0.35, alpha=0.8)
    if show_mask and pv.excluded is not None:
        ax.contour(pv.x, pv.v, pv.excluded.astype(float), levels=[0.5], colors="c",
                   linewidths=0.5)
    if trace is not None and trace.n:
        ax.errorbar(trace.x, trace.v, yerr=trace.err, fmt="+", ms=3, color="k", lw=0.5)
    for m in ([] if model is None else (model if isinstance(model, (list, tuple)) else [model])):
        c = m.curve(pv.y0) if hasattr(m, "curve") else m
        draw_curve(ax, c, color="lime", lw=1.0)
    ax.text(0.03, 0.95, cut_label(pv.y0, pv.zc_au), transform=ax.transAxes, color="w",
            va="top", fontsize=7.5)
    if xlim:
        ax.set_xlim(*xlim)
    if vlim:
        ax.set_ylim(*vlim)
    ax.set_xlabel("X (arcsec)")
    ax.set_ylabel(r"$V_{\rm LSR}$ (km s$^{-1}$)")
    finish(ax.figure, save)
    return ax


def plot_pv_grid(cube, y_offsets=None, result=None, show_mask=False, ncols=None, save=None):
    """Transverse PV diagrams of many cuts (Louvet+2018 Fig. 7 style).

    Parameters
    ----------
    cube : WindCube
    y_offsets : array-like, optional
        Default: ``plot.pv_grid_offsets_arcsec`` or the fit cuts.
    result : FitResult, optional
        Model curve to overlay on the cuts with y0 > 0.
    show_mask : bool
    ncols : int, optional

    Returns
    -------
    Figure
    """
    cfg = cube.config
    pc = cfg.plot
    from ..config import _offsets
    ys = (np.atleast_1d(y_offsets) if y_offsets is not None else
          (_offsets(pc.pv_grid_offsets_arcsec) if pc.pv_grid_offsets_arcsec is not None
           else cfg.y_offsets))
    fig, axs = grid(len(ys), ncols or pc.ncols)
    for ax, y0 in zip(axs.flat, ys):
        pv = cube.pv(y0, with_mask=show_mask)
        model = [result] if (result is not None and y0 > 0) else None
        plot_pv(pv, ax=ax, show_mask=show_mask, cmap=pc.cmap, start_sigma=pc.contour_start_sigma,
                step_sigma=pc.contour_step_sigma, model=model)
        ax.plot(0, cube.vsys, "+", color="cyan", ms=8, mew=1.2)
        ax.set_xlabel("")
        ax.set_ylabel("")
    a0 = axs.flat[0]
    a0.set_xlim(*(pc.xlim or (cube.x[0], cube.x[-1])))
    if pc.vlim:
        a0.set_ylim(*pc.vlim)
    tidy_grid_ticks(axs)
    fig.supxlabel("X: distance to the wind axis (arcsec)")
    fig.supylabel(r"$V_{\rm LSR}$ (km s$^{-1}$)")
    fig.suptitle(f"{cfg.source.name}: transverse PV diagrams (contours from "
                 f"{pc.contour_start_sigma:g} sigma; + = (0, Vsys))", fontsize=9)
    fig.subplots_adjust(wspace=0.04, hspace=0.04)
    return finish(fig, save, cfg.output.dpi)


def plot_gallery(cube, y_offsets=None, result=None, n_boot_curves=0, gamma=0.5,
                 contour_sigma=3.0, cmap="Greys", common_scale=False, show_mask=False,
                 ncols=4, panel=(3.0, 2.5), seed=42, save=None):
    """High-contrast gallery of transverse PV diagrams along the wind.

    Parameters
    ----------
    cube : WindCube
    y_offsets : array-like, optional
        Default: the fit cuts.
    result : FitResult, optional
        Model curve (red) and ring center (x); ``n_boot_curves`` bootstrap curves (gold).
    gamma : float
        Power-law stretch.
    contour_sigma : float
        First contour and half the step, in units of the rms.
    common_scale : bool
        One intensity scale for all panels (default: each panel to its peak).

    Returns
    -------
    Figure
    """
    from matplotlib.colors import PowerNorm
    cfg = cube.config
    ys = np.atleast_1d(y_offsets) if y_offsets is not None else cfg.y_offsets
    pvs = {y: cube.pv(y, with_mask=show_mask) for y in ys}
    gmax = max(np.nanmax(p.data) for p in pvs.values())
    fig, axs = grid(len(ys), ncols, panel)
    xl = cfg.plot.xlim or (cube.x[0], cube.x[-1])
    rng = np.random.default_rng(seed)
    draws = result.draws() if (result is not None and n_boot_curves) else None
    for ax, y0 in zip(axs.flat, ys):
        pv = pvs[y0]
        vmax = 0.6 * (gmax if common_scale else np.nanmax(pv.data))
        ax.imshow(np.clip(pv.data, 0, None), origin="lower", aspect="auto", extent=pv.extent,
                  cmap=cmap, norm=PowerNorm(gamma, vmin=0, vmax=vmax), interpolation="nearest",
                  rasterized=True)
        lev = contour_levels(pv.rms, np.nanmax(pv.data), contour_sigma, 2 * contour_sigma)
        ax.contour(pv.x, pv.v, pv.data, levels=lev, colors="0.45", linewidths=0.35)
        if show_mask and pv.excluded is not None:
            ax.contour(pv.x, pv.v, pv.excluded.astype(float), levels=[0.5], colors="c",
                       linewidths=0.5)
        if result is not None and y0 > 0:
            if draws is not None:
                for k in rng.choice(len(draws), size=n_boot_curves,
                                    replace=len(draws) < n_boot_curves):
                    draw_curve(ax, result.curve(y0, draws[k]), color="gold", lw=0.4, alpha=0.4)
            draw_curve(ax, result.curve(y0), color="tab:red", lw=1.2)
            rg = result.ring_at(y0)
            ax.plot(rg["xc"] / cube.dist, rg["vcent"], "x", color="tab:red", ms=6, mew=1.5)
        ax.axvline(0, color="0.6", lw=0.5, ls=":")
        ax.axhline(cube.vsys, color="0.6", lw=0.5, ls=":")
        ax.text(0.04, 0.95, f"y = {y0:+.2f}″\nz ≈ {cube.zc_au(y0):.0f} au",
                transform=ax.transAxes, va="top", fontsize=8,
                bbox=dict(fc="w", ec="none", alpha=0.75, pad=1.5))
        ax.set_xlim(*xl)
        if cfg.plot.vlim:
            ax.set_ylim(*cfg.plot.vlim)
        ax.tick_params(top=False, right=False, direction="out", labelsize=8)
    tidy_grid_ticks(axs)
    fig.subplots_adjust(wspace=0.05, hspace=0.07)
    label_grid(fig, axs, "distance to the wind axis, X (arcsec)", r"velocity (km s$^{-1}$)",
               f"{cfg.source.name}: transverse PV diagrams along the wind")
    return finish(fig, save, cfg.output.dpi)
