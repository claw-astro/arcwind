"""Sky maps: the layout of the cuts on the moment-0 / continuum images."""
from __future__ import annotations

import warnings

import numpy as np

from ..io import CelestialWCS, parse_angle, read_fits
from ._style import finish, plt

__all__ = ["plot_slit_layout", "load_image"]


def load_image(path):
    """2-D image (mJy/beam [km/s]), its WCS and beam (arcsec, arcsec, deg)."""
    hdr, data, beams = read_fits(path)
    img = np.asarray(data, float)
    while img.ndim > 2:
        img = img[0]
    wcs = CelestialWCS(hdr, 1, 2)
    if beams is not None and len(beams):
        bmaj, bmin, bpa = np.median(np.array([list(r)[:3] for r in beams], float), axis=0)
    else:
        bmaj, bmin, bpa = hdr["BMAJ"] * 3600.0, hdr["BMIN"] * 3600.0, hdr["BPA"]
    fac = 1.0 if str(hdr.get("BUNIT", "Jy/beam")).lower().startswith("mjy") else 1e3
    return img * fac, wcs, (bmaj, bmin, bpa)


def _extent(wcs, shape, ra0, de0, off):
    if abs(wcs.cd[0, 1]) > 1e-3 * abs(wcs.cd[0, 0]):
        warnings.warn("image axes are rotated w.r.t. RA/Dec; extent is approximate")
    px, py = wcs.world_to_pixel0(ra0, de0)
    pc = np.array([px, py]) + wcs.cdinv @ (np.asarray(off, float) / 3600.0)
    ny, nx = shape
    return [wcs.cd[0, 0] * (-0.5 - pc[0]) * 3600, wcs.cd[0, 0] * (nx - 0.5 - pc[0]) * 3600,
            wcs.cd[1, 1] * (-0.5 - pc[1]) * 3600, wcs.cd[1, 1] * (ny - 0.5 - pc[1]) * 3600]


def plot_slit_layout(cube, y_offsets=None, mom0=None, continuum=None, half_size=3.5,
                     gamma=0.75, mom0_min=None, cmap="bone", cont_min=None, cont_gamma=0.5,
                     figsize=(5.6, 4.8), save=None):
    """The transverse cuts on the sky, with the wind-frame axes and the beams.

    Parameters
    ----------
    cube : WindCube
    y_offsets : array-like, optional
        Cuts to draw (default: ``plot.pv_grid_offsets_arcsec`` or the fit cuts).
    mom0 : str, optional
        Moment-0 FITS file; default: the integrated intensity of the cube.
    continuum : str, optional
        Continuum FITS file drawn on top (inferno).
    half_size : float
        Half-width of the field (arcsec).
    gamma, mom0_min, cmap : stretch, floor and colour map of the moment-0 image.
    cont_min, cont_gamma : floor and stretch of the continuum.

    Returns
    -------
    Figure

    Examples
    --------
    >>> plot_slit_layout(cube, mom0="Elias20_CO_M0.fits",
    ...                  continuum="Elias20_continuum.fits")      # doctest: +SKIP
    """
    from matplotlib.colors import PowerNorm
    from matplotlib.patches import Ellipse
    from ..config import _offsets
    p = plt()
    cfg = cube.config
    src = cfg.source
    try:
        cm = p.get_cmap(cmap).with_extremes(under="k", bad="k")
    except ValueError:
        cm = p.get_cmap("bone").with_extremes(under="k", bad="k")
    fig, ax = p.subplots(figsize=figsize)
    ax.set_facecolor("k")
    beams = []
    ra0, de0 = parse_angle(src.center_ra, True), parse_angle(src.center_dec, False)
    if mom0 is not None:
        m0, mw, mb = load_image(mom0)
        if ra0 is None or de0 is None:
            ra0, de0 = mw.crval
        ext = _extent(mw, m0.shape, ra0, de0, src.center_offset_arcsec)
        label = r"integrated intensity (mJy beam$^{-1}$ km s$^{-1}$)"
    else:
        e, n, m0 = cube.sky_map(half_size, func="sum")
        st = cube.pixscale
        ext = [e[0] + st / 2, e[-1] - st / 2, n[0] - st / 2, n[-1] + st / 2]
        mb = (cube.bmaj, cube.bmin, cube.bpa)
        label = f"integrated intensity ({cube.units} km/s)"
    beams.append((mb, "w", "line"))
    vmax = np.nanpercentile(m0, 99.9)
    vmin = mom0_min if mom0_min is not None else max(np.nanpercentile(m0, 50), 0)
    im = ax.imshow(m0, origin="lower", cmap=cm, extent=ext, aspect="equal",
                   norm=PowerNorm(gamma, vmin=vmin, vmax=vmax), interpolation="nearest",
                   rasterized=True)
    if continuum is not None:
        c, cw, cb = load_image(continuum)
        if ra0 is None or de0 is None:
            ra0, de0 = cw.crval
        cmin = cont_min if cont_min is not None else np.nanpercentile(c, 90)
        c = np.where(c >= cmin, c, np.nan)
        ax.imshow(c ** cont_gamma, origin="lower", cmap="inferno", vmin=0, aspect="equal",
                  extent=_extent(cw, c.shape, ra0, de0, src.center_offset_arcsec),
                  interpolation="nearest", rasterized=True)
        beams.append((cb, "#ff8c1a", "cont."))
    cbar = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
    cbar.set_label(label, fontsize=8)
    ys = (np.atleast_1d(y_offsets) if y_offsets is not None else
          (_offsets(cfg.plot.pv_grid_offsets_arcsec)
           if cfg.plot.pv_grid_offsets_arcsec is not None else cfg.y_offsets))
    xh = cfg.grid.x_half_width_arcsec
    for y0 in ys:
        dE, dN = cube.xy_to_sky(np.array([-xh, xh]), np.array([y0, y0]))
        ax.plot(dE, dN, color="w", lw=0.6, alpha=0.8, zorder=4)
    for (Xa, Ya), col, lab in (((0.0, 1.6), "#ffdc00", r"$+Y$ (wind)"),
                               ((1.1, 0.0), "#2ecc40", r"$+X$")):
        dE, dN = cube.xy_to_sky(np.array([0.0, Xa]), np.array([0.0, Ya]))
        ax.annotate("", xy=(dE[1], dN[1]), xytext=(dE[0], dN[0]), zorder=7,
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.5, mutation_scale=10))
        u = np.array([dE[1] - dE[0], dN[1] - dN[0]])
        u /= np.hypot(*u)
        ax.text(dE[1] + 0.25 * u[0], dN[1] + 0.25 * u[1], lab, color=col, fontsize=9,
                ha="center", va="center", zorder=7, clip_on=True)
    pa = np.radians(src.disk_pa_deg)
    ax.plot([-1.2 * np.sin(pa), 1.2 * np.sin(pa)], [-1.2 * np.cos(pa), 1.2 * np.cos(pa)],
            ls="--", color="0.85", lw=0.8, zorder=6)
    ax.plot(0, 0, "+", color="w", ms=6, zorder=7)
    ax.set_xlim(half_size, -half_size)
    ax.set_ylim(-half_size, half_size)
    bx, by = half_size * 0.86, -half_size * 0.86
    for (bmaj, bmin, bpa), col, lab in beams:
        ax.add_patch(Ellipse((bx, by), width=bmin, height=bmaj, angle=-bpa, ec=col, fc=col,
                             zorder=10))
        ax.text(bx, by + 0.6 * max(bmaj, 0.25), lab, color=col, fontsize=7, ha="center",
                va="bottom", zorder=10)
        bx -= 0.26 * half_size
    ax.text(0.97, 0.03, src.name, transform=ax.transAxes, color="w", ha="right", fontsize=9)
    ax.set_xlabel(r"$\Delta\alpha\,\cos\delta$ ($^{\prime\prime}$)")
    ax.set_ylabel(r"$\Delta\delta$ ($^{\prime\prime}$)")
    ax.tick_params(top=False, right=False, direction="out")
    return finish(fig, save, 300)
