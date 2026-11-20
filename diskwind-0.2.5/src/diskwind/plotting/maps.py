"""Sky maps: the layout of the cuts on the moment-0 / continuum images."""
from __future__ import annotations

import warnings

import numpy as np

from ..io import CelestialWCS, parse_angle, read_fits
from ._style import finish, plt

__all__ = ["plot_slit_layout", "load_image"]


def image_scale(bunit):
    """Factor converting an image in ``bunit`` to Jy/beam (continuum) or
    Jy/beam km/s (moment 0), and whether the unit includes a velocity.

    >>> image_scale("Jy/beam m/s")
    (0.001, True)
    >>> image_scale("mJy/beam km/s")
    (0.001, True)
    >>> image_scale("Jy/beam")
    (1.0, False)
    """
    u = (bunit or "Jy/beam").replace(" ", "").lower()
    flux = 1e-3 if u.startswith("mjy") else (1e-6 if u.startswith("ujy") else 1.0)
    if "km/s" in u or "kms-1" in u or "km.s-1" in u:
        return flux, True
    if "m/s" in u or "ms-1" in u or "m.s-1" in u:
        return flux * 1e-3, True
    return flux, False


def load_image(path):
    """2-D image converted to Jy/beam (continuum) or Jy/beam km/s (moment 0, from
    the BUNIT of the file, e.g. 'Jy/beam m/s'), its WCS and beam (arcsec, arcsec,
    deg)."""
    hdr, data, beams = read_fits(path)
    img = np.asarray(data, float)
    while img.ndim > 2:
        img = img[0]
    wcs = CelestialWCS(hdr, 1, 2)
    if beams is not None and len(beams):
        bmaj, bmin, bpa = np.median(np.array([list(r)[:3] for r in beams], float), axis=0)
    else:
        bmaj, bmin, bpa = hdr["BMAJ"] * 3600.0, hdr["BMIN"] * 3600.0, hdr["BPA"]
    bunit = str(hdr.get("BUNIT", "Jy/beam"))
    if "jy" not in bunit.lower():
        warnings.warn(f"{path}: BUNIT = {bunit!r}, assumed Jy/beam")
    fac, _ = image_scale(bunit)
    return img * fac, wcs, (bmaj, bmin, bpa)


def _extent(wcs, shape, ra0, de0, off):
    if abs(wcs.cd[0, 1]) > 1e-3 * abs(wcs.cd[0, 0]):
        warnings.warn("image axes are rotated w.r.t. RA/Dec; extent is approximate")
    px, py = wcs.world_to_pixel0(ra0, de0)
    pc = np.array([px, py]) + wcs.cdinv @ (np.asarray(off, float) / 3600.0)
    ny, nx = shape
    return [wcs.cd[0, 0] * (-0.5 - pc[0]) * 3600, wcs.cd[0, 0] * (nx - 0.5 - pc[0]) * 3600,
            wcs.cd[1, 1] * (-0.5 - pc[1]) * 3600, wcs.cd[1, 1] * (ny - 0.5 - pc[1]) * 3600]


def _get_cmap(name):
    """A matplotlib color map by name; 'cmr.*' names load cmasher if installed."""
    import matplotlib
    if name.startswith("cmr."):
        try:
            import cmasher  # noqa: F401  (registers the cmr.* maps)
        except ImportError:
            warnings.warn(f"{name} needs the cmasher package (pip install cmasher); "
                          "using 'bone'")
            name = "bone"
    try:
        return matplotlib.colormaps[name]
    except KeyError:
        warnings.warn(f"unknown color map {name!r}; using 'bone'")
        return matplotlib.colormaps["bone"]


def plot_slit_layout(cube, y_offsets=None, mom0=None, continuum=None, half_size=3.5,
                     gamma=0.75, mom0_min=None, cmap=None, cont_min=None, cont_gamma=0.5,
                     line_label=None, figsize=(5.6, 4.8), save=None):
    """The transverse cuts on the sky, with the wind-frame axes and the beams.

    Parameters
    ----------
    cube : WindCube
    y_offsets : array-like, optional
        Cuts to draw (default: ``plot.pv_grid_offsets_arcsec`` or the fit cuts).
    mom0 : str, optional
        Moment-0 FITS file; default: ``plot.mom0_fits``, else the integrated
        intensity of the cube.
    continuum : str or False, optional
        Continuum FITS file drawn on top (inferno); default:
        ``source.continuum_fits``. ``False``: no continuum.
    half_size : float
        Half-width of the field (arcsec).
    gamma, mom0_min, cmap : float, float, str
        Stretch, floor and color map of the moment-0 image (defaults:
        ``plot.map_min``, ``plot.map_cmap`` = cmasher 'cmr.freeze').
    cont_min, cont_gamma : float
        Floor (default ``plot.continuum_min``) and stretch of the continuum.
    line_label : str, optional
        Name of the line, written on the map and the color bar (default
        ``plot.line_label``).

    Returns
    -------
    Figure

    Examples
    --------
    >>> plot_slit_layout(cube, mom0="Elias20_CO_M0.fits",
    ...                  continuum="Elias20_continuum.fits",
    ...                  line_label=r"$^{12}$CO $J$=2$-$1")       # doctest: +SKIP
    """
    from matplotlib.colors import PowerNorm
    from matplotlib.patches import Ellipse
    from ..config import _offsets
    p = plt()
    cfg = cube.config
    src, pc = cfg.source, cfg.plot
    mom0 = cfg.resolve(mom0 if mom0 is not None else pc.mom0_fits)
    if continuum is None:
        continuum = cfg.resolve(src.continuum_fits)
    elif continuum is False:
        continuum = None
    import os
    for name, f in (("moment-0", mom0), ("continuum", continuum)):
        if f is not None and not os.path.exists(f):
            warnings.warn(f"{name} image {f} not found: not used")
    if mom0 is not None and not os.path.exists(mom0):
        mom0 = None
    if continuum is not None and not os.path.exists(continuum):
        continuum = None
    mom0_min = pc.map_min if mom0_min is None else mom0_min
    cont_min = pc.continuum_min if cont_min is None else cont_min
    line_label = pc.line_label if line_label is None else line_label
    cm = _get_cmap(cmap or pc.map_cmap).with_extremes(under="k", bad="k")
    fig, ax = p.subplots(figsize=figsize)
    ax.set_facecolor("k")
    beams = []
    ra0, de0 = parse_angle(src.center_ra, True), parse_angle(src.center_dec, False)
    if mom0 is not None:
        m0, mw, mb = load_image(mom0)
        m0 = m0 * 1e3                               # Jy -> mJy/beam km/s
        if ra0 is None or de0 is None:
            ra0, de0 = mw.crval
        ext = _extent(mw, m0.shape, ra0, de0, src.center_offset_arcsec)
        label = r"(mJy beam$^{-1}$ km s$^{-1}$)"
    else:
        e, n, m0 = cube.sky_map(half_size, func="sum")
        if cube.units != "K":
            m0 = m0 / cube.ufac * 1e3               # -> mJy/beam km/s
        st = cube.pixscale
        ext = [e[0] + st / 2, e[-1] - st / 2, n[0] - st / 2, n[-1] + st / 2]
        mb = (cube.bmaj, cube.bmin, cube.bpa)
        label = (r"(K km s$^{-1}$)" if cube.units == "K"
                 else r"(mJy beam$^{-1}$ km s$^{-1}$)")
    beams.append((mb, "w", "line"))
    label = (f"{line_label} integrated intensity " if line_label else "integrated intensity ") \
        + label
    # an image from the cube is dominated by the disk: saturate it so the wind shows
    vmax = np.nanpercentile(m0, 99.9 if mom0 is not None else 98.0)
    vmin = mom0_min if mom0_min is not None else max(np.nanpercentile(m0, 50), 0)
    im = ax.imshow(m0, origin="lower", cmap=cm, extent=ext, aspect="equal",
                   norm=PowerNorm(gamma, vmin=vmin, vmax=vmax), interpolation="nearest",
                   rasterized=True)
    if continuum is not None:
        c, cw, cb = load_image(continuum)
        c = c * 1e3                                 # Jy -> mJy/beam
        if ra0 is None or de0 is None:
            ra0, de0 = cw.crval
        cmin = cont_min if cont_min is not None else np.nanpercentile(c, 90)
        c = np.where(c >= cmin, c, np.nan)
        ax.imshow(c ** cont_gamma, origin="lower", cmap="inferno", vmin=0, aspect="equal",
                  extent=_extent(cw, c.shape, ra0, de0, src.center_offset_arcsec),
                  interpolation="nearest", rasterized=True)
        beams.append((cb, "#ff8c1a", "cont."))
    cbar = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
    cbar.set_label(label, fontsize="x-small")
    ys = (np.atleast_1d(y_offsets) if y_offsets is not None else
          (_offsets(cfg.plot.pv_grid_offsets_arcsec)
           if cfg.plot.pv_grid_offsets_arcsec is not None else cfg.y_offsets))
    xh = cfg.grid.x_half_width_arcsec
    both = cfg.source.lobes == "both"
    for y0 in ys:
        for sgn in ((1, -1) if both else (1,)):     # with both lobes: cuts on both sides
            dE, dN = cube.xy_to_sky(np.array([-xh, xh]), np.array([sgn * y0, sgn * y0]))
            ax.plot(dE, dN, color="w", lw=0.6, alpha=0.8, zorder=4)
    arrows = [((0.0, 1.6), "#ffdc00", r"$+Y$ (wind)"), ((1.1, 0.0), "#2ecc40", r"$+X$")]
    if both:
        arrows.append(((0.0, -1.6), "#ffdc00", "counter lobe"))
    for (Xa, Ya), col, lab in arrows:
        dE, dN = cube.xy_to_sky(np.array([0.0, Xa]), np.array([0.0, Ya]))
        ax.annotate("", xy=(dE[1], dN[1]), xytext=(dE[0], dN[0]), zorder=7,
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.5, mutation_scale=10))
        u = np.array([dE[1] - dE[0], dN[1] - dN[0]])
        u /= np.hypot(*u)
        ax.text(dE[1] + 0.25 * u[0], dN[1] + 0.25 * u[1], lab, color=col, fontsize="small",
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
        ax.text(bx, by + 0.6 * max(bmaj, 0.25), lab, color=col, fontsize="x-small", ha="center",
                va="bottom", zorder=10)
        bx -= 0.26 * half_size
    ax.text(0.97, 0.03, src.name, transform=ax.transAxes, color="w", ha="right", fontsize="small")
    ytxt = 0.96
    if line_label:
        ax.text(0.04, ytxt, line_label, transform=ax.transAxes, color=cm(0.85), ha="left",
                va="top", fontsize="small", zorder=11)
        ytxt -= 0.065
    if continuum is not None:
        ax.text(0.04, ytxt, "continuum", transform=ax.transAxes, color="#ff8c1a", ha="left",
                va="top", fontsize="small", zorder=11)
    ax.set_xlabel(r"$\Delta\alpha\,\cos\delta$ ($^{\prime\prime}$)")
    ax.set_ylabel(r"$\Delta\delta$ ($^{\prime\prime}$)")
    ax.tick_params(top=False, right=False, direction="out")
    return finish(fig, save, 300)
