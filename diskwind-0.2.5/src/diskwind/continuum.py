"""
Disk center (and orientation) from a continuum image.

The disk center sets the origin of the wind frame, so it should come from the
continuum, which traces the disk midplane. :func:`fit_disk_center` fits an
elliptical 2-D Gaussian to the continuum around its peak and returns the
center as absolute coordinates and as an offset from the image reference pixel
(``CRVAL``), ready for the ``[source]`` block of a configuration.

For disks with rings, an ellipse fitted to a ring (as in Huang et al. 2018 for
the DSHARP disks) is more precise; enter that center directly
(``center_offset_arcsec``).

Examples
--------
>>> from diskwind.continuum import fit_disk_center
>>> c = fit_disk_center("Elias20_continuum.fits")          # doctest: +SKIP
>>> c.offset_arcsec, c.offset_err_arcsec                    # doctest: +SKIP
>>> print(c.toml())                 # lines to paste into [source]   # doctest: +SKIP

or ``diskwind center Elias20_continuum.fits`` on the command line, or
``center_from_continuum = true`` with ``continuum_fits`` in the configuration.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from .io import CelestialWCS, read_fits

__all__ = ["fit_disk_center", "ContinuumCenter"]


@dataclass
class ContinuumCenter:
    """Result of :func:`fit_disk_center`.

    Attributes
    ----------
    ra_deg, dec_deg : float
        Disk center (deg).
    offset_arcsec, offset_err_arcsec : tuple
        Center relative to the reference pixel (CRVAL) of the continuum image, as
        (dRA cos dec, dDec), and its uncertainty (arcsec).
    fwhm_major_arcsec, fwhm_minor_arcsec, pa_deg : float
        Size and position angle (E of N) of the fitted Gaussian.
    inclination_deg : float
        arccos(minor / major) -- a rough estimate (beam-convolved, thick disk).
    peak : float
        Peak of the Gaussian (image units).
    """
    ra_deg: float
    dec_deg: float
    offset_arcsec: tuple
    offset_err_arcsec: tuple
    fwhm_major_arcsec: float
    fwhm_minor_arcsec: float
    pa_deg: float
    inclination_deg: float
    peak: float
    path: str = ""

    def toml(self):
        """Lines for the ``[source]`` block: the absolute center, so it is right for
        any image of the source (cubes and continuum may have different reference
        pixels)."""
        e = self.offset_err_arcsec
        return (f"center_ra = {self.ra_deg:.8f}      # disk center from the continuum fit\n"
                f"center_dec = {self.dec_deg:.8f}    # (+/- {e[0] * 1e3:.1f}, "
                f"{e[1] * 1e3:.1f} mas)\n"
                f"center_offset_arcsec = [0.0, 0.0]")

    def summary(self):
        o, e = self.offset_arcsec, self.offset_err_arcsec
        return (f"continuum center: RA = {self.ra_deg:.7f}, Dec = {self.dec_deg:.7f} deg\n"
                f"offset from the reference pixel: dRA cos dec = {o[0] * 1e3:.1f} +/- "
                f"{e[0] * 1e3:.1f} mas, dDec = {o[1] * 1e3:.1f} +/- {e[1] * 1e3:.1f} mas\n"
                f"Gaussian FWHM {self.fwhm_major_arcsec:.3f}\" x {self.fwhm_minor_arcsec:.3f}\", "
                f"PA {self.pa_deg:.1f} deg (inclination ~ {self.inclination_deg:.0f} deg, rough)")


def _gauss2d(p, X, Y):
    a, x0, y0, sx, sy, th, c = p
    ct, st = np.cos(th), np.sin(th)
    u = (X - x0) * ct + (Y - y0) * st
    v = -(X - x0) * st + (Y - y0) * ct
    return a * np.exp(-0.5 * ((u / sx) ** 2 + (v / sy) ** 2)) + c


def fit_disk_center(path, box_arcsec=None, guess_offset_arcsec=None):
    """Fit an elliptical 2-D Gaussian to a continuum image.

    Parameters
    ----------
    path : str
        Continuum FITS image (any extra axes of length 1 are dropped).
    box_arcsec : float, optional
        Half-size of the fitted region around the peak; default: 3 x the
        half-width at 10% of the peak, at least 5 beams.
    guess_offset_arcsec : tuple, optional
        Start the search here (offset from CRVAL) instead of at the brightest pixel.

    Returns
    -------
    ContinuumCenter
    """
    hdr, data, _ = read_fits(path)
    img = np.asarray(data, float)
    while img.ndim > 2:
        img = img[0]
    wcs = CelestialWCS(hdr, 1, 2)
    pix = wcs.pixscale_arcsec
    ny, nx = img.shape
    ra0, de0 = wcs.crval
    cx, cy = wcs.world_to_pixel0(ra0, de0)          # reference pixel (0-based)
    # offsets (arcsec) of every pixel from the reference pixel, east positive
    J, I = np.mgrid[0:ny, 0:nx]
    d = np.tensordot(wcs.cd, np.array([I - cx, J - cy], float), axes=1) * 3600.0
    E, N = d[0], d[1]
    good = np.isfinite(img)
    if guess_offset_arcsec is None:
        k = np.nanargmax(np.where(good, img, -np.inf))
        e0, n0 = E.flat[k], N.flat[k]
    else:
        e0, n0 = guess_offset_arcsec
    peak = np.nanmax(img)
    bmaj = hdr.get("BMAJ", 3 * pix / 3600) * 3600
    if box_arcsec is None:
        bright = good & (img > 0.1 * peak) & (np.hypot(E - e0, N - n0) < 50 * bmaj)
        r10 = np.max(np.hypot(E[bright] - e0, N[bright] - n0)) if bright.any() else bmaj
        box_arcsec = max(3 * r10, 5 * bmaj)
    sel = good & (np.abs(E - e0) <= box_arcsec) & (np.abs(N - n0) <= box_arcsec)
    x, y, z = E[sel], N[sel], img[sel]
    w = np.clip(z, 0, None)
    xm, ym = np.sum(w * x) / np.sum(w), np.sum(w * y) / np.sum(w)
    sx0 = np.sqrt(np.sum(w * (x - xm) ** 2) / np.sum(w)) + 0.5 * bmaj
    p0 = [peak, xm, ym, sx0, sx0, 0.0, np.median(z)]
    res = least_squares(lambda p: _gauss2d(p, x, y) - z, p0,
                        bounds=([0, xm - box_arcsec, ym - box_arcsec, pix / 3, pix / 3, -np.pi,
                                 -np.inf],
                                [np.inf, xm + box_arcsec, ym + box_arcsec, 10 * box_arcsec,
                                 10 * box_arcsec, np.pi, np.inf]))
    a, xc, yc, sx, sy, th, c = res.x
    r = res.fun
    dof = max(len(z) - 7, 1)
    cov = np.linalg.pinv(res.jac.T @ res.jac) * np.sum(r ** 2) / dof
    # pixels are correlated over a beam: scale the errors by sqrt(pixels per beam)
    bmin = hdr.get("BMIN", hdr.get("BMAJ", pix / 3600)) * 3600
    npb = max(1.1331 * bmaj * bmin / pix ** 2, 1.0)
    err = np.sqrt(np.clip(np.diag(cov), 0, None) * npb)
    if sx < sy:
        sx, sy, th = sy, sx, th + np.pi / 2
    # theta is measured from +E toward +N; PA from N toward E
    pa = (90.0 - np.degrees(th)) % 180.0
    ra = ra0 + xc / 3600.0 / np.cos(np.radians(de0))
    dec = de0 + yc / 3600.0
    f = 2 * np.sqrt(2 * np.log(2))
    return ContinuumCenter(ra_deg=float(ra), dec_deg=float(dec), offset_arcsec=(float(xc), float(yc)),
                           offset_err_arcsec=(float(err[1]), float(err[2])),
                           fwhm_major_arcsec=float(f * sx), fwhm_minor_arcsec=float(f * sy),
                           pa_deg=float(pa), inclination_deg=float(np.degrees(np.arccos(sy / sx))),
                           peak=float(a), path=str(path))
