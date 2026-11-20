"""
Contamination masks: which pixels of a PV diagram the trace may use.

Three masks are combined (each can be switched off in ``[mask]``):

1. a **velocity window** ``vrel_range_kms`` and excluded ranges ``v_exclude_kms``
   (relative to Vsys), e.g. to avoid cloud absorption;
2. a **spatial window** in X: ``x_cone`` (follows the opening of the wind),
   ``x_range_arcsec`` (fixed) or ``|X| <= x_max_arcsec``;
3. a **Keplerian disk mask** (``kepler = true``): channels where disk gas between
   the emitting surfaces ``|z| <= disk_z_over_r * r`` (``r <= disk_rout_au``) can
   emit anywhere along the line of sight within the beam, padded by
   ``kepler_dv_kms``.

Examples
--------
>>> from diskwind.masks import ContaminationMask
>>> m = ContaminationMask(cube)                                  # doctest: +SKIP
>>> parts = m.components(1.0)    # dict: 'velocity', 'spatial', 'kepler'   # doctest: +SKIP
>>> excluded = m.excluded(1.0)   # their union, shape (nv, nx)    # doctest: +SKIP
"""
from __future__ import annotations

import warnings

import numpy as np

from .geometry import keplerian_speed

__all__ = ["ContaminationMask", "x_window"]


def x_window(cfg, y0):
    """Allowed X range (arcsec) at ``y0`` for a (per-cut) configuration ``cfg``.

    Priority: ``mask.x_cone`` > ``mask.x_range_arcsec`` > ``|X| <= mask.x_max_arcsec``
    (a per-cut ``x_range_arcsec`` override replaces ``x_cone``, see
    :meth:`WindConfig.for_cut`).
    """
    m = cfg.mask
    if m.x_cone is not None:
        xc, w0, slope = m.x_cone[:3]
        xcs = m.x_cone[3] if len(m.x_cone) > 3 else 0.0
        hw = w0 + slope * max(y0, 0.0)
        xc = xc + xcs * y0
        return (xc - hw, xc + hw)
    if m.x_range_arcsec is not None:
        return tuple(m.x_range_arcsec)
    return (-m.x_max_arcsec, m.x_max_arcsec)


class ContaminationMask:
    """Builds the masks of the cuts of a :class:`~diskwind.cube.WindCube`."""

    def __init__(self, cube):
        self.cube = cube

    def kepler_band(self, y0, x=None, width=None, nlos=301):
        """Range of ``v - Vsys`` of Keplerian disk gas along each line of sight.

        For each X of the cut at ``y0`` the line of sight is sampled at nine
        points spanning the beam and the cut width; points inside the disk volume
        get the projected Keplerian velocity.

        Returns
        -------
        vmin, vmax : ndarray
            NaN where no disk gas lies on the line of sight.
        """
        c = self.cube
        cfg = c.config
        x = c.x if x is None else x
        w = c.slit_width if width is None else width
        bx, by = 0.5 * c.beam_x, 0.5 * (c.beam_y + w)
        offs = [(dx, dy) for dx in (-bx, 0.0, bx) for dy in (-by, 0.0, by)]
        rout = cfg.mask.disk_rout_au
        W = np.linspace(-1.05 * rout, 1.05 * rout, nlos)
        vmin = np.full(len(x), np.inf)
        vmax = np.full(len(x), -np.inf)
        for dx, dy in offs:
            Xa = ((x + dx) * c.dist)[:, None]
            Ya = (y0 + dy) * c.dist
            yy = Ya * c.s * c.cosi - W[None, :] * c.sini
            zz = Ya * c.sini + W[None, :] * c.s * c.cosi
            r = np.hypot(Xa, yy)
            inside = (r <= rout) & (np.abs(zz) <= c.disk_height(r))
            v = c.dsense * keplerian_speed(r, c.mstar, zz) * c.sini * Xa / np.maximum(r, 1e-3)
            v = np.where(inside, v, np.nan)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                vmin = np.fmin(vmin, np.nanmin(v, axis=1))
                vmax = np.fmax(vmax, np.nanmax(v, axis=1))
        vmin[~np.isfinite(vmin)] = np.nan
        vmax[~np.isfinite(vmax)] = np.nan
        return vmin, vmax

    def components(self, y0, width=None):
        """The individual masks of the cut at ``y0`` (True = excluded).

        Returns
        -------
        dict
            ``velocity``, ``spatial`` and ``kepler`` boolean arrays of shape (nv, nx).
        """
        c = self.cube
        cfg = c.cut_config(y0)
        vrel = (c.v - c.vsys)[:, None]
        shape = (len(c.v), len(c.x))
        mv = np.zeros(shape, bool)
        for lo, hi in cfg.mask.v_exclude_kms or []:
            mv |= (vrel >= lo) & (vrel <= hi)
        vr = cfg.mask.vrel_range_kms
        if vr is not None:
            mv |= (vrel < vr[0]) | (vrel > vr[1])
        xlo, xhi = x_window(cfg, y0)
        ms = np.broadcast_to(((c.x < xlo) | (c.x > xhi))[None, :], shape).copy()
        mk = np.zeros(shape, bool)
        if cfg.mask.kepler:
            vmin, vmax = self.kepler_band(y0, width=width)
            dv = cfg.mask.kepler_dv_kms
            ok = np.isfinite(vmin)
            mk = ok[None, :] & (vrel >= (vmin - dv)[None, :]) & (vrel <= (vmax + dv)[None, :])
        return dict(velocity=mv, spatial=ms, kepler=mk)

    def excluded(self, y0, width=None):
        """Union of all masks of the cut at ``y0`` (True = excluded)."""
        comp = self.components(y0, width=width)
        return comp["velocity"] | comp["spatial"] | comp["kepler"]
