"""Transverse position-velocity diagrams."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

__all__ = ["PVDiagram", "contour_levels"]


@dataclass
class PVDiagram:
    """A transverse PV diagram (cut perpendicular to the wind axis).

    Attributes
    ----------
    y0 : float
        Offset of the cut along the wind axis (arcsec).
    x : ndarray, shape (nx,)
        Position along the cut (arcsec, +X toward PA_w - 90).
    v : ndarray, shape (nv,)
        LSRK velocity of each channel (km/s).
    data : ndarray, shape (nv, nx)
        Intensity (``units``).
    rms : float
        Noise of this diagram, measured from noise channels processed like the
        data (same resampling, slit averaging and smoothing).
    width : float
        Width of the cut along Y (arcsec).
    units : str
    zc_au : float
        Height of the ring centered on the cut (au).
    excluded : ndarray of bool, shape (nv, nx), optional
        PV pixels not used by the trace (masks). Filled by
        :meth:`WindCube.pv` when ``with_mask=True``.

    Examples
    --------
    >>> pv = cube.pv(1.0)                                    # doctest: +SKIP
    >>> pv.data.shape, pv.rms                                # doctest: +SKIP
    ((86, 201), 0.97)
    >>> spec = pv.spectrum(0.3)          # spectrum at X = 0.3"   # doctest: +SKIP
    >>> pv.plot(contours=True)                               # doctest: +SKIP
    """
    y0: float
    x: np.ndarray
    v: np.ndarray
    data: np.ndarray
    rms: float
    width: float
    units: str = "mJy/beam"
    zc_au: float = np.nan
    excluded: Optional[np.ndarray] = None
    smoothing: tuple = field(default=(1, 1))

    @property
    def extent(self):
        """``imshow`` extent ``[x0, x1, v0, v1]``."""
        return [self.x[0], self.x[-1], self.v[0], self.v[-1]]

    @property
    def snr(self):
        """Data in units of the rms."""
        return self.data / self.rms

    def spectrum(self, x):
        """Spectrum at the position nearest to ``x`` (arcsec)."""
        return self.data[:, int(np.argmin(np.abs(self.x - x)))]

    def plot(self, ax=None, contours=True, show_mask=False, cmap="magma", **kw):
        """Show the diagram (see :func:`diskwind.plotting.plot_pv`)."""
        from .plotting import plot_pv
        return plot_pv(self, ax=ax, contours=contours, show_mask=show_mask, cmap=cmap, **kw)


def contour_levels(rms, vmax, start=5.0, step=5.0):
    """Contours from ``start`` x rms in steps of ``step`` x rms up to ``vmax``."""
    if not np.isfinite(vmax) or vmax < start * rms:
        return np.array([start * rms])
    return np.arange(start * rms, vmax + step * rms, step * rms)
