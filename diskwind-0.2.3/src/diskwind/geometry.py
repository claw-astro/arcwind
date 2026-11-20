"""
Viewing geometry and the projection equations of the conical-wind model.

Conventions (see :doc:`/concepts`):

* ``Y`` is the sky offset along the projected wind axis (toward the lobe,
  position angle ``PA_w``); ``X`` is perpendicular to it, toward ``PA_w - 90``.
* ``z`` is the height along the wind axis, ``phi`` the azimuth about ``+z``
  counted counter-clockwise from ``+X`` seen from the tip of the lobe.
* ``s = +1`` if the lobe axis is tilted toward the observer, ``-1`` if away.
* A point at height ``z`` and azimuth ``phi`` of a ring of radius ``R`` is seen at

  .. math::

     X = x_c + R\\cos\\phi, \\qquad Y = z\\sin i + sR\\sin\\phi\\cos i,

  with line-of-sight velocity

  .. math::

     v - V_{\\rm sys} = -sV_z\\cos i + V_r\\sin i\\sin\\phi + V_\\phi\\sin i\\cos\\phi .
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .constants import AU_CM, G_CGS, MSUN_G

__all__ = ["ViewGeometry", "cone_cut_heights", "los_velocity", "wind_to_sky",
           "keplerian_speed"]


@dataclass
class ViewGeometry:
    """How the wind is seen: inclination, lobe orientation, distance, Vsys.

    Parameters
    ----------
    inclination_deg : float
        Inclination of the wind axis to the line of sight (0-90 deg).
    s : int
        +1 lobe tilted toward the observer, -1 away.
    distance_pc : float
        Distance (1 arcsec = ``distance_pc`` au).
    vsys : float
        Systemic velocity (km/s).
    vsys_err : float
        Its uncertainty (km/s).
    disk_sense : int
        +1 if the disk rotates counter-clockwise about the wind axis.
    mstar_msun : float, optional
        Stellar mass.

    Examples
    --------
    >>> g = ViewGeometry(inclination_deg=49.0, s=-1, distance_pc=137.5, vsys=1.9)
    >>> round(float(g.zc_au(1.0)), 1)          # central height of the cut at y0 = 1"
    182.2
    >>> round(float(np.degrees(np.arctan(g.tan_theta_max))), 1)   # widest allowed cone
    47.5
    """
    inclination_deg: float
    s: int
    distance_pc: float
    vsys: float
    vsys_err: float = 0.1
    disk_sense: int = 1
    mstar_msun: Optional[float] = None

    @property
    def inc(self):
        """Inclination in radians."""
        return np.radians(self.inclination_deg)

    @property
    def sini(self):
        return float(np.sin(self.inc))

    @property
    def cosi(self):
        return float(np.cos(self.inc))

    @property
    def edge_on(self):
        """True if cos i < 0.05: Vz is then not constrained by the data."""
        return self.cosi < 0.05

    @property
    def tan_theta_max(self):
        """Largest cone opening (tan theta) for which every cut is a closed curve
        (0.95 tan i; 10 for edge-on views)."""
        return 0.95 * self.sini / max(self.cosi, 1e-6) if not self.edge_on else 10.0

    def zc_au(self, y0):
        """Height (au) of the ring centered on the cut at ``y0`` arcsec."""
        return np.asarray(y0) / self.sini * self.distance_pc

    def au(self, arcsec):
        return np.asarray(arcsec) * self.distance_pc

    def arcsec(self, au):
        return np.asarray(au) / self.distance_pc

    def vcent(self, Vz):
        """Velocity of a ring center, ``Vsys - s Vz cos i``."""
        return self.vsys - self.s * np.asarray(Vz) * self.cosi

    def mirrored(self):
        """The geometry of the opposite (counter) lobe.

        The counter lobe has its own wind frame: ``+Y`` along ``PA_w + 180``, so
        ``X`` is reversed, and its axis is tilted the other way (``s -> -s``).
        Seen about its own axis the disk rotates the other way
        (``disk_sense -> -disk_sense``). Applying it twice gives the original.
        """
        from dataclasses import replace
        return replace(self, s=-int(self.s), disk_sense=-int(self.disk_sense))

    def as_dict(self):
        return dict(inclination_deg=self.inclination_deg, s=int(self.s),
                    distance_pc=self.distance_pc, vsys=self.vsys, vsys_err=self.vsys_err,
                    disk_sense=int(self.disk_sense), mstar_msun=self.mstar_msun)


def cone_cut_heights(y0_au, R0, tan_theta, phi, geom):
    """Height at which a cone ``R = R0 + z tan(theta)`` crosses the cut ``Y = y0``.

    .. math::

       z(\\phi) = \\frac{y_0 - sR_0\\sin\\phi\\cos i}{\\sin i + s\\tan\\theta\\sin\\phi\\cos i}

    Parameters
    ----------
    y0_au : float
        Cut offset in au.
    R0, tan_theta : float
        Cone radius at z = 0 (au) and opening.
    phi : ndarray
        Azimuths (rad).
    geom : ViewGeometry

    Returns
    -------
    z, denom : ndarray
        Heights (au) and the denominator (the solution is valid where denom > 0).
    """
    sp = np.sin(phi)
    denom = geom.sini + geom.s * tan_theta * sp * geom.cosi
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (y0_au - geom.s * sp * geom.cosi * R0) / denom
    return z, denom


def los_velocity(Vz, Vr, Vphi, phi, geom):
    """Line-of-sight velocity ``v - Vsys`` of gas at azimuth ``phi``."""
    return (-geom.s * Vz * geom.cosi + Vr * geom.sini * np.sin(phi)
            + Vphi * geom.sini * np.cos(phi))


def wind_to_sky(X, Y, wind_pa_deg):
    """Wind-frame offsets (X, Y) -> sky offsets (dRA cos dec, dDec), same units."""
    p = np.radians(wind_pa_deg)
    return -X * np.cos(p) + Y * np.sin(p), X * np.sin(p) + Y * np.cos(p)


def keplerian_speed(r_au, mstar_msun, h_au=0.0):
    """Keplerian speed (km/s) at radius ``r_au`` and height ``h_au`` above the midplane."""
    r = np.maximum(r_au, 0.05)
    return np.sqrt(G_CGS * mstar_msun * MSUN_G * r ** 2 /
                   ((r ** 2 + h_au ** 2) ** 1.5 * AU_CM)) / 1e5
