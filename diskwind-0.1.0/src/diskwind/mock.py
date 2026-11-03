"""
Synthetic cubes with known parameters, to test the method and your settings.

The forward model is independent of the fitting code: particles are placed on
a conical shell and on a Keplerian disk (two emitting surfaces and the
midplane), projected, gridded, convolved with an elliptical beam, and noise is
added.

Examples
--------
>>> from diskwind.mock import MockWind, make_mock_cube
>>> truth = MockWind(theta_deg=30, vp=8.0, j=-40.0, approaching=True)
>>> path = make_mock_cube("mock.fits", truth, n_wind=300_000)       # doctest: +SKIP
>>> cfg = truth.config("mock.fits")      # a matching WindConfig      # doctest: +SKIP
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy.signal import fftconvolve

from .constants import C_KMS
from .io import write_fits

__all__ = ["MockWind", "make_mock_cube"]


@dataclass
class MockWind:
    """Parameters of a synthetic disk + conical wind.

    The wind has constant poloidal speed ``vp`` along a cone of half-angle
    ``theta_deg`` (``Vz = vp cos theta``, ``Vr = vp sin theta``) and specific
    angular momentum ``j`` (au km/s; sign: + = counter-clockwise about the axis).
    """
    inclination_deg: float = 49.0
    disk_pa_deg: float = 153.2
    wind_pa_deg: float = 63.2
    approaching: bool = True
    disk_sense: int = -1
    vsys: float = 3.0
    mstar_msun: float = 0.48
    distance_pc: float = 138.0
    R0_au: float = 8.0
    theta_deg: float = 30.0
    vp: float = 8.0
    j: float = -40.0
    zmax_au: float = 450.0
    shell_sigma_au: float = 2.0
    star_offset_arcsec: tuple = (0.05, -0.03)
    pix_arcsec: float = 0.03
    npix: int = 320
    beam: tuple = (0.12, 0.10, 80.0)
    dv: float = 0.35
    v0: float = -14.0
    nchan: int = 100
    rms_jy: float = 1.0e-3
    restfreq_hz: float = 230.538e9
    ra_deg: float = 246.5786169
    dec_deg: float = -24.4722560
    disk_rout_au: float = 200.0
    disk_z_over_r: float = 0.3

    @property
    def Vz(self):
        return self.vp * np.cos(np.radians(self.theta_deg))

    @property
    def Vr(self):
        return self.vp * np.sin(np.radians(self.theta_deg))

    def truth(self):
        """The true model parameters in the fitting code's names."""
        return dict(R0_au=self.R0_au, tan_theta=float(np.tan(np.radians(self.theta_deg))),
                    theta_deg=self.theta_deg, Vz=self.Vz, Vr=self.Vr, j=self.j, xoff_au=0.0,
                    xslope=0.0)

    def config(self, fits_file, **options):
        """A :class:`~diskwind.config.WindConfig` with this geometry."""
        from .config import WindConfig
        cfg = WindConfig.minimal(fits_file=fits_file, distance_pc=self.distance_pc,
                                 mstar_msun=self.mstar_msun,
                                 inclination_deg=self.inclination_deg,
                                 disk_pa_deg=self.disk_pa_deg, wind_pa_deg=self.wind_pa_deg)
        cfg.source.name = "mock"
        cfg.source.center_offset_arcsec = tuple(self.star_offset_arcsec)
        cfg.cuts.y_offsets_arcsec = dict(start=0.15, stop=2.0, step=0.1)
        cfg.grid.v_range_kms = (-12.0, 18.0)
        cfg.mask.v_exclude_kms = [(-1.0, 1.0)]
        return cfg.set(**options) if options else cfg


def _sky_from_wind(x, y, z, s, inc):
    ci, si = np.cos(np.radians(inc)), np.sin(np.radians(inc))
    return x, s * y * ci + z * si


def make_mock_cube(path, mock=None, n_wind=1_500_000, n_disk=2_000_000, seed=1):
    """Write a synthetic FITS cube.

    Parameters
    ----------
    path : str
        Output FITS file.
    mock : MockWind, optional
        Parameters (default: ``MockWind()``).
    n_wind, n_disk : int
        Number of particles (fewer = faster, noisier).
    seed : int

    Returns
    -------
    str
        ``path``.
    """
    m = mock or MockWind()
    rng = np.random.default_rng(seed)
    ci, si = np.cos(np.radians(m.inclination_deg)), np.sin(np.radians(m.inclination_deg))
    s = 1 if m.approaching else -1
    G, MS, AU = 6.674e-8, 1.989e33, 1.496e13
    z = m.zmax_au * rng.random(n_wind) ** 1.6 + 1.0
    phi = 2 * np.pi * rng.random(n_wind)
    R = m.R0_au + z * np.tan(np.radians(m.theta_deg)) + rng.normal(0, m.shell_sigma_au, n_wind)
    x, y = R * np.cos(phi), R * np.sin(phi)
    vphi = m.j / np.maximum(R, 1)
    vlos = -s * m.Vz * ci + m.Vr * si * np.sin(phi) + vphi * si * np.cos(phi)
    Xw, Yw = _sky_from_wind(x, y, z, s, m.inclination_deg)
    ww = (z / 50.0) ** 0.3 / R
    ww *= 30.0 / ww.sum() * n_wind / 1e6
    r = m.disk_rout_au * np.sqrt(rng.random(n_disk)) + 2.0
    ph = 2 * np.pi * rng.random(n_disk)
    h = m.disk_z_over_r * r * rng.choice([-1.0, 1.0, 0.0], n_disk, p=[0.4, 0.4, 0.2])
    vk = np.sqrt(G * m.mstar_msun * MS * r ** 2 / ((r ** 2 + h ** 2) ** 1.5 * AU)) / 1e5
    vlos_d = m.disk_sense * vk * si * np.cos(ph)
    Xd, Yd = _sky_from_wind(r * np.cos(ph), r * np.sin(ph), h, s, m.inclination_deg)
    wd = (r / 50.0) ** -0.5
    wd *= 60.0 / wd.sum() * n_disk / 1e6
    X = np.r_[Xw, Xd] / m.distance_pc
    Y = np.r_[Yw, Yd] / m.distance_pc
    V = np.r_[vlos, vlos_d] + m.vsys + rng.normal(0, 0.25, n_wind + n_disk)
    W = np.r_[ww, wd]
    p = np.radians(m.wind_pa_deg)
    dE = -X * np.cos(p) + Y * np.sin(p) + m.star_offset_arcsec[0]
    dN = X * np.sin(p) + Y * np.cos(p) + m.star_offset_arcsec[1]
    crpix = m.npix / 2 + 1
    px = (crpix - 1) - dE / m.pix_arcsec
    py = (crpix - 1) + dN / m.pix_arcsec
    vgrid = m.v0 + m.dv * np.arange(m.nchan)
    H, _ = np.histogramdd(np.c_[(V - m.v0) / m.dv, py, px], bins=(m.nchan, m.npix, m.npix),
                          range=((-0.5, m.nchan - 0.5), (-0.5, m.npix - 0.5),
                                 (-0.5, m.npix - 0.5)), weights=W)
    k = np.arange(-15, 16) * m.pix_arcsec
    KE, KN = np.meshgrid(-k, k)
    bmaj, bmin, bpa = m.beam
    b = np.radians(bpa)
    u = KE * np.sin(b) + KN * np.cos(b)
    w = KE * np.cos(b) - KN * np.sin(b)
    ker = np.exp(-0.5 * (u ** 2 / (bmaj / 2.3548) ** 2 + w ** 2 / (bmin / 2.3548) ** 2))
    cube = fftconvolve(H, ker[None], mode="same", axes=(1, 2))
    noise = fftconvolve(rng.standard_normal(cube.shape), ker[None], mode="same", axes=(1, 2))
    cube = (cube + noise * m.rms_jy / noise.std()).astype(np.float32)
    freqs = m.restfreq_hz * (1 - vgrid / C_KMS)
    cube, freqs = cube[::-1], freqs[::-1]
    hdr = {"SIMPLE": True, "BITPIX": -32, "NAXIS": 4, "NAXIS1": m.npix, "NAXIS2": m.npix,
           "NAXIS3": m.nchan, "NAXIS4": 1, "BUNIT": "Jy/beam", "BMAJ": bmaj / 3600,
           "BMIN": bmin / 3600, "BPA": bpa, "CTYPE1": "RA---SIN", "CRVAL1": m.ra_deg,
           "CDELT1": -m.pix_arcsec / 3600, "CRPIX1": crpix, "CUNIT1": "deg",
           "CTYPE2": "DEC--SIN", "CRVAL2": m.dec_deg, "CDELT2": m.pix_arcsec / 3600,
           "CRPIX2": crpix, "CUNIT2": "deg", "CTYPE3": "FREQ", "CRVAL3": freqs[0],
           "CDELT3": freqs[1] - freqs[0], "CRPIX3": 1.0, "CUNIT3": "Hz", "CTYPE4": "STOKES",
           "CRVAL4": 1.0, "CDELT4": 1.0, "CRPIX4": 1.0, "RESTFRQ": m.restfreq_hz,
           "SPECSYS": "LSRK"}
    write_fits(path, cube[None], hdr)
    return path
