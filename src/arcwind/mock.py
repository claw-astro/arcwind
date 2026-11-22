"""
Synthetic cubes with known parameters, to test the method and your settings.

The forward model is independent of the fitting code: particles are placed on
a conical shell and on a Keplerian disk (two emitting surfaces and the
midplane), projected, gridded, convolved with an elliptical beam, and noise is
added.

Examples
--------
>>> from arcwind.mock import MockWind, make_mock_cube
>>> truth = MockWind(theta_deg=30, vp=8.0, j=-40.0, approaching=True)
>>> path = make_mock_cube("mock.fits", truth, n_wind=300_000)       # doctest: +SKIP
>>> cfg = truth.config("mock.fits")      # a matching WindConfig      # doctest: +SKIP

Two presets:

* :meth:`MockWind.elias20_like` (the default): inclined (49 deg), one lobe;
* :meth:`MockWind.hh30_like`: HH 30-like -- nearly edge-on, a wide cone with a
  constant poloidal speed, slow rotation with the disk, and **both lobes**
  (``lobes="both"``; the real CO outflow is one-sided: use ``lobes="main"``).

>>> hh30 = MockWind.hh30_like()
>>> hh30.inclination_deg, hh30.theta_deg, hh30.vp, hh30.lobes
(89.0, 35.0, 9.3, 'both')
"""
from __future__ import annotations

from dataclasses import dataclass

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
    angular momentum ``j`` (au km/s; sign: + = counter-clockwise about the axis of
    the main lobe, like ``disk_sense``).

    ``lobes = "main"`` makes one lobe (along ``wind_pa_deg``); ``"both"`` adds the
    mirror-symmetric counter lobe (along ``wind_pa_deg + 180``), with the same
    cone, speeds and physical sense of rotation.
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
    lobes: str = "main"
    name: str = "mock"

    @classmethod
    def elias20_like(cls, **changes):
        """Inclined (i = 49 deg) one-sided wind with the geometry of Elias 2-20
        (the default mock)."""
        return cls(**changes)

    @classmethod
    def hh30_like(cls, **changes):
        """A disk and bipolar wind like HH 30 (Louvet et al. 2018).

        Distance 140 pc, M* = 0.45 Msun, Vsys = 6.9 km/s; the wind axis at PA 31.2
        deg, inclined by 91 deg to the line of sight (here: i = 89 deg with the
        main lobe tilted away, the same geometry); a cone of half-angle 35 deg from
        R0 = 20 au (<= 22 au), constant poloidal speed 9.3 km/s, specific angular
        momentum 38 au km/s in the sense of the disk. The beam
        (0.26" x 0.19", PA 29 deg) and channels (0.3 km/s) are those of the ALMA
        data; the CO disk extends to 150 au.

        The CO outflow of HH 30 is monopolar (north-east lobe only); this mock
        adds a mirror-image counter lobe (``lobes="both"``) to demonstrate fitting
        two lobes. Use ``MockWind.hh30_like(lobes="main")`` for one lobe.
        """
        p = dict(inclination_deg=89.0, disk_pa_deg=121.2, wind_pa_deg=31.2,
                 approaching=False, disk_sense=+1, vsys=6.9, mstar_msun=0.45,
                 distance_pc=140.0, R0_au=20.0, theta_deg=35.0, vp=9.3, j=38.0,
                 zmax_au=500.0, shell_sigma_au=3.0, star_offset_arcsec=(0.03, -0.02),
                 pix_arcsec=0.04, npix=256, beam=(0.26, 0.19, 29.0), dv=0.3, v0=-5.1,
                 nchan=80, rms_jy=1.0e-3, ra_deg=67.906250, dec_deg=18.206611,
                 disk_rout_au=150.0, disk_z_over_r=0.2, lobes="both", name="HH 30-like mock")
        p.update(changes)
        return cls(**p)

    @classmethod
    def tall_ring_like(cls, **changes):
        """A narrow, fast cone seen at i = 60 deg: tall, narrow rings in the PV
        diagrams (velocity extent >> spatial extent), whose steep walls need
        ``trace.mode = "ring"``.

        Geometry and star loosely follow CW Tau (131.6 pc, 0.71 Msun, PA 240 deg,
        Vsys 6.3 km/s); a cone of half-angle 18 deg from R0 = 5 au, poloidal speed
        24 km/s and j = -120 au km/s, observed with a 0.28" x 0.23" beam.
        """
        p = dict(inclination_deg=60.0, disk_pa_deg=240.0, wind_pa_deg=150.0,
                 approaching=True, disk_sense=-1, vsys=6.3, mstar_msun=0.71,
                 distance_pc=131.6, R0_au=5.0, theta_deg=18.0, vp=24.0, j=-120.0,
                 zmax_au=400.0, shell_sigma_au=3.0, pix_arcsec=0.04, npix=200,
                 beam=(0.28, 0.23, 0.0), dv=0.32, v0=-20.0, nchan=110, rms_jy=1.0e-3,
                 disk_rout_au=120.0, name="tall-ring mock")
        p.update(changes)
        return cls(**p)

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
        """A :class:`~arcwind.config.WindConfig` with this geometry."""
        from .config import WindConfig
        cfg = WindConfig.minimal(fits_file=fits_file, distance_pc=self.distance_pc,
                                 mstar_msun=self.mstar_msun,
                                 inclination_deg=self.inclination_deg,
                                 disk_pa_deg=self.disk_pa_deg, wind_pa_deg=self.wind_pa_deg)
        cfg.source.name = self.name
        cfg.source.center_offset_arcsec = tuple(self.star_offset_arcsec)
        cfg.source.lobes = self.lobes
        if self.inclination_deg > 80:
            # nearly edge-on: the tilt and the disk sense are hard to estimate,
            # so state them (as one would from the jet or scattered-light images)
            cfg.source.wind_approaching = bool(self.approaching)
            cfg.source.disk_rotation_sense = int(self.disk_sense)
            cfg.source.vsys_kms = float(self.vsys)
            zc = 0.95 * self.zmax_au * np.sin(np.radians(self.inclination_deg))
            ymax = round(min(zc / self.distance_pc, 0.45 * self.npix * self.pix_arcsec), 1)
            cfg.cuts.y_offsets_arcsec = dict(start=0.4, stop=ymax, step=0.2)
            cfg.grid.v_range_kms = (self.v0 - self.dv / 2, self.v0 + self.dv * (self.nchan - 0.5))
            cfg.grid.y_range_arcsec = (-1.0, ymax + 0.5)
            cfg.mask.disk_rout_au = self.disk_rout_au
            cfg.mask.disk_z_over_r = self.disk_z_over_r
        else:
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
    # lobe of each particle: +1 main (along +z), -1 counter (along -z). A separate
    # generator keeps the one-lobe cube identical to earlier versions.
    if m.lobes == "both":
        lobe = np.where(np.random.default_rng(seed + 7).random(n_wind) < 0.5, 1.0, -1.0)
    elif m.lobes == "main":
        lobe = np.ones(n_wind)
    else:
        raise ValueError("MockWind.lobes must be 'main' or 'both'")
    # velocity components in the frame of the main lobe: vz along +z, vy along y
    vz = lobe * m.Vz
    vy = m.Vr * np.sin(phi) + vphi * np.cos(phi)
    vlos = -s * vz * ci + vy * si
    Xw, Yw = _sky_from_wind(x, y, lobe * z, s, m.inclination_deg)
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
