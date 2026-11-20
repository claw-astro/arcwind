"""
Checks to run before fitting.

* :func:`check_noise` -- are the noise channels line-free? (per-channel rms,
  peak S/N, fraction of pixels above 4 sigma; per-cut PV rms from the blue end,
  the red end and the current channels).
* :func:`check_disk_rotation` -- which sense of disk rotation fits the data?
  (emission along projected Keplerian curves for both senses, in cuts parallel
  to the major axis on the side without wind).

Examples
--------
>>> rep = cube.check_noise()                  # doctest: +SKIP
>>> print(rep.summary())                      # doctest: +SKIP
>>> rep.plot()                                # doctest: +SKIP
>>> chk = cube.check_disk_rotation()          # doctest: +SKIP
>>> chk.best_sense, chk.scores                # doctest: +SKIP
(1, {1: 5230.0, -1: 812.0})
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import map_coordinates

from .cube import robust_std
from .geometry import keplerian_speed

__all__ = ["check_noise", "NoiseReport", "check_disk_rotation", "DiskRotationCheck"]


@dataclass
class NoiseReport:
    """Result of :func:`check_noise`.

    Attributes
    ----------
    v : ndarray
        Velocity of every channel of the cube (sorted).
    rms_ratio, peak_snr, frac_4sigma : ndarray
        Per channel: rms / median rms, peak / median rms, fraction of pixels > 4 sigma
        (pure Gaussian noise: 3.2e-5); NaN for blank channels.
    noise_v : ndarray
        Velocities of the channels used for the noise.
    cuts : dict
        Per cut ``y0 -> (current, blue_only, red_only)`` PV rms.
    flagged : ndarray
        Velocities of channels that look like they contain emission.
    """
    v: np.ndarray
    rms: np.ndarray
    rms_ratio: np.ndarray
    peak_snr: np.ndarray
    frac_4sigma: np.ndarray
    noise_v: np.ndarray
    cuts: dict
    vsys: float
    wind_window: tuple = None
    n_edge: int = 5
    units: str = ""
    flagged: np.ndarray = field(default=None)
    blank_v: np.ndarray = field(default=None)

    @property
    def red_blue_ratio(self):
        """Per cut, red-end over blue-end PV rms."""
        return {y: c[2] / c[1] for y, c in self.cuts.items()}

    @property
    def ok(self):
        """True if no noise channel is flagged and red/blue is within 15% at every cut."""
        used = np.isin(np.round(self.v, 4), np.round(self.noise_v, 4))
        bad = (self.frac_4sigma > 1e-3) | (self.peak_snr > 6)
        rb = np.array(list(self.red_blue_ratio.values()))
        return bool(not np.any(bad & used) and np.all(np.abs(rb - 1) < 0.15))

    def summary(self):
        L = [f"cube: {len(self.v)} channels, {self.v[0]:.2f} to {self.v[-1]:.2f} km/s; "
             f"noise channels at {np.array2string(np.sort(self.noise_v), precision=2)}"]
        if self.wind_window is not None:
            L.append(f"wind window {self.wind_window[0]:.2f}-{self.wind_window[1]:.2f} km/s; "
                     f"red end of the cube {self.v[-1] - self.wind_window[1]:+.2f} km/s beyond")
        if self.blank_v is not None and len(self.blank_v):
            L.append(f"blank channels (all zero/NaN, never used): "
                     f"{np.array2string(self.blank_v, precision=2)}")
        if self.flagged is not None and len(self.flagged):
            L.append(f"channels with emission?: {np.array2string(self.flagged, precision=2)}")
        L.append(f"{'y0':>6s} {'current':>10s} {'blue':>10s} {'red':>10s} {'red/blue':>9s}")
        for y, (a, b, c) in self.cuts.items():
            L.append(f"{y:6.2f} {a:10.4g} {b:10.4g} {c:10.4g} {c / b:9.3f}")
        L.append("OK: noise channels look line-free" if self.ok else
                 "CHECK: red/blue differs by >15% or a noise channel is flagged; consider "
                 "noise.vel_ranges_kms")
        return "\n".join(L)

    def plot(self, **kw):
        from .plotting import plot_noise_check
        return plot_noise_check(self, **kw)


def check_noise(cube, n_edge=None):
    """Check the noise channels of a :class:`~diskwind.cube.WindCube`.

    Parameters
    ----------
    n_edge : int, optional
        Channels per end for the blue/red comparison (default ``noise.n_edge_channels``).

    Returns
    -------
    NoiseReport
    """
    cfg = cube.config
    order = np.argsort(cube.full_vel)
    vs = cube.full_vel[order]
    nedge = n_edge or int(cfg.noise.n_edge_channels)
    allc = cube.regrid(order)
    blank = np.isin(order, cube.blank_channels)
    rms_ch = np.array([robust_std(c) for c in allc])
    rms_ch[blank] = np.nan
    ref = np.nanmedian(rms_ch)
    usable = np.where(~blank)[0]          # sorted-velocity positions with data
    with np.errstate(all="ignore"):
        peak = np.array([np.nanmax(c) for c in allc]) / ref
        f4 = np.array([np.nanmean(c > 4 * ref) for c in allc])
    cuts = {}
    for y0 in (cfg.y_offsets if cfg.y_offsets is not None else [1.0]):
        w = cfg.for_cut(y0).cuts.slit_width_arcsec or cube.slit_width
        sel = cube.rows(y0, w)
        r = []
        for arr in (cube.noise_cube, allc[usable[:nedge]], allc[usable[-nedge:]]):
            with np.errstate(all="ignore"):
                r.append(robust_std(np.nanmean(arr[:, sel, :], axis=1)))
        cuts[float(y0)] = tuple(r)
    vr = cfg.mask.vrel_range_kms
    peak[blank] = np.nan
    f4[blank] = np.nan
    return NoiseReport(v=vs, rms=rms_ch, rms_ratio=rms_ch / ref, peak_snr=peak,
                       frac_4sigma=f4, noise_v=np.asarray(cube.v_noise), cuts=cuts,
                       vsys=cube.vsys, n_edge=nedge, units=cube.units,
                       wind_window=None if vr is None else (vr[0] + cube.vsys, vr[1] + cube.vsys),
                       flagged=vs[(f4 > 1e-3) | (peak > 6)], blank_v=vs[blank])


@dataclass
class DiskRotationCheck:
    """Result of :func:`check_disk_rotation`.

    Attributes
    ----------
    scores : dict
        Summed emission (S/N units) along the projected Keplerian curves of each
        sense (+1, -1), at |v - Vsys| > 1.5 km/s.
    best_sense : int
        The sense the data favor.
    config_sense : int or None
        The value in the configuration (or estimated).
    pvs, curves : list
        PV diagrams of the cuts and, per cut, ``{(sense, h): v(x)}``.
    """
    y_offsets: list
    scores: dict
    best_sense: int
    config_sense: int
    pvs: list
    curves: list
    x: np.ndarray
    vsys: float

    @property
    def agrees(self):
        return self.config_sense is None or self.config_sense == self.best_sense

    def summary(self):
        return (f"emission along Keplerian curves (S/N units): +1: {self.scores[1]:.0f}, "
                f"-1: {self.scores[-1]:.0f} -> data favor disk_rotation_sense = "
                f"{self.best_sense:+d}; configuration: {self.config_sense}"
                + ("" if self.agrees else "  *** DISAGREES: check the figure"))

    def plot(self, **kw):
        from .plotting import plot_disk_rotation_check
        return plot_disk_rotation_check(self, **kw)


def check_disk_rotation(cube, y_offsets=(0.0, -0.15, -0.3)):
    """Compare the disk emission with Keplerian curves of both rotation senses.

    Cuts parallel to the major axis at ``y_offsets`` (<= 0: the side without wind)
    are compared with the projected Keplerian velocity of gas in the midplane and
    on the emitting surfaces ``z = +-H(r)``.

    Returns
    -------
    DiskRotationCheck
    """
    d = cube.dist
    xs = cube.x
    scores = {+1: 0.0, -1: 0.0}
    pvs, curves = [], []
    for y0 in y_offsets:
        pv = cube.pv(y0, with_mask=False)
        pvs.append(pv)
        cc = {}
        for sense in (+1, -1):
            for hs in (0.0, +1.0, -1.0):
                X = xs * d
                Y = y0 * d
                yy = np.full_like(X, cube.s * Y / max(cube.cosi, 1e-3))
                for _ in range(30):
                    r = np.hypot(X, yy)
                    h = hs * cube.disk_height(r)
                    yy = cube.s * (Y - h * cube.sini) / max(cube.cosi, 1e-3)
                r = np.hypot(X, yy)
                v = cube.vsys + sense * keplerian_speed(r, cube.mstar, h) * cube.sini * X \
                    / np.maximum(r, 1e-3)
                v = np.where((r <= cube.config.mask.disk_rout_au) & (r > 3), v, np.nan)
                cc[(sense, hs)] = v
                ok = np.isfinite(v) & (np.abs(v - cube.vsys) > 1.5)
                if ok.any():
                    vi = np.interp(v[ok], cube.v, np.arange(len(cube.v)))
                    xi = np.arange(len(xs))[ok].astype(float)
                    vals = map_coordinates(np.nan_to_num(pv.data), [vi, xi], order=1)
                    scores[sense] += np.sum(np.clip(vals, 0, None)) / max(pv.rms, 1e-30)
        curves.append(cc)
    best = +1 if scores[+1] > scores[-1] else -1
    return DiskRotationCheck(list(y_offsets), scores, best,
                             cube.config.source.disk_rotation_sense, pvs, curves, xs, cube.vsys)
