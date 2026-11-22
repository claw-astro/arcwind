"""
Trace extraction: the locus of the wind emission in each PV diagram.

At every position along a cut, a Gaussian is fitted to the spectrum around a
peak (Louvet et al. 2018, Sect. 4.3). Two modes:

* ``mode = "wings"`` (full rings): the blue and the red wing are fitted
  separately, split at the mean velocity (or a supplied mid-line);
* ``mode = "single"`` (partial rings): one peak per position.
* ``mode = "ring"`` (full rings with steep sides, e.g. tall, narrow rings): the
  "wings" points (one velocity per position, where the ring is flat) plus the
  walls traced channel by channel (one position per velocity, left and right of
  the ring), where the ring is steep and a column cuts it lengthwise.

In **guided** extraction a model curve (e.g. a previous global fit) restricts
each position to channels within ``guide_dv_kms`` of the model, so emission off
the wind's curve is ignored; ``peak_choice = "nearest"`` then traces the
significant peak closest to the model rather than the brightest.

Examples
--------
>>> from arcwind import TraceExtractor
>>> ext = TraceExtractor(cube)                                   # doctest: +SKIP
>>> tr = ext.extract(1.0)                    # unguided            # doctest: +SKIP
>>> tr.n, tr.x[:3], tr.v[:3], tr.err[:3]                         # doctest: +SKIP
>>> tr2 = ext.extract(1.0, guide=result.curve(1.0))   # guided by a fit   # doctest: +SKIP
>>> traces = ext.extract_all([0.5, 0.65, 0.8])        # dict y0 -> Trace  # doctest: +SKIP
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import curve_fit

from .constants import FWHM2SIG

__all__ = ["Trace", "TraceExtractor"]


@dataclass
class Trace:
    """Trace points of one cut.

    Attributes
    ----------
    y0 : float
        Offset of the cut (arcsec).
    x, v, err : ndarray
        Position (arcsec), fitted velocity (km/s) and velocity error (km/s).
    side : ndarray of int
        +1 red wing, -1 blue wing (points traced along velocity at a fixed position);
        +2 right wall, -2 left wall (points traced along position in a fixed channel,
        ``mode = "ring"``).
    snr, fwhm : ndarray
        Peak S/N and fitted FWHM (km/s for wing points, arcsec for wall points).
    xerr : ndarray
        Position error (arcsec) of wall points; NaN for wing points, whose position
        error is the fit's ``sig_x`` (see :meth:`x_errors`).
    """
    y0: float
    x: np.ndarray = field(default_factory=lambda: np.zeros(0))
    v: np.ndarray = field(default_factory=lambda: np.zeros(0))
    err: np.ndarray = field(default_factory=lambda: np.zeros(0))
    side: np.ndarray = field(default_factory=lambda: np.zeros(0, int))
    snr: np.ndarray = field(default_factory=lambda: np.zeros(0))
    fwhm: np.ndarray = field(default_factory=lambda: np.zeros(0))
    xerr: np.ndarray = field(default_factory=lambda: np.zeros(0))

    FIELDS = ("x", "v", "err", "side", "snr", "fwhm", "xerr")

    def __post_init__(self):
        self.xerr = np.asarray(self.xerr, float)
        if self.xerr.shape != np.shape(self.x):
            self.xerr = np.full(np.shape(self.x), np.nan)

    def x_errors(self, sig_x):
        """Position error of every point: ``xerr`` where set (wall points), else
        ``sig_x``."""
        return np.where(np.isfinite(self.xerr), self.xerr, sig_x)

    @property
    def n(self):
        """Number of trace points."""
        return len(self.x)

    def __len__(self):
        return self.n

    def thinned(self, dx, x0=0.0, step=None):
        """Keep about one point every ``dx`` arcsec along the cut (points on a
        regular grid of spacing ``step`` starting at ``x0``)."""
        if not dx or self.n == 0:
            return self
        step = step or float(np.min(np.diff(np.unique(self.x)))) if self.n > 1 else dx
        stride = max(1, int(round(dx / step)))
        if stride == 1:
            return self
        idx = np.rint((self.x - x0) / step).astype(int)
        # wall points (|side| = 2) are sampled in velocity, already one per smoothing length
        return self.subset((idx % stride == 0) | (np.abs(self.side) == 2))

    def __getitem__(self, key):            # dict-style access, e.g. tr["x"]
        return getattr(self, key)

    def subset(self, mask):
        """A new Trace with the points selected by ``mask`` (boolean or index array)."""
        return Trace(self.y0, *(getattr(self, k)[mask] for k in Trace.FIELDS))

    def as_dict(self):
        return {"y0": self.y0, **{k: getattr(self, k) for k in Trace.FIELDS}}

    @staticmethod
    def save_all(traces, path):
        """Save a dict ``{y0: Trace}`` to a ``.npz`` file."""
        d = {}
        for y0, t in traces.items():
            for k in Trace.FIELDS:
                d[f"{y0:.4f}__{k}"] = getattr(t, k)
        np.savez(path, **d)

    @staticmethod
    def load_all(path):
        """Read a ``.npz`` written by :meth:`save_all`; returns ``{y0: Trace}``."""
        z = np.load(path)
        ys = sorted({float(k.split("__")[0]) for k in z.files})
        return {y: Trace(y, *(z[f"{y:.4f}__{k}"] for k in Trace.FIELDS
                              if f"{y:.4f}__{k}" in z.files)) for y in ys}


def _gauss(v, a, mu, sig):
    return a * np.exp(-0.5 * ((v - mu) / sig) ** 2)


def _peak_fit(u, prof, allowed, rms, snr_min, fit_frac, drop_clipped, du, sig_max,
              near=None, floor_snr=0.0):
    """Gaussian fit of a 1-D profile ``prof(u)`` around its peak in ``allowed``.

    ``near``: if given, take the significant local maximum closest to any of these
    values of ``u`` instead of the brightest. The fitted window extends while the
    profile is above ``fit_frac`` x the peak and above ``floor_snr`` x the rms (at low
    S/N, ``fit_frac`` x the peak alone lets the window run into the noise).
    Returns ``(mu, sigma, snr)`` or None.
    """
    if allowed.sum() < 3:
        return None
    sp = np.where(allowed, prof, -np.inf)
    ip = int(np.argmax(sp))
    if near is not None:
        left = np.r_[-np.inf, sp[:-1]]
        right = np.r_[sp[1:], -np.inf]
        cand = np.where((sp >= left) & (sp >= right) & (sp > snr_min * rms))[0]
        if cand.size:
            dist = np.min(np.abs(u[cand][:, None] - np.asarray(near)[None, :]), axis=1)
            ip = int(cand[np.argmin(dist)])
    peak = prof[ip]
    if peak < snr_min * rms:
        return None
    n = len(u)
    edge = max(fit_frac * peak, floor_snr * rms)
    lo = hi = ip
    while lo - 1 >= 0 and allowed[lo - 1] and prof[lo - 1] > edge:
        lo -= 1
    while hi + 1 < n and allowed[hi + 1] and prof[hi + 1] > edge:
        hi += 1
    clipped = ((lo - 1 >= 0 and not allowed[lo - 1] and prof[lo - 1] > edge) or
               (hi + 1 < n and not allowed[hi + 1] and prof[hi + 1] > edge))
    if clipped and drop_clipped:
        return None
    if hi - lo < 2:
        if lo - 1 >= 0 and allowed[lo - 1]:
            lo -= 1
        if hi + 1 < n and allowed[hi + 1]:
            hi += 1
    if hi - lo < 2:
        return None
    uu, ss = u[lo:hi + 1], prof[lo:hi + 1]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p, _ = curve_fit(_gauss, uu, ss,
                             p0=[peak, u[ip], max(du, (hi - lo) * du / 2.355)],
                             bounds=([0, uu.min() - du, 0.3 * du],
                                     [10 * peak, uu.max() + du, sig_max]),
                             maxfev=2000)
    except (RuntimeError, ValueError):
        return None
    a, mu, sg = p
    if not (uu.min() - 0.5 * du <= mu <= uu.max() + 0.5 * du):
        return None
    return mu, sg, peak / rms


def _floor(t):
    """Noise floor (x rms) of the fitted window: ``trace.fit_floor_snr``; auto: 2 in
    ``mode = "ring"``, none otherwise (the behavior of earlier versions)."""
    if t.fit_floor_snr is not None:
        return float(t.fit_floor_snr)
    return 2.0 if t.mode == "ring" else 0.0


class TraceExtractor:
    """Extracts traces from the PV diagrams of a :class:`~arcwind.cube.WindCube`.

    All settings come from ``config.trace`` (with per-cut overrides); keyword
    arguments of :meth:`extract` override them for one call.
    """

    def __init__(self, cube):
        self.cube = cube

    def extract(self, y0, guide=None, v_split=None, pv=None, **options):
        """Trace of the cut at ``y0``.

        Parameters
        ----------
        y0 : float
            Offset of the cut (arcsec).
        guide : ModelCurve or dict, optional
            Model curve with ``X``, ``v`` and ``valid`` arrays (e.g.
            ``result.curve(y0)``): guided extraction.
        v_split : float or ndarray, optional
            Velocity separating the blue and red wings (``mode='wings'``); default:
            intensity-weighted mean velocity of the cut.
        pv : PVDiagram, optional
            Use this PV diagram (with its ``excluded`` mask) instead of computing it.
        **options
            Any ``[trace]`` option for this call only, e.g. ``snr_min=3``.

        Returns
        -------
        Trace
        """
        cube = self.cube
        cfg = cube.cut_config(y0)
        t = cfg.trace
        if options:
            from copy import copy
            t = copy(t)
            for k, val in options.items():
                if not hasattr(t, k):
                    raise ValueError(f"unknown trace option {k!r}")
                setattr(t, k, val)
        if pv is None:
            pv = cube.pv(y0)
        data, rms = pv.data, pv.rms
        excl = pv.excluded if pv.excluded is not None else cube.excluded_mask(y0)
        v, dv, xs = cube.v, cube.dv, cube.x
        snr_min = t.snr_min
        if v_split is None:
            ok = (data > 3 * rms) & ~excl
            w = np.where(ok, data, 0.0)
            v_split = (np.nansum(w * v[:, None]) / np.nansum(w)) if np.nansum(w) > 0 else cube.vsys
        v_split = np.broadcast_to(np.asarray(v_split, float), xs.shape)
        stride = 1
        if t.xsmooth_arcsec:
            stride = max(stride, max(1, int(round(t.xsmooth_arcsec / cube.step))))
        out = {k: [] for k in Trace.FIELDS}
        if guide is not None:
            gok = np.asarray(guide["valid"], bool)
            gX, gv = np.asarray(guide["X"])[gok], np.asarray(guide["v"])[gok]
            gdv = t.guide_dv_kms
            dlo, dhi = (gdv, gdv) if np.isscalar(gdv) else (gdv[0], gdv[1])
        nearest = guide is not None and t.peak_choice == "nearest"
        floor = _floor(t)
        for j in range(0, len(xs), stride):
            if excl[:, j].all():
                continue
            spec = data[:, j]
            if not np.any(np.isfinite(spec)):
                continue
            gmask, gvn = None, None
            if guide is not None:
                near = np.abs(gX - xs[j]) <= 0.5 * cube.beam_x + cube.step
                if not near.any():
                    continue
                gvn = gv[near]
                dvg = v[:, None] - gvn[None, :]
                gmask = np.any((dvg >= -dlo) & (dvg <= dhi), axis=1)
            for side in ((0,) if t.mode == "single" else (-1, +1)):
                allowed = (~excl[:, j]) & np.isfinite(spec)
                if gmask is not None:
                    allowed &= gmask
                if side < 0:
                    allowed &= v < v_split[j]
                elif side > 0:
                    allowed &= v > v_split[j]
                r = _peak_fit(v, spec, allowed, rms, snr_min, t.fit_frac, t.drop_clipped,
                              dv, 20.0, gvn if nearest else None, floor)
                if r is None:
                    continue
                mu, sg, snr = r
                out["x"].append(xs[j])
                out["v"].append(mu)
                out["err"].append(max(t.err_inflate * sg / snr, 0.05 * dv))
                out["side"].append(side if side != 0 else (1 if mu >= v_split[j] else -1))
                out["snr"].append(snr)
                out["fwhm"].append(sg / FWHM2SIG)
                out["xerr"].append(np.nan)
        if t.mode == "ring":
            self._walls(out, data, rms, excl, t, guide, v_split)
        return Trace(float(y0), *(np.array(out[k], int if k == "side" else float)
                                  for k in Trace.FIELDS))

    def _walls(self, out, data, rms, excl, t, guide, v_split):
        """``mode = "ring"``: trace the left and right walls of the ring channel by
        channel (one position per velocity), keep the wall points where the ring is
        steep (a narrow profile along X), and drop the column points whose velocity
        profile runs along such a wall (they are broad or clipped there)."""
        cube = self.cube
        v, dv, xs, step = cube.v, cube.dv, cube.x, cube.step
        bx = cube.beam_x
        sbeam = bx * FWHM2SIG
        nv = max(1, int(t.vsmooth_chan or 1))     # independent channels after smoothing
        dv_eff = nv * dv
        wmax = t.wall_max_fwhm_beam * bx
        floor = _floor(t)
        if guide is not None:
            gok = np.asarray(guide["valid"], bool)
            gX, gv = np.asarray(guide["X"])[gok], np.asarray(guide["v"])[gok]
            gdx = t.guide_dx_beam * bx
        ok = (data > 3 * rms) & ~excl
        w = np.where(ok, data, 0.0)
        wall = {k: [] for k in Trace.FIELDS}
        for k in range(0, len(v), nv):
            row = data[k, :]
            base = (~excl[k, :]) & np.isfinite(row)
            if base.sum() < 3:
                continue
            gxn = None
            if guide is not None:
                near = np.abs(gv - v[k]) <= 0.5 * dv_eff + dv
                if not near.any():
                    continue
                gxn = gX[near]
                base &= np.any(np.abs(xs[:, None] - gxn[None, :]) <= gdx, axis=1)
                if base.sum() < 3:
                    continue
            # left/right split: intensity-weighted centre of the channel
            # (or the middle of the model curve in this channel)
            if gxn is not None and gxn.size >= 2:
                xsplit = 0.5 * (gxn.min() + gxn.max())
            elif w[k].sum() > 0:
                xsplit = float(np.sum(w[k] * xs) / np.sum(w[k]))
            else:
                continue
            for side in (-2, +2):
                allowed = base & ((xs < xsplit) if side < 0 else (xs > xsplit))
                r = _peak_fit(xs, row, allowed, rms, t.snr_min, t.fit_frac, t.drop_clipped,
                              step, 10 * bx, gxn if (gxn is not None and
                                                    t.peak_choice == "nearest") else None,
                              floor)
                if r is None:
                    continue
                mu, sg, snr = r
                if sg / FWHM2SIG > wmax:          # flat part of the ring: columns trace it
                    continue
                wall["x"].append(mu)
                wall["v"].append(v[k])
                wall["err"].append(max(0.5 * dv_eff, 0.05 * dv))
                wall["xerr"].append(max(t.err_inflate * sg / snr, t.sigx_beam * sbeam))
                wall["side"].append(side)
                wall["snr"].append(snr)
                wall["fwhm"].append(sg / FWHM2SIG)
        if not wall["x"]:
            return
        # drop column points whose profile runs along a kept wall
        wx, wv = np.array(wall["x"]), np.array(wall["v"])
        keep = []
        for i in range(len(out["x"])):
            half = 0.5 * out["fwhm"][i]
            hit = (np.abs(wx - out["x"][i]) <= 0.5 * bx) & \
                  (np.abs(wv - out["v"][i]) <= max(half, dv_eff))
            keep.append(not hit.any())
        for key in Trace.FIELDS:
            out[key] = [val for val, kp in zip(out[key], keep) if kp] + wall[key]

    def extract_all(self, y_offsets=None, guide=None, min_points=None, **options):
        """Traces of several cuts.

        Parameters
        ----------
        y_offsets : array-like, optional
            Default: ``config.cuts.y_offsets_arcsec``.
        guide : FitResult, ConicalWindModel or callable, optional
            Anything with a ``curve(y0)`` method (or a function ``y0 -> curve``).
        min_points : int, optional
            Drop cuts with fewer points (default ``trace.min_points``; 0 keeps all).

        Returns
        -------
        dict
            ``{y0: Trace}`` in increasing ``y0``.
        """
        cube = self.cube
        ys = cube.config.y_offsets if y_offsets is None else np.atleast_1d(y_offsets)
        if ys is None:
            raise ValueError("no cuts: set cuts.y_offsets_arcsec or pass y_offsets")
        out = {}
        for y0 in ys:
            g = None
            if guide is not None:
                g = guide(y0) if callable(guide) and not hasattr(guide, "curve") \
                    else guide.curve(y0)
            tr = self.extract(float(y0), guide=g, **options)
            tc = cube.cut_config(y0).trace
            nmin = tc.min_points if min_points is None else min_points
            if tr.n >= nmin:              # counted before thinning (every pixel)
                dx = 0.5 * cube.beam_x if tc.dx_arcsec is None else tc.dx_arcsec
                out[float(y0)] = tr.thinned(dx, cube.x[0], cube.step)
            else:
                cube._log(f"cut y0={y0:.2f}\": only {tr.n} trace points, skipped")
        return out
