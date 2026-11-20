"""
Rotation measured cut by cut, with the rest of a global fit held fixed.

The global fit describes the specific angular momentum with one law for all
heights (a constant j, or j proportional to a power of z). When the rotation
signal is confined to some heights, such a law can extrapolate badly. Here the
geometry and the velocity law of a global fit are held fixed and a single j is
fitted to each cut on its own, which measures j(z) height by height without
assuming a law (the analogue of Louvet et al. 2018, Fig. 14).

Examples
--------
>>> rot = wa.rotation_per_cut(bpl)            # doctest: +SKIP
>>> print(rot.summary())                      # doctest: +SKIP
>>> wa.plot_rotation(rot)                     # doctest: +SKIP
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .constants import GM_SUN_AU_KMS2
from .model import ConicalWindModel

__all__ = ["RotationProfile", "rotation_per_cut"]


@dataclass
class RotationProfile:
    """Specific angular momentum measured in each cut (geometry and velocities fixed).

    Attributes
    ----------
    label : str
        Label of the global fit it is based on.
    y0, z, R : ndarray
        Cut offsets (arcsec), central heights (au) and cone radii there (au).
    j, j_err : ndarray
        Specific angular momentum fitted to each cut alone (au km/s) and its
        covariance error (scaled by sqrt(chi2_red) when > 1).
    n, d2 : ndarray
        Trace points and mean squared normalized distance of each cut.
    j_model : ndarray
        j of the global fit at the same heights.
    mstar_msun, sini : float
    nsigma : float
        Detection threshold used by :attr:`detected`.
    max_d2 : float
        Cuts whose mean squared distance exceeds this are flagged as poorly fitted
        (:attr:`poor`) and left out of the detections and joint fits.
    joint : dict
        One constant j fitted jointly to the detected cuts (``'detected'``) and to
        the others (``'undetected'``): ``{'j', 'err', 'y0'}``.
    """
    label: str
    y0: np.ndarray
    z: np.ndarray
    R: np.ndarray
    j: np.ndarray
    j_err: np.ndarray
    n: np.ndarray
    d2: np.ndarray
    j_model: np.ndarray
    mstar_msun: Optional[float]
    sini: float
    nsigma: float = 2.0
    max_d2: float = 2.0
    joint: dict = field(default_factory=dict)
    seed: int = 42

    @property
    def poor(self):
        """Cuts the model fits poorly even with their own j (mean d^2 > max_d2)."""
        return self.d2 > self.max_d2

    @property
    def detected(self):
        """Well-fitted cuts where |j| >= nsigma x its error."""
        return (np.abs(self.j) >= self.nsigma * self.j_err) & ~self.poor

    def median(self, quantity="j", n_boot=2000):
        """Median of ``quantity`` ('j' or 'Vphi') over the cuts that are not poorly
        fitted, with a 16-84% range from a bootstrap over those cuts.

        Returns
        -------
        dict
            ``{'value', 'lo', 'hi', 'n'}`` (``lo``/``hi``: the 16th/84th percentiles
            of the bootstrap medians).
        """
        v = (self.j if quantity == "j" else self.Vphi)[~self.poor]
        if len(v) == 0:
            return {"value": np.nan, "lo": np.nan, "hi": np.nan, "n": 0}
        rng = np.random.default_rng(self.seed)
        b = np.median(rng.choice(v, size=(n_boot, len(v)), replace=True), axis=1)
        lo, hi = np.percentile(b, [16, 84])
        return {"value": float(np.median(v)), "lo": float(lo), "hi": float(hi), "n": int(len(v))}

    @property
    def contiguous(self):
        """True if the detected cuts are adjacent (no undetected cut between them,
        poorly fitted cuts ignored)."""
        idx = np.where(self.detected)[0]
        if len(idx) < 2:
            return True
        keep = np.where(~self.poor)[0]
        pos = np.searchsorted(keep, idx)
        return bool(np.all(np.diff(pos) == 1))

    @property
    def Vphi(self):
        """Rotation velocity j/R (km/s)."""
        return self.j / self.R

    @property
    def Vphi_err(self):
        return self.j_err / self.R

    @property
    def VK(self):
        """Keplerian velocity at the cone radius (km/s), or NaN without M*."""
        if not self.mstar_msun:
            return np.full_like(self.R, np.nan)
        return np.sqrt(GM_SUN_AU_KMS2 * self.mstar_msun / self.R)

    @property
    def jK(self):
        """Keplerian specific angular momentum at the cone radius (au km/s)."""
        if not self.mstar_msun:
            return np.full_like(self.R, np.nan)
        return np.sqrt(GM_SUN_AU_KMS2 * self.mstar_msun * self.R)

    @property
    def r_L(self):
        """Launch radius without magnetic torque, j^2/GM* (au)."""
        if not self.mstar_msun:
            return np.full_like(self.R, np.nan)
        return self.j ** 2 / (GM_SUN_AU_KMS2 * self.mstar_msun)

    def table(self):
        """One line per cut."""
        L = ["# y0(\")  z(au)  R(au)  n   j(au km/s)   err   Vphi(km/s)  Vphi_sini  |Vphi|/VK  "
             "r_L(au)  j_model  <d2>  detected"]
        for k in range(len(self.y0)):
            L.append(f"{self.y0[k]:6.2f} {self.z[k]:6.0f} {self.R[k]:6.0f} {self.n[k]:3d} "
                     f"{self.j[k]:10.1f} {self.j_err[k]:6.1f} {self.Vphi[k]:10.2f} "
                     f"{self.Vphi[k] * self.sini:9.2f} {abs(self.Vphi[k]) / self.VK[k]:9.2f} "
                     f"{self.r_L[k]:8.1f} {self.j_model[k]:8.1f} {self.d2[k]:5.2f}  "
                     f"{'poor fit' if self.poor[k] else 'yes' if self.detected[k] else 'no'}")
        return "\n".join(L)

    def summary(self):
        """Table, joint values and an upper limit above the detections."""
        L = [f"# rotation per cut, geometry and velocity law of '{self.label}' held fixed",
             f"# a single j fitted to each cut alone; detected: |j| >= {self.nsigma:g} sigma; "
             f"poor fit: <d2> > {self.max_d2:g} (left out)",
             self.table(), ""]
        det = self.detected
        if det.any():
            zr = self.z[det]
            L.append(f"rotation detected in {det.sum()} of {len(det)} cuts, z = "
                     f"{zr.min():.0f}-{zr.max():.0f} au")
        else:
            L.append("rotation not detected in any single cut")
        m = self.median()
        L.append(f"median j of the {m['n']} well-fitted cuts: {m['value']:.1f} au km/s "
                 f"(16-84% bootstrap range {m['lo']:.1f} to {m['hi']:.1f})")
        for key, txt in (("detected", "detected cuts"), ("undetected", "other cuts")):
            d = self.joint.get(key)
            if d:
                L.append(f"one j fitted to the {txt} ({', '.join(f'{y:.2f}' for y in d['y0'])}): "
                         f"{d['j']:.1f} +/- {d['err']:.1f} au km/s")
        d = self.joint.get("undetected")
        if d:
            lim = abs(d["j"]) + self.nsigma * d["err"]
            L.append(f"  -> |j| < {lim:.0f} au km/s ({self.nsigma:g} sigma) where rotation is "
                     "not detected")
        sense = None
        if det.any():
            s = np.sign(self.j[det])
            sense = "all negative" if (s < 0).all() else "all positive" if (s > 0).all() \
                else "mixed signs"
            L.append(f"sign of the detections: {sense} (j < 0: against a disk of sense +1)")
            if not self.contiguous:
                L.append("the detected cuts are NOT adjacent: cuts between them show no "
                         "rotation, so these may be isolated excursions rather than a "
                         "rotating zone")
        L.append("NB: the axis (x0, xi) is held fixed, so j of a cut also absorbs any error of "
                 "its ring-center position.")
        return "\n".join(L)

    def save(self, path):
        with open(path, "w") as fh:
            fh.write(self.summary() + "\n")


def _fit_constant_j(result, traces, j0, loss, equal):
    """Best constant j (and error) for some cuts, everything else fixed."""
    from .fitting import GlobalFitter
    gm = result.model
    m = ConicalWindModel(gm.law, radial=gm.radial, rotation="constant", zref_au=gm.zref_au)
    m.fix(**{k: v for k, v in result.params.items() if k in m.names and k != "j"})
    m.set(j=float(j0))
    f = GlobalFitter(m, result.geometry, traces, result.sig_x, loss=loss,
                     equal_cut_weight=equal, j_multistart=True, verbose=False, n_jobs=1)
    r = f.fit(init={"j": float(j0)})
    return float(r.params["j"]), float(r.errors("max")["j"]), r


def rotation_per_cut(result, traces=None, nsigma=2.0, max_d2=2.0, min_points=4):
    """Measure j in each cut with the geometry and velocity law of ``result`` fixed.

    Parameters
    ----------
    result : FitResult
        A global fit (any velocity law and rotation law).
    traces : dict, optional
        ``{y0: Trace}``; default: the traces of ``result``.
    nsigma : float
        Detection threshold (|j| / error).
    max_d2 : float
        Cuts with a larger mean squared distance (with their own j) are flagged as
        poorly fitted and left out of the detections and joint fits.
    min_points : int
        Cuts with fewer trace points are skipped.

    Returns
    -------
    RotationProfile
    """
    g = result.geometry
    traces = result.traces if traces is None else traces
    st = result.settings or {}
    loss = st.get("loss", "soft_l1")
    equal = st.get("equal_cut_weight", True)
    p = result.params
    rows = []
    for y0, tr in sorted(traces.items()):
        if tr.n < min_points:
            continue
        z = g.zc_au(y0)
        jm = float(result.model.j_of_z(z))
        try:
            j, e, r = _fit_constant_j(result, {y0: tr}, jm, loss, False)
        except Exception:                      # e.g. a cut the model does not cross
            continue
        R = p["R0_au"] + z * p["tan_theta"]
        rows.append((y0, z, R, j, e, tr.n, r.chi2 / max(tr.n, 1), jm))
    if not rows:
        raise ValueError("no cut could be fitted")
    a = np.array(rows, float)
    prof = RotationProfile(result.label, a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4],
                           a[:, 5].astype(int), a[:, 6], a[:, 7], g.mstar_msun,
                           float(np.sin(np.radians(g.inclination_deg))), float(nsigma),
                           float(max_d2))
    det = prof.detected
    for key, sel in (("detected", det), ("undetected", ~det & ~prof.poor)):
        ys = prof.y0[sel]
        if len(ys):
            tsub = {float(y): traces[float(y)] for y in ys}
            j0 = float(np.average(prof.j[sel], weights=prof.j_err[sel] ** -2.0))
            try:
                j, e, _ = _fit_constant_j(result, tsub, j0, loss, equal)
                prof.joint[key] = {"j": j, "err": e, "y0": [float(y) for y in ys]}
            except Exception:
                pass
    return prof
