"""
Independent ring fits to each cut (Louvet et al. 2018, Sect. 4.4 and App. B),
generalized to any inclination.

Each cut is fitted on its own with a locally conical ring
``R(z) = Rc + (z - zc) tan(theta)`` of center ``x_off`` and velocities
``vcent`` (ring center), ``A = Vr sin i`` and ``B = Vphi sin i``. This is the
method of the paper for **full rings** (both halves of every ring detected),
e.g. edge-on sources like HH 30; with partial rings the per-cut fit is
degenerate and the global fit (:mod:`arcwind.fitting`) should be used.

Examples
--------
>>> from arcwind import RingFitter
>>> rf = RingFitter(cube)                                  # doctest: +SKIP
>>> fit = rf.fit_cut(0.75)              # one cut          # doctest: +SKIP
>>> fit.R_au, fit.errors["R_au"], fit.Vphi, fit.j           # doctest: +SKIP
>>> rings = rf.fit_all()                # all cuts, two-pass cone   # doctest: +SKIP
>>> rings.cone                          # R0, theta from R(z)       # doctest: +SKIP
>>> rings.table()                       # one row per cut            # doctest: +SKIP
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import curve_fit, least_squares

from .constants import FWHM2SIG
from .model import ModelCurve
from .trace import TraceExtractor

__all__ = ["RingFitter", "RingFit", "RingFitSet", "ring_curve", "RING_PARAMS"]

RING_PARAMS = ["Rc", "xoff", "vcent", "A", "B", "T"]
DERIVED_KEYS = ["z_au", "R_au", "xoff_au", "vcent_rel", "Vr_sini", "Vphi_sini", "tan_theta",
                "Vz", "Vr", "Vphi", "Vp", "theta_v_deg", "j", "j_sini"]


def ring_curve(p, y0, geom, nphi=721, rot_law="j"):
    """Locus in (X, v) of a locally conical ring cut by the cut at ``y0``.

    Parameters
    ----------
    p : sequence
        ``(Rc [arcsec], xoff [arcsec], vcent [km/s], A, B, T)``.
    y0 : float
    geom : ViewGeometry
    rot_law : {'j', 'const'}
        ``'j'``: V_phi ~ 1/R across the cut; ``'const'``: constant (paper App. B).

    Returns
    -------
    ModelCurve
    """
    Rc, xoff, vcent, A, B, T = p
    si, ci, s = geom.sini, geom.cosi, geom.s
    phi = np.linspace(0.0, 2 * np.pi, nphi)
    sp, cp = np.sin(phi), np.cos(phi)
    zc = y0 / si
    denom = si + s * T * sp * ci
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (y0 - s * sp * ci * (Rc - zc * T)) / denom
        R = Rc + (z - zc) * T
        X = xoff + R * cp
        v = vcent + A * sp + (B * (Rc / R) if rot_law == "j" else B) * cp
    valid = (denom > 0) & (z > 0) & (R > 0)
    return ModelCurve(phi=phi, X=X, v=v, z=z * geom.distance_pc, R=R * geom.distance_pc,
                      valid=valid)


@dataclass
class RingFit:
    """Fit of one cut.

    Attributes
    ----------
    y0 : float
    ok : bool
    params : ndarray
        ``Rc, xoff, vcent, A, B, T`` (arcsec, km/s).
    values, errors : dict
        Derived quantities (``z_au, R_au, xoff_au, vcent_rel, Vr_sini, Vphi_sini,
        tan_theta, Vz, Vr, Vphi, Vp, theta_v_deg, j, j_sini``) and their 1-sigma
        Monte Carlo errors.
    n, chi2r, coverage : number of points, reduced chi2, azimuthal coverage.
    good : bool
        Coverage >= ``ring.min_phi_coverage``.
    """
    y0: float
    ok: bool
    reason: str = ""
    params: np.ndarray = None
    free: np.ndarray = None
    cov: np.ndarray = None
    draws: np.ndarray = None
    values: dict = field(default_factory=dict)
    errors: dict = field(default_factory=dict)
    n: int = 0
    chi2r: float = np.nan
    coverage: float = 0.0
    good: bool = False
    flat: bool = False
    z_range_au: tuple = (np.nan, np.nan)
    trace: object = None

    def __getattr__(self, k):
        v = self.__dict__.get("values", {})
        if k in v:
            return v[k]
        raise AttributeError(k)

    def curve(self, geom, rot_law="j"):
        return ring_curve(self.params, self.y0, geom, rot_law=rot_law)


@dataclass
class RingFitSet:
    """Ring fits of all cuts, plus the cone fitted to R(z)."""
    fits: list
    cone: dict = None

    def good(self):
        return [f for f in self.fits if f.ok and f.good]

    def table(self):
        """Dict of columns (``y0``, ``n``, ``chi2r``, ``coverage``, ``good`` and every
        derived quantity with ``_err``)."""
        rows = [f for f in self.fits if f.ok]
        out = dict(y0=np.array([f.y0 for f in rows]), n=np.array([f.n for f in rows]),
                   chi2r=np.array([f.chi2r for f in rows]),
                   coverage=np.array([f.coverage for f in rows]),
                   good=np.array([f.good for f in rows]))
        for k in DERIVED_KEYS:
            out[k] = np.array([f.values[k] for f in rows])
            out[k + "_err"] = np.array([f.errors[k] for f in rows])
        return out

    def save_csv(self, path):
        t = self.table()
        cols = list(t)
        with open(path, "w") as fh:
            fh.write(",".join(cols) + "\n")
            for i in range(len(t["y0"])):
                fh.write(",".join(f"{float(t[c][i]):.6g}" for c in cols) + "\n")
        return path


class RingFitter:
    """Per-cut ring fits on a :class:`~arcwind.cube.WindCube` (settings: ``[ring]``)."""

    def __init__(self, cube, rng=None):
        self.cube = cube
        self.cfg = cube.config
        self.geom = cube.geometry
        self.rng = rng or np.random.default_rng(self.cfg.fit.seed)
        self.extractor = TraceExtractor(cube)

    @property
    def sig_x(self):
        return self.cfg.trace.sigx_beam * self.cube.beam_x * FWHM2SIG

    def _nearest(self, p, y0, tr, nphi=361):
        c = ring_curve(p, y0, self.geom_arcsec, nphi=nphi, rot_law=self.cfg.ring.rotation_law)
        ok = c.valid
        if ok.sum() < 10:
            return None, None
        X, v, phi = c.X[ok], c.v[ok], c.phi[ok]
        d2 = ((tr.x[:, None] - X[None, :]) / self.sig_x) ** 2 + \
             ((tr.v[:, None] - v[None, :]) / tr.err[:, None]) ** 2
        j = np.argmin(d2, axis=1)
        return np.sqrt(d2[np.arange(len(j)), j]), phi[j]

    @property
    def geom_arcsec(self):
        """Geometry with distance 1 (ring parameters are in arcsec)."""
        from .geometry import ViewGeometry
        g = self.geom
        return ViewGeometry(g.inclination_deg, g.s, 1.0, g.vsys, g.vsys_err, g.disk_sense)

    def coverage(self, p, y0, tr, nbin=12):
        """Fraction of 30-degree azimuth bins of the ring containing trace points."""
        _, ph = self._nearest(p, y0, tr)
        if ph is None:
            return 0.0
        return len(np.unique((ph / (2 * np.pi) * nbin).astype(int) % nbin)) / nbin

    def derive(self, p, y0, vsys):
        """Physical quantities from ring parameters (works on arrays of draws)."""
        Rc, xoff, vcent, A, B, T = p
        g, d = self.geom, self.geom.distance_pc
        vrel = vcent - vsys
        out = dict(z_au=np.broadcast_to(y0 / g.sini * d, np.shape(Rc)) * 1.0, R_au=Rc * d,
                   xoff_au=xoff * d, vcent_rel=vrel, Vr_sini=A, Vphi_sini=B, tan_theta=T,
                   j_sini=Rc * d * B, Vr=A / g.sini, Vphi=B / g.sini, j=Rc * d * B / g.sini)
        if not g.edge_on:
            Vz = -g.s * vrel / g.cosi
            out["Vz"] = Vz
            out["Vp"] = np.sqrt(Vz ** 2 + (A / g.sini) ** 2)
            with np.errstate(divide="ignore", invalid="ignore"):
                out["theta_v_deg"] = np.degrees(np.arctan2(A / g.sini, Vz))
        else:
            for k in ("Vz", "Vp", "theta_v_deg"):
                out[k] = np.full(np.shape(Rc), np.nan)
        return out

    def fit_trace(self, y0, tr, T_fixed=None, fix=None, prior=None, nmc=None):
        """Fit the ring model to a trace.

        Parameters
        ----------
        y0 : float
        tr : Trace
        T_fixed : float, optional
            Fix tan(theta).
        fix : dict, optional
            Fix ring parameters, e.g. ``{'vcent': 4.0, 'A': 2.1}``.
        prior : dict, optional
            Gaussian priors ``{name: (mean, sigma)}`` on ring parameters.
        nmc : int, optional
            Monte Carlo draws (default ``ring.n_mc``).

        Returns
        -------
        RingFit
        """
        cfg, cube, g = self.cfg.for_cut(y0), self.cube, self.geom
        nmc = cfg.ring.n_mc if nmc is None else nmc
        n = tr.n
        single = cfg.trace.mode == "single"
        if n < cfg.trace.min_points or (not single and ((tr.side > 0).sum() < 2
                                                        or (tr.side < 0).sum() < 2)):
            return RingFit(y0, False, reason=f"too few trace points ({n})", trace=tr)
        if g.edge_on:
            T_fixed = 0.0 if T_fixed is None else T_fixed
        Tmax = 0.95 * g.sini / max(g.cosi, 1e-6)
        free = np.array([True] * 5 + [T_fixed is None])
        fix = fix or {}
        for k in fix:
            free[RING_PARAMS.index(k)] = False
        red, blue = tr.side > 0, tr.side < 0
        R0 = max(0.5 * (np.percentile(tr.x, 97) - np.percentile(tr.x, 3)), 2 * cube.step)
        xoff0 = 0.5 * (np.percentile(tr.x, 97) + np.percentile(tr.x, 3))
        vlo, vhi = np.percentile(tr.v, 5), np.percentile(tr.v, 95)
        starts_v = [(0.5 * (vlo + vhi), max(0.5 * (vhi - vlo), 0.2)),
                    (vlo, max(vhi - vlo, 0.2)), (vhi, max(vhi - vlo, 0.2))]
        if red.sum() >= 2 and blue.sum() >= 2 and not single:
            starts_v.insert(0, (0.5 * (np.median(tr.v[red]) + np.median(tr.v[blue])),
                                max(0.5 * (np.percentile(tr.v[red], 90)
                                           - np.percentile(tr.v[blue], 10)), 0.2)))
        vc0, A0 = starts_v[0]
        T0 = T_fixed if T_fixed is not None else min(np.tan(np.radians(30.0)), 0.8 * Tmax)
        lo = np.array([0.3 * cube.step, xoff0 - 3 * R0 - 1, vlo - 30, 0.0, -30.0, 0.0])
        hi = np.array([10 * R0 + 2.0, xoff0 + 3 * R0 + 1, vhi + 30, 30.0, 30.0,
                       Tmax if not g.edge_on else 10.0])

        def unpack(q):
            p = np.array([0, 0, 0, 0, 0, T0 if T_fixed is None else T_fixed], float)
            for k, val in fix.items():
                p[RING_PARAMS.index(k)] = val
            p[free] = q
            return p

        starts = [(vc0, A0, B0) for B0 in (0.0, 0.3 * A0, -0.3 * A0)] + \
                 [(vc, A, 0.0) for vc, A in starts_v[1:]]
        prior = prior or {}
        pri = [(RING_PARAMS.index(k), m, s) for k, (m, s) in prior.items()]
        if prior:
            vcp = prior.get("vcent", (vc0, 0))[0]
            Ap = max(prior.get("A", (A0, 0))[0], 0.2)
            starts = [(vcp, Ap, B0) for B0 in (0.0, 0.3 * Ap, -0.3 * Ap)] + starts
        if fix:
            vcf, Af = fix.get("vcent", vc0), fix.get("A", A0)
            starts = [(vcf, Af, B0) for B0 in (0.0, 0.3 * Af, -0.3 * Af)]

        def resid(q):
            p = unpack(q)
            r, _ = self._nearest(p, y0, tr)
            r = np.full(n, 1e3) if r is None else r
            if pri:
                r = np.concatenate([r, [(p[k] - m) / s for k, m, s in pri]])
            return r
        best = None
        for vc_s, A_s, B_s in starts:
            p0 = np.clip(np.array([R0, xoff0, vc_s, A_s, B_s, T0])[free],
                         lo[free] + 1e-6, hi[free] - 1e-6)
            try:
                res = least_squares(resid, p0, bounds=(lo[free], hi[free]),
                                    loss=self.cfg.fit.loss, x_scale="jac", max_nfev=4000)
            except ValueError:
                continue
            if best is None or res.cost < best.cost:
                best = res
        if best is None:
            return RingFit(y0, False, reason="optimizer failed", trace=tr)
        p = unpack(best.x)
        rr, _ = self._nearest(p, y0, tr)
        dof = max(n - free.sum(), 1)
        chi2r = float(np.sum(rr ** 2) / dof)
        cov = np.linalg.pinv(best.jac.T @ best.jac)
        if self.cfg.fit.scale_errors_by_chi2 and chi2r > 1:
            cov = cov * chi2r
        ev, evec = np.linalg.eigh(0.5 * (cov + cov.T))
        ev = np.clip(ev, 0, None)
        dq = best.x + (self.rng.standard_normal((nmc, len(best.x))) * np.sqrt(ev)) @ evec.T
        draws = np.tile(p, (nmc, 1))
        draws[:, free] = dq
        vs = g.vsys + g.vsys_err * self.rng.standard_normal(nmc)
        best_d = self.derive(p, y0, g.vsys)
        mc = self.derive(draws.T, y0, vs)
        cover = self.coverage(p, y0, tr)
        vspan = np.percentile(tr.v, 95) - np.percentile(tr.v, 5)
        c = ring_curve(p, y0, self.geom_arcsec)
        zv = c.z[c.valid] * g.distance_pc
        return RingFit(y0, True, params=p, free=free, cov=cov, draws=draws,
                       values={k: float(v) for k, v in best_d.items()},
                       errors=_mc_errors(mc, best_d),
                       n=n, chi2r=chi2r, coverage=cover,
                       good=bool(cover >= cfg.ring.min_phi_coverage),
                       flat=bool(p[3] < 0.35 * vspan and single),
                       z_range_au=(float(zv.min()), float(zv.max())) if zv.size else
                       (np.nan, np.nan), trace=tr)

    def fit_cut(self, y0, T_fixed=None, fix=None, prior=None, guide=None, nmc=None):
        """Extract the trace of the cut and fit it; for ``trace.mode='wings'`` the
        wings are re-split along the fitted mid-line and the fit repeated.

        Returns
        -------
        RingFit
        """
        pv = self.cube.pv(y0)
        tr = self.extractor.extract(y0, guide=guide, pv=pv)
        fit = self.fit_trace(y0, tr, T_fixed=T_fixed, fix=fix, prior=prior, nmc=nmc)
        if self.cfg.for_cut(y0).trace.mode != "single" and fit.ok:
            Rc, xoff, vcent, A, B, T = fit.params
            mid = vcent + B * np.clip((self.cube.x - xoff) / Rc, -1, 1)
            tr2 = self.extractor.extract(y0, guide=guide, v_split=mid, pv=pv)
            fit2 = self.fit_trace(y0, tr2, T_fixed=T_fixed, fix=fix, prior=prior, nmc=nmc)
            if fit2.ok:
                fit = fit2
        return fit

    def fit_all(self, y_offsets=None, prior_from=None, fix_from=None, guide=None):
        """Fit every cut.

        Without a fixed cone angle (``ring.cone_half_angle_deg``), the fits are done
        in ``ring.n_passes`` passes: fit all cuts, fit R(z) = R0 + z tan(theta) to the
        good cuts, refit with theta fixed.

        Parameters
        ----------
        y_offsets : array-like, optional
        prior_from, fix_from : FitResult, optional
            Gaussian priors on (or fixed values of) ``vcent`` and ``A`` and a fixed
            theta from a global fit, at each cut's height.
        guide : FitResult, optional
            Guided trace extraction.

        Returns
        -------
        RingFitSet
        """
        cfg, g = self.cfg, self.geom
        ys = cfg.y_offsets if y_offsets is None else np.atleast_1d(y_offsets)
        ca = cfg.ring.cone_half_angle_deg
        T_fixed = None if ca is None else np.tan(np.radians(ca))
        src = prior_from or fix_from
        npass = 1 if (T_fixed is not None or g.edge_on or src is not None) else cfg.ring.n_passes
        cone = None
        for ipass in range(npass):
            fits = []
            for y0 in ys:
                kw = dict(T_fixed=T_fixed)
                if src is not None:
                    ra = src.ring_at(y0)
                    kw["T_fixed"] = ra["T"]
                    if fix_from is not None:
                        kw["fix"] = dict(vcent=ra["vcent"], A=ra["A"])
                    else:
                        s = cfg.ring.prior_sigma_kms
                        kw["prior"] = dict(vcent=(ra["vcent"], s), A=(ra["A"], s))
                gcurve = guide.curve(y0) if guide is not None else None
                f = self.fit_cut(float(y0), guide=gcurve, **kw)
                fits.append(f)
                self.cube._log(_line(f))
            fs = RingFitSet(fits)
            cone = self.fit_cone(fs)
            if cone is None:
                break
            self.cube._log(f"cone: R0 = {cone['R0']:.1f} +/- {cone['R0_err']:.1f} au, theta = "
                           f"{cone['theta']:.1f} +/- {cone['theta_err']:.1f} deg")
            if T_fixed is None and ipass < npass - 1:
                T_fixed = float(np.clip(cone["tan"], 0, 0.95 * g.sini / max(g.cosi, 1e-6)))
        return RingFitSet(fits, cone)

    def fit_cone(self, fitset):
        """Weighted straight-line fit R(z) = R0 + z tan(theta) to the good cuts."""
        t = fitset.table()
        if not len(t["y0"]):
            return None
        z, R, eR = t["z_au"], t["R_au"], t["R_au_err"]
        lo, hi = self.cfg.ring.cone_fit_zrange_au
        m = (z >= lo) & (z <= hi) & np.isfinite(eR) & (eR > 0) & t["good"]
        if m.sum() < 3:
            return None
        p, c = curve_fit(lambda zz, r0, tt: r0 + tt * zz, z[m], R[m], p0=[R[m].min(), 0.5],
                         sigma=eR[m], absolute_sigma=True)
        e = np.sqrt(np.diag(c))
        return dict(R0=p[0], R0_err=e[0], tan=p[1], tan_err=e[1],
                    theta=np.degrees(np.arctan(p[1])),
                    theta_err=np.degrees(e[1] / (1 + p[1] ** 2)))


def _mc_errors(mc, keys):
    """Standard deviations of the Monte Carlo draws (NaN for undefined quantities,
    e.g. Vz of an edge-on ring)."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return {k: float(np.nanstd(mc[k])) if np.isfinite(mc[k]).any() else np.nan
                for k in keys}


def _line(f):
    if not f.ok:
        return f"y0={f.y0:5.2f}\"  no fit: {f.reason}"
    v, e = f.values, f.errors
    return (f"y0={f.y0:5.2f}\" z={v['z_au']:6.1f} au  R={v['R_au']:6.1f}+/-{e['R_au']:.1f}  "
            f"Vc-Vsys={v['vcent_rel']:+5.2f}  Vr sin i={v['Vr_sini']:5.2f}  "
            f"Vphi sin i={v['Vphi_sini']:+5.2f}  chi2r={f.chi2r:.2f} (n={f.n}) "
            f"cov={f.coverage:.2f}" + ("" if f.good else "  [partial ring]"))
