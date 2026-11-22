"""
Global fit of a conical wind to the traces of all cuts at once.

The model curve of every cut (:meth:`ConicalWindModel.curve`) is compared with
the trace points of that cut through their normalized orthogonal distance

.. math::

   d_k^2 = \\min_\\phi\\left[\\left(\\frac{x_k - X(\\phi)}{\\sigma_x}\\right)^2
           + \\left(\\frac{v_k - v(\\phi)}{\\sigma_k}\\right)^2\\right],

optionally weighted so that every cut counts equally
(``fit.equal_cut_weight``), and :math:`\\sum_k \\rho(w_k^2 d_k^2)` is minimized
with ``scipy.optimize.least_squares`` (trust-region reflective, loss
``fit.loss``) from several starting points. Uncertainties come from the
covariance matrix and from a bootstrap over cuts.

Examples
--------
>>> from arcwind import GlobalFitter, ConicalWindModel
>>> model = ConicalWindModel("power-law", radial="cone")
>>> fitter = GlobalFitter(model, cube.geometry, traces, sig_x=0.02)   # doctest: +SKIP
>>> res = fitter.fit(n_boot=0)                                        # doctest: +SKIP
>>> print(res.summary())                                              # doctest: +SKIP

Usually you call :meth:`arcwind.WindAnalysis.fit`, which builds the fitter
from the configuration.
"""
from __future__ import annotations

import json
import warnings

import numpy as np
from scipy.optimize import least_squares

from .geometry import ViewGeometry
from .model import ConicalWindModel
from .trace import Trace

__all__ = ["GlobalFitter", "FitResult", "odr_distances", "compare_fits"]


def odr_distances(curve, trace, sig_x):
    """Normalized orthogonal distance of every trace point to a model curve.

    Parameters
    ----------
    curve : ModelCurve
    trace : Trace
    sig_x : float
        Position uncertainty (arcsec) of the points without their own
        (``trace.xerr``; wall points of ``trace.mode = "ring"`` have one).

    Returns
    -------
    ndarray or None
        ``d_k`` for each point (None if the curve has < 10 valid points).
    """
    ok = curve.valid
    if ok.sum() < 10:
        return None
    X, v = curve.X[ok], curve.v[ok]
    sx = trace.x_errors(sig_x)
    d2 = ((trace.x[:, None] - X[None, :]) / sx[:, None]) ** 2 + \
         ((trace.v[:, None] - v[None, :]) / trace.err[:, None]) ** 2
    return np.sqrt(d2.min(axis=1))


def cut_weights(traces, equal=True):
    """Per-cut weights ``sqrt(<n>/n_c)`` (every cut counts equally), or None."""
    if not equal:
        return None
    n = np.array([max(t.n, 1) for t in traces], float)
    return np.sqrt(n.mean() / n)


class GlobalFitter:
    """Fits a :class:`ConicalWindModel` to traces of several cuts.

    Parameters
    ----------
    model : ConicalWindModel
        Model with its starting values, fixed parameters and bounds (not modified).
    geometry : ViewGeometry
    traces : dict
        ``{y0: Trace}``.
    sig_x : float
        Position uncertainty of trace points (arcsec), typically
        ``trace.sigx_beam`` x the beam dispersion along X.
    loss : str
        ``least_squares`` loss ('linear', 'soft_l1', ...).
    equal_cut_weight : bool
        Weight cuts equally.
    max_nfev, nphi : int
        Function evaluations per start; azimuths per model curve.
    scale_errors_by_chi2 : bool
        Scale covariance errors by chi2_red when > 1.
    check_vr_sign : bool or 'prefer_positive'
        With ``radial='proportional'``: if Vr < 0, refit with Vr >= 0.
    j_multistart : bool
        Restart with j = +-j_s, +-3 j_s after the first stage.
    seed : int
    verbose : bool
    counter_traces : dict, optional
        ``{y0: Trace}`` of the opposite lobe, extracted in its own frame
        (``WindCube(..., lobe="counter")``). Both lobes are then fitted with one
        mirror-symmetric wind: the counter lobe sees the model through
        :meth:`ConicalWindModel.mirrored` and :meth:`ViewGeometry.mirrored`.
    """

    def __init__(self, model, geometry, traces, sig_x, loss="linear", equal_cut_weight=False,
                 max_nfev=3000, nphi=181, scale_errors_by_chi2=True, check_vr_sign=True,
                 j_multistart=True, seed=42, verbose=True, counter_traces=None,
                 error_method="bootstrap", n_jobs=-1, resample_tol=1e-6,
                 screen_tol=1e-5, screen_nfev=600, n_refine=4):
        self.model = model.copy()
        self.geom = geometry
        self.traces = {float(k): v for k, v in sorted(traces.items())}
        self.counter_traces = {float(k): v for k, v in sorted((counter_traces or {}).items())}
        if len(self.traces) + len(self.counter_traces) == 0:
            raise ValueError("no traces to fit")
        # every cut of every lobe: (lobe index, y0, trace); lobe 1 = counter lobe
        self.groups = [(geometry, False), (geometry.mirrored(), True)]
        self.entries = ([(0, y, t) for y, t in self.traces.items()]
                        + [(1, y, t) for y, t in self.counter_traces.items()])
        self.msign = self.model.mirror_signs()
        self.sig_x = float(sig_x)
        self.loss = loss
        self.equal = bool(equal_cut_weight)
        self.max_nfev = int(max_nfev)
        self.nphi = int(nphi)
        self.scale = bool(scale_errors_by_chi2)
        self.check_vr_sign = check_vr_sign
        self.j_multistart = bool(j_multistart)
        self.rng = np.random.default_rng(seed)
        self.seed = seed
        self.verbose = verbose
        if error_method not in ("bootstrap", "jackknife"):
            raise ValueError("error_method must be 'bootstrap' or 'jackknife'")
        self.error_method = error_method
        self.n_jobs = n_jobs
        self.resample_tol = float(resample_tol)
        self.screen_tol = float(screen_tol)
        self.screen_nfev = int(screen_nfev)
        self.n_refine = max(1, int(n_refine))
        self.names = list(self.model.names)
        self.free = np.array([k not in self.model.fixed for k in self.names])
        zc = np.array([geometry.zc_au(y) for _, y, _ in self.entries])
        self.z_cuts = zc
        b = self.model.default_bounds(geometry, zc)
        self.lo = np.array([b[k][0] for k in self.names], float)
        self.hi = np.array([b[k][1] for k in self.names], float)
        # a fixed parameter keeps its value, even outside the default bounds (e.g. a
        # break height fixed from a fit to more cuts than this one)
        for i, k in enumerate(self.names):
            if k in self.model.fixed:
                v = float(self.model.fixed[k])
                self.lo[i], self.hi[i] = min(self.lo[i], v), max(self.hi[i], v)

    @classmethod
    def from_config(cls, model, geometry, traces, sig_x, fit_config, verbose=True,
                    counter_traces=None):
        """Fitter with the options of a :class:`~arcwind.config.FitConfig`."""
        f = fit_config
        return cls(model, geometry, traces, sig_x, loss=f.loss,
                   equal_cut_weight=f.equal_cut_weight, max_nfev=f.max_nfev, nphi=f.nphi,
                   scale_errors_by_chi2=f.scale_errors_by_chi2, check_vr_sign=f.check_vr_sign,
                   j_multistart=f.j_multistart, seed=f.seed, verbose=verbose,
                   counter_traces=counter_traces, error_method=f.error_method,
                   n_jobs=f.n_jobs, resample_tol=f.resample_tol, n_refine=f.n_refine)

    # ----------------------------------------------------------- internals
    def _log(self, msg):
        if self.verbose:
            print(f"[fit] {msg}")

    def _full(self, q, base):
        p = base.copy()
        p[self.free] = q
        return p

    def _curve(self, gi, y0, p):
        """Model curve of a cut of lobe ``gi`` (parameters given in the main frame)."""
        g, mirror = self.groups[gi]
        return self.model.curve(y0, g, p * self.msign if mirror else p, nphi=self.nphi)

    def _resid(self, p, entries):
        w = cut_weights([t for _, _, t in entries], self.equal)
        out = []
        for k, (gi, y0, tr) in enumerate(entries):
            r = odr_distances(self._curve(gi, y0, p), tr, self.sig_x)
            r = np.full(tr.n, 1e3) if r is None else r
            out.append(r * (w[k] if w is not None else 1.0))
        return np.concatenate(out)

    def residuals(self, p):
        """Weighted distances of all trace points (all cuts, all lobes) for the full
        parameter vector ``p``."""
        return self._resid(p, self.entries)

    def _run(self, entries, starts, lo, hi, base):
        return self._run_list(entries, starts, base, lo, hi)

    def _vec(self, d):
        v = self.model.vector()
        for i, k in enumerate(self.names):
            if k in d:
                v[i] = d[k]
            if k in self.model.fixed:
                v[i] = self.model.fixed[k]
        return v

    def velocity_guess(self):
        """Typical axial and radial speeds of the traces (km/s): the visible half of
        each arch is the one farthest from Vsys and its legs are near the ring-center
        velocity."""
        g = self.geom
        trs = self.traces or self.counter_traces
        if not self.traces:
            g = g.mirrored()
        v = np.concatenate([t.v for t in trs.values()]) - g.vsys
        lo, hi = np.percentile(v, 5), np.percentile(v, 95)
        vc, A = (lo, hi - lo) if abs(hi) >= abs(lo) else (hi, hi - lo)
        vr = max(A, 0.3) / g.sini
        vz = -g.s * vc / g.cosi if not g.edge_on else vr
        return float(vz), float(vr)

    # ----------------------------------------------------------------- fit
    def fit(self, init=None, start_from=None, n_boot=0):
        """Fit the model.

        Parameters
        ----------
        init : dict, optional
            Extra starting values (tried first).
        start_from : FitResult or ConicalWindModel, optional
            A previous (e.g. simpler) fit: start from its geometry and rotation, with
            starting points of this velocity law that follow its velocity profile
            (no cone-angle grid or j restarts).
        n_boot : int
            Bootstrap resamplings of the cuts (0: covariance errors only).

        Returns
        -------
        FitResult
        """
        m, g = self.model, self.geom
        base = self._vec(init or {})
        lo, hi = self.lo.copy(), self.hi.copy()
        base = np.clip(base, lo, hi)
        if start_from is not None:
            starts = [self._vec(s) for s in m.starts(g, self.z_cuts, 0, 0,
                                                     previous=start_from)]
            if init:
                starts.insert(0, base)
            fresh = False
        else:
            vz, vr = self.velocity_guess()
            xs = np.concatenate([t.x if gi == 0 else -t.x for gi, _, t in self.entries])
            xref = float(np.median(xs) * g.distance_pc)
            starts = [base]
            if "tan_theta" not in m.fixed:
                starts += [self._vec(s) for s in m.starts(g, self.z_cuts, vz, xref,
                                                          vr_ref=vr)]
            fresh = True
        npts = sum(t.n for _, _, t in self.entries)
        lobes = (f" ({len(self.traces)} main + {len(self.counter_traces)} counter lobe)"
                 if self.counter_traces else "")
        self._log(f"{m!r}: {len(self.entries)} cuts{lobes}, {npts} points, {self.free.sum()} "
                  f"free parameters, {len(starts)} starting points")
        best = self._run(self.entries, starts, lo, hi, base)
        if best is None:
            raise RuntimeError("the optimizer failed for every starting point")
        if fresh and self.j_multistart and "j" not in m.fixed:
            p1 = self._full(best.x, base)
            P1 = dict(zip(self.names, p1))
            js = 0.5 * (P1["R0_au"] + m.zref_au * P1["tan_theta"]) + 20.0
            st2 = []
            for jj in (js, -js, 3 * js, -3 * js):
                st = p1.copy()
                st[self.names.index("j")] = jj
                st2.append(st)
            b2 = self._run(self.entries, st2, lo, hi, base)
            if b2 is not None and b2.cost < best.cost:
                best = b2
        vr_check = None
        if (m.radial == "proportional" and self.check_vr_sign and "Vr" not in m.fixed
                and "Vr" not in m.bounds):
            iv = self.names.index("Vr")
            if self._full(best.x, base)[iv] < 0:
                lo2 = lo.copy()
                lo2[iv] = 0.0
                st = self._full(best.x, base).copy()
                st[iv] = abs(st[iv])
                alt = self._run(self.entries, [st] + [s for s in starts if s[iv] >= 0], lo2,
                                hi, base)
                if alt is not None:
                    vr_check = dict(chi2_neg=float(2 * best.cost), chi2_pos=float(2 * alt.cost),
                                    p_neg=self._full(best.x, base).tolist(),
                                    p_pos=self._full(alt.x, base).tolist())
                    self._log(f"Vr < 0 in the free fit (gas moving toward the axis): "
                              f"chi2(Vr<0) = {vr_check['chi2_neg']:.1f}, chi2(Vr>=0) = "
                              f"{vr_check['chi2_pos']:.1f}")
                    if self.check_vr_sign == "prefer_positive" or alt.cost <= best.cost:
                        best = alt
        p = self._full(best.x, base)
        res = self._result(p, best)
        res.vr_check = vr_check
        self._log(f"chi2 = {res.chi2:.1f} for {res.n} points (chi2_red = {res.chi2r:.2f})")
        if n_boot:
            self.bootstrap(res, n_boot)
        return res

    def _result(self, p, best):
        r = self.residuals(p)
        n, k = len(r), int(self.free.sum())
        chi2 = float(np.sum(r ** 2))
        chi2r = chi2 / max(n - k, 1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cov = np.linalg.pinv(best.jac.T @ best.jac)
        if self.scale and chi2r > 1:
            cov = cov * chi2r
        err = np.zeros(len(p))
        err[self.free] = np.sqrt(np.clip(np.diag(cov), 0, None))
        amp = self.model.law.amplitude
        if not self.geom.edge_on and amp not in self.model.fixed:
            ia = self.names.index(amp)
            err[ia] = np.hypot(err[ia], self.geom.vsys_err / self.geom.cosi)
        m = self.model.copy()
        m.set(**{kk: v for kk, v in zip(self.names, p) if kk not in m.fixed})
        per_cut, per_cut_c = {}, {}
        for gi, y0, tr in self.entries:
            d = odr_distances(self._curve(gi, y0, p), tr, self.sig_x)
            (per_cut if gi == 0 else per_cut_c)[y0] = \
                float(np.mean(d ** 2)) if d is not None else np.nan
        return FitResult(model=m, geometry=self.geom, err=dict(zip(self.names, err)),
                         cov=cov, free=[kk for kk, f in zip(self.names, self.free) if f],
                         chi2=chi2, chi2r=chi2r, n=n, traces=dict(self.traces),
                         other_traces=dict(self.counter_traces),
                         sig_x=self.sig_x, per_cut_d2=per_cut, other_per_cut_d2=per_cut_c,
                         settings=dict(loss=self.loss, equal_cut_weight=self.equal,
                                       error_method=self.error_method,
                                       nphi=self.nphi))

    def bootstrap(self, result, n_boot):
        """Uncertainties from refits to resampled cuts; fills ``result.boot``
        (samples x k) in place.

        ``error_method = 'bootstrap'``: ``n_boot`` refits to cuts drawn with
        replacement. ``'jackknife'``: one refit per left-out cut (``n_boot`` only
        switches it on); the samples are stored as pseudo-values,
        ``best + sqrt(n - 1) (p_i - mean)``, whose spread is the jackknife error,
        so every later use (errors, bands, histograms) works unchanged.

        Every refit starts from the best fit, with a looser tolerance
        (``resample_tol``), and the refits run in parallel on ``n_jobs`` processes
        (-1: all cores; 1: serial). If the parallel run fails, they run serially.
        The Vsys uncertainty is added to the velocity amplitude of every sample.
        """
        ent = self.entries
        if len(ent) < 3:
            self._log("resampling needs >= 3 cuts; skipped")
            return result
        p = result.model.vector()
        n = len(ent)
        if self.error_method == "jackknife":
            samples = [np.delete(np.arange(n), k) for k in range(n)]
        else:
            samples = [self.rng.integers(0, n, n) for _ in range(int(n_boot))]
        bs = self._resample(samples, p)
        if len(bs) and self.error_method == "jackknife":
            bs = p[None, :] + np.sqrt(len(bs) - 1) * (bs - bs.mean(axis=0))
        if len(bs) and not self.geom.edge_on:
            ia = self.names.index(self.model.law.amplitude)
            if self.model.law.amplitude not in self.model.fixed:
                bs[:, ia] += self.rng.normal(0, self.geom.vsys_err / self.geom.cosi, len(bs))
        result.boot = bs
        result.settings["error_method"] = self.error_method
        return result

    def _refit(self, idx, p):
        # duplicate cuts are kept as separate entries (cuts of both lobes pooled)
        rb = self._run_list([self.entries[j] for j in idx], [p], p, tol=self.resample_tol)
        return None if rb is None else self._full(rb.x, p)

    def _resample(self, samples, p):
        """Refit every index set in ``samples``; in parallel when possible."""
        nproc = self._nproc(len(samples))
        out = None
        if nproc > 1:
            try:
                done = _parallel_map(self, "_refit", [(idx, p) for idx in samples], nproc,
                                     progress=self.error_method)
                out = [done[i] for i in range(len(samples))]
                self._log(f"{self.error_method}: {len(samples)} refits on {nproc} processes")
            except Exception as e:                      # fall back to serial refits
                self._log(f"parallel refits failed ({type(e).__name__}: {e}); "
                          f"running them serially")
        if out is None:
            out = []
            for k, idx in enumerate(samples):
                out.append(self._refit(idx, p))
                if (k + 1) % 10 == 0:
                    self._log(f"{self.error_method} {k + 1}/{len(samples)}")
        return np.array([r for r in out if r is not None])

    def _lsq(self, entries, q0, base, lo, hi, tol, max_nfev):
        """One bounded least-squares solve from the free-parameter vector ``q0``."""
        fr = self.free

        def resid(q):
            return self._resid(self._full(q, base), entries)
        q0 = np.clip(q0, lo[fr] + 1e-6, hi[fr] - 1e-6)
        try:
            return least_squares(resid, q0, bounds=(lo[fr], hi[fr]), loss=self.loss,
                                 x_scale="jac", max_nfev=max_nfev, ftol=tol, xtol=tol, gtol=tol)
        except ValueError:
            return None

    def _screen_one(self, entries, q0, base, lo, hi):
        """Coarse solve of one starting point (multistart screening)."""
        return self._lsq(entries, q0, base, lo, hi, self.screen_tol,
                         min(self.max_nfev, self.screen_nfev))

    def _run_list(self, entries, starts, base, lo=None, hi=None, tol=1e-8):
        """Best least-squares solution over ``starts`` for a list of
        (lobe, y0, trace) (cuts may repeat, as in bootstrap samples).

        With several starts, every start is first solved coarsely (tolerance
        ``screen_tol``, at most ``screen_nfev`` evaluations; in parallel on
        ``n_jobs`` processes when possible), and only the ``n_refine`` best are
        solved to the full tolerance."""
        lo = self.lo if lo is None else lo
        hi = self.hi if hi is None else hi
        fr = self.free
        qs = [np.asarray(st, float)[fr] for st in starts]
        if len(qs) > 1:
            coarse = self._map("_screen_one", [(entries, q, base, lo, hi) for q in qs])
            coarse = sorted([r for r in coarse if r is not None], key=lambda r: r.cost)
            qs = [r.x for r in coarse[:self.n_refine]]
        if len(qs) > 1:
            fine = self._map("_lsq", [(entries, q, base, lo, hi, tol, self.max_nfev)
                                      for q in qs])
        else:
            fine = [self._lsq(entries, q, base, lo, hi, tol, self.max_nfev) for q in qs]
        best = None
        for r in fine:
            if r is not None and (best is None or r.cost < best.cost):
                best = r
        return best

    def _nproc(self, n_tasks):
        import os
        n = self.n_jobs if self.n_jobs and self.n_jobs > 0 else (os.cpu_count() or 1)
        return max(1, min(n, n_tasks))

    def _map(self, method, arglist):
        """``[self.<method>(*args) for args in arglist]``, in parallel when possible
        (forked processes; serial if that fails)."""
        nproc = self._nproc(len(arglist))
        if nproc > 1 and not _WORKER:                # never nest pools
            try:
                done = _parallel_map(self, method, arglist, nproc)
                return [done[i] for i in range(len(arglist))]
            except Exception as e:
                self._log(f"parallel run failed ({type(e).__name__}: {e}); running serially")
        return [getattr(self, method)(*a) for a in arglist]


_WORKER = {}


def _worker_call(i):
    f = _WORKER["fitter"]
    return i, getattr(f, _WORKER["method"])(*_WORKER["args"][i])


def _parallel_map(fitter, method, arglist, nproc, progress=None):
    """Run ``fitter.<method>(*args)`` for every ``args`` in ``arglist`` on ``nproc``
    forked processes (the fitter and arguments are inherited, not pickled).
    Returns ``{index: result}``. Raises if fork is unavailable (e.g. Windows) or a
    worker fails."""
    import multiprocessing as mp
    import warnings as _w
    from concurrent.futures import ProcessPoolExecutor, as_completed
    if "fork" not in mp.get_all_start_methods():
        raise RuntimeError("the 'fork' start method is not available")
    _WORKER.update(fitter=fitter, method=method, args=arglist)
    done = {}
    try:
        with _w.catch_warnings(), \
                ProcessPoolExecutor(max_workers=nproc, mp_context=mp.get_context("fork")) as ex:
            _w.simplefilter("ignore", DeprecationWarning)   # fork with threads (py >= 3.12)
            futs = [ex.submit(_worker_call, i) for i in range(len(arglist))]
            for k, fu in enumerate(as_completed(futs)):
                i, r = fu.result()
                done[i] = r
                if progress and (k + 1) % 10 == 0:
                    fitter._log(f"{progress} {k + 1}/{len(arglist)}")
    finally:
        _WORKER.clear()
    return done


class FitResult:
    """Result of a global fit.

    Attributes
    ----------
    model : ConicalWindModel
        The model at the best-fit values (with its fixed parameters and bounds).
    geometry : ViewGeometry
    err : dict
        Covariance uncertainties (including the Vsys error on the velocity amplitude).
    boot : ndarray or None
        Bootstrap samples, shape (n_boot, k), columns in ``names`` order.
    free : list of str
        Fitted parameters.
    chi2, chi2r : float
        Sum of squared (weighted) distances, and per degree of freedom.
    n : int
        Number of trace points.
    traces : dict
        ``{y0: Trace}`` used in the fit (lobe ``frame``).
    other_traces : dict
        Traces of the opposite lobe in a two-lobe fit (in that lobe's own frame).
    frame : {'main', 'counter'}
        The lobe whose frame the parameters refer to.
    per_cut_d2 : dict
        Mean squared distance of the points of each cut (goodness of fit per cut).
    vr_check : dict or None
        Outcome of the Vr-sign test.

    A two-lobe fit describes both lobes with one model. ``res.view("counter")``
    gives the same result in the counter lobe's frame (mirrored model and geometry,
    counter-lobe traces), which every plotting function accepts with the
    counter-lobe cube.

    Examples
    --------
    >>> res = wa.fit()                                     # doctest: +SKIP
    >>> res.params["tan_theta"], res.errors()["j"]         # doctest: +SKIP
    >>> res.theta_deg                                      # doctest: +SKIP
    (18.9, 0.4)
    >>> prof = res.profiles(np.linspace(50, 600, 100))     # with 16-84% bands  # doctest: +SKIP
    >>> res.save("bpl_fit.json"); FitResult.load("bpl_fit.json")              # doctest: +SKIP
    """

    def __init__(self, model, geometry, err, free, chi2, chi2r, n, cov=None, boot=None,
                 traces=None, sig_x=None, per_cut_d2=None, settings=None, vr_check=None,
                 label=None, other_traces=None, other_per_cut_d2=None, frame="main"):
        self.model = model
        self.geometry = geometry
        self.err = dict(err)
        self.cov = cov
        self.free = list(free)
        self.chi2 = float(chi2)
        self.chi2r = float(chi2r)
        self.n = int(n)
        self.boot = boot
        self.traces = traces or {}
        self.sig_x = sig_x
        self.per_cut_d2 = per_cut_d2 or {}
        self.settings = settings or {}
        self.vr_check = vr_check
        self.label = label or model.law.name
        self.other_traces = other_traces or {}
        self.other_per_cut_d2 = other_per_cut_d2 or {}
        self.frame = frame

    # ------------------------------------------------------------- lobes
    @property
    def lobes(self):
        """'both' for a two-lobe fit, else the lobe that was fitted."""
        if self.traces and self.other_traces:
            return "both"
        if self.other_traces and not self.traces:
            return "counter" if self.frame == "main" else "main"
        return self.frame

    def view(self, lobe):
        """This result in the frame of ``lobe`` ('main' or 'counter').

        For the other lobe, the model and geometry are mirrored
        (:meth:`ConicalWindModel.mirrored`, :meth:`ViewGeometry.mirrored`), the
        bootstrap samples and covariance are transformed accordingly, and the traces
        are those of that lobe. Use it to plot the counter lobe of a two-lobe fit:
        ``plot_cuts(counter_cube, res.view("counter"))``.
        """
        if lobe not in ("main", "counter"):
            raise ValueError("lobe must be 'main' or 'counter'")
        if lobe == self.frame:
            return self
        sg = self.model.mirror_signs()                       # +1/-1 per parameter
        sg_free = sg[[self.names.index(k) for k in self.free]]
        cov = None if self.cov is None else self.cov * np.outer(sg_free, sg_free)
        return FitResult(model=self.model.mirrored(), geometry=self.geometry.mirrored(),
                         err=self.err, free=self.free, chi2=self.chi2, chi2r=self.chi2r,
                         n=self.n, cov=cov,
                         boot=None if self.boot is None else np.asarray(self.boot) * sg,
                         traces=self.other_traces, other_traces=self.traces,
                         sig_x=self.sig_x, per_cut_d2=self.other_per_cut_d2,
                         other_per_cut_d2=self.per_cut_d2, settings=self.settings,
                         vr_check=self.vr_check, label=self.label, frame=lobe)

    # ------------------------------------------------------------ values
    @property
    def names(self):
        return list(self.model.names)

    @property
    def params(self):
        """Best-fit values ``{name: value}``."""
        return dict(self.model.values)

    @property
    def k(self):
        """Number of free parameters."""
        return len(self.free)

    @property
    def bic(self):
        """Bayesian information criterion, chi2 + k ln n."""
        return self.chi2 + self.k * np.log(max(self.n, 1))

    @property
    def boot_err(self):
        """Bootstrap standard deviations ``{name: err}`` (or None)."""
        if self.boot is None or not len(self.boot):
            return None
        return dict(zip(self.names, np.std(self.boot, axis=0)))

    def errors(self, kind="max"):
        """Uncertainties: ``'cov'``, ``'boot'`` or ``'max'`` (the larger of the two)."""
        be = self.boot_err
        if kind == "cov" or be is None:
            return {k: float(v) for k, v in self.err.items()}
        if kind == "boot":
            return {k: float(v) for k, v in be.items()}
        return {k: float(max(self.err[k], be[k])) for k in self.names}

    @property
    def theta_deg(self):
        """Cone half-opening angle and its uncertainty (deg)."""
        T = self.params["tan_theta"]
        e = self.errors()["tan_theta"]
        return float(np.degrees(np.arctan(T))), float(np.degrees(e / (1 + T ** 2)))

    def draws(self, n=300, seed=42):
        """Parameter samples (n, k): the bootstrap if available, else Gaussian draws
        from the covariance errors."""
        if self.boot is not None and len(self.boot):
            return np.asarray(self.boot)
        rng = np.random.default_rng(seed)
        p = self.model.vector()
        e = np.array([self.err[k] for k in self.names])
        return p + e * rng.standard_normal((n, len(p)))

    # ------------------------------------------------------- evaluations
    def curve(self, y0, p=None, nphi=361):
        """Model curve of the cut at ``y0`` (see :meth:`ConicalWindModel.curve`)."""
        return self.model.curve(y0, self.geometry, p, nphi=nphi)

    def ring_at(self, y0, p=None):
        """Ring quantities of the model at the cut ``y0``."""
        return self.model.ring_at(y0, self.geometry, p)

    def velocities(self, z, p=None):
        """(Vz, Vr) at heights ``z`` (au)."""
        return self.model.velocities(z, p)

    def profiles(self, z, band=(16, 84)):
        """Model quantities along the flow with uncertainty bands.

        Returns
        -------
        dict
            For each quantity of :meth:`ConicalWindModel.profiles` the best value,
            and ``<name>_lo`` / ``<name>_hi`` percentiles over :meth:`draws`.
        """
        z = np.asarray(z, float)
        best = self.model.profiles(z)
        D = self.model.profiles(z, self.draws())
        out = dict(best)
        for k, v in D.items():
            if k == "z":
                continue
            lo, hi = np.nanpercentile(v, band, axis=0)
            out[k + "_lo"], out[k + "_hi"] = lo, hi
        return out

    def chi2_on(self, traces, equal_cut_weight=None, nphi=None, other_traces=None):
        """chi2 of this model on other traces (e.g. to compare models on the same
        points); ``other_traces`` are cuts of the opposite lobe."""
        eq = self.settings.get("equal_cut_weight", False) if equal_cut_weight is None \
            else equal_cut_weight
        nphi = self.settings.get("nphi", 181) if nphi is None else nphi
        tl = [(self, y, t) for y, t in sorted(traces.items())]
        if other_traces:
            o = self.view("counter" if self.frame == "main" else "main")
            tl += [(o, y, t) for y, t in sorted(other_traces.items())]
        w = cut_weights([t for _, _, t in tl], eq)
        tot = 0.0
        for k, (res, y0, tr) in enumerate(tl):
            r = odr_distances(res.curve(y0, nphi=nphi), tr, self.sig_x)
            r = np.full(tr.n, 1e3) if r is None else r
            tot += float(np.sum((r * (w[k] if w is not None else 1.0)) ** 2))
        return tot

    def derived(self, mass=None, n_draw=300):
        """Derived physical quantities (:class:`~arcwind.derived.WindProperties`)."""
        from .derived import WindProperties
        return WindProperties(self, mass=mass)

    # ------------------------------------------------------------ report
    def summary(self):
        """Text summary (parameters, errors, fit statistics, rotation sense)."""
        g = self.geometry
        e, be = self.err, self.boot_err
        lob = {"both": f"both lobes ({len(self.traces)} + {len(self.other_traces)} cuts), "
                       f"parameters in the {self.frame} lobe's frame",
               "main": "main lobe", "counter": "counter lobe"}[self.lobes]
        L = [f"# {self.label}: {self.model!r}; {lob}",
             f"# {len(self.traces) + len(self.other_traces)} cuts, {self.n} points; "
             f"chi2 = {self.chi2:.1f}, "
             f"chi2_red = {self.chi2r:.2f}, BIC = {self.bic:.1f} (k = {self.k})",
             f"# i = {g.inclination_deg:.2f} deg, s = {g.s:+d}, Vsys = {g.vsys:.3f} +/- "
             f"{g.vsys_err} km/s, disk sense = {g.disk_sense:+d}",
             "# errors: cov = covariance (x sqrt(chi2_red) if > 1; + Vsys on the amplitude); "
             f"boot = {self.settings.get('error_method', 'bootstrap')} over cuts "
             "(quote the larger)"]
        for k in self.names:
            v = self.params[k]
            if k in self.model.fixed:
                L.append(f"{k:10s} = {v:10.4g}   (fixed)")
                continue
            s = f"{k:10s} = {v:10.4g}   cov +/-{e[k]:.3g}"
            if be is not None:
                s += f"   boot +/-{be[k]:.3g}"
            L.append(s)
        th, eth = self.theta_deg
        L.append(f"theta      = {th:.2f} +/- {eth:.2f} deg")
        if self.boot is not None and len(self.boot) and "j" in self.free:
            jb = self.boot[:, self.names.index("j")]
            L.append(f"sign of j kept in {100 * np.mean(np.sign(jb) == np.sign(self.params['j'])):.0f}"
                     f"% of {self.settings.get('error_method', 'bootstrap')} samples")
        j = self.params["j"]
        L.append(f"wind rotation (j {'>' if j > 0 else '<'} 0) is in the "
                 f"{'SAME' if j * g.disk_sense > 0 else 'OPPOSITE'} sense as the disk")
        if self.vr_check:
            L.append(f"Vr sign check: chi2(Vr<0) = {self.vr_check['chi2_neg']:.1f}, "
                     f"chi2(Vr>=0) = {self.vr_check['chi2_pos']:.1f}")
        L.append("NB: trace points are ~2 per beam (trace.dx_arcsec): chi2/BIC differences "
                 "are overestimated by about that factor.")
        return "\n".join(L)

    def table(self, fmt="text", kind="max"):
        """Parameter table: ``'text'``, ``'markdown'`` or ``'latex'`` rows."""
        err = self.errors(kind)
        rows = []
        for k in self.names:
            st = "fixed" if k in self.model.fixed else f"{err[k]:.3g}"
            rows.append((k, f"{self.params[k]:.4g}", st))
        if fmt == "markdown":
            return "\n".join(["| parameter | value | error |", "|---|---|---|"] +
                             [f"| {a} | {b} | {c} |" for a, b, c in rows])
        if fmt == "latex":
            esc = [(a.replace("_", "\\_"), b, c) for a, b, c in rows]
            return "\n".join(f"{a} & {b} & {c} \\\\" for a, b, c in esc)
        return "\n".join(f"{a:10s} {b:>12s} {c:>10s}" for a, b, c in rows)

    def __repr__(self):
        th, eth = self.theta_deg
        return (f"FitResult({self.model.law.name}, radial={self.model.radial}: "
                f"theta={th:.1f}+/-{eth:.1f} deg, chi2_red={self.chi2r:.2f}, n={self.n}, "
                f"boot={0 if self.boot is None else len(self.boot)})")

    # ------------------------------------------------------------- I/O
    def to_dict(self, include_traces=True):
        m = self.model
        d = dict(format="arcwind.FitResult/1", label=self.label,
                 model=dict(velocity_law=m.law.name, radial=m.radial, rotation=m.rotation,
                            zref_au=m.zref_au, bpl_delta=m.law.delta, values=m.values,
                            fixed=m.fixed, bounds=m.bounds),
                 geometry=self.geometry.as_dict(), names=self.names, err=self.err,
                 free=self.free, chi2=self.chi2, chi2r=self.chi2r, n=self.n, bic=self.bic,
                 boot=None if self.boot is None else np.asarray(self.boot).tolist(),
                 sig_x=self.sig_x, per_cut_d2={f"{k:.4f}": v for k, v in self.per_cut_d2.items()},
                 settings=self.settings, vr_check=self.vr_check, frame=self.frame,
                 other_per_cut_d2={f"{k:.4f}": v for k, v in self.other_per_cut_d2.items()})
        if include_traces:
            for key, trs in (("traces", self.traces), ("other_traces", self.other_traces)):
                d[key] = {f"{y:.4f}": {k: np.asarray(getattr(t, k)).tolist()
                                       for k in Trace.FIELDS}
                          for y, t in trs.items()}
        return d

    def save(self, path, include_traces=True):
        """Write the result (including bootstrap samples and traces) to JSON."""
        with open(path, "w") as fh:
            json.dump(_jsonable(self.to_dict(include_traces)), fh, indent=1)
        return path

    @classmethod
    def load(cls, path):
        """Read a result written by :meth:`save`."""
        with open(path) as fh:
            d = json.load(fh)
        md = d["model"]
        m = ConicalWindModel(md["velocity_law"], radial=md["radial"], rotation=md["rotation"],
                             zref_au=md["zref_au"], bpl_delta=md["bpl_delta"])
        m.set(**md["values"])
        m.fixed = {k: float(v) for k, v in md["fixed"].items()}
        m.bounds = {k: tuple(v) for k, v in md["bounds"].items()}
        g = ViewGeometry(**d["geometry"])
        def _tr(key):
            return {float(y): Trace(float(y), *(np.asarray(t[k], int if k == "side" else float)
                                                for k in Trace.FIELDS if k in t))
                    for y, t in (d.get(key) or {}).items()}
        tr = _tr("traces")
        return cls(model=m, geometry=g, err=d["err"], free=d["free"], chi2=d["chi2"],
                   chi2r=d["chi2r"], n=d["n"],
                   boot=None if d.get("boot") is None else np.asarray(d["boot"], float),
                   traces=tr, sig_x=d.get("sig_x"),
                   per_cut_d2={float(k): v for k, v in (d.get("per_cut_d2") or {}).items()},
                   settings=d.get("settings"), vr_check=d.get("vr_check"), label=d.get("label"),
                   other_traces=_tr("other_traces"), frame=d.get("frame", "main"),
                   other_per_cut_d2={float(k): v for k, v in
                                     (d.get("other_per_cut_d2") or {}).items()})


def compare_fits(results, traces=None, other_traces=None):
    """Compare fits on the same trace points (chi2, BIC, Delta BIC).

    Parameters
    ----------
    results : dict or list of FitResult
        ``{label: result}`` or a list (labels from ``result.label``).
    traces : dict, optional
        Common traces; default: those of the last result.
    other_traces : dict, optional
        Common traces of the opposite lobe (two-lobe fits); default: those of the
        last result when ``traces`` is not given.

    Returns
    -------
    list of dict
        One row per model: ``label, k, n, chi2, bic, dbic`` (``dbic`` relative to
        the best BIC), sorted by BIC.
    """
    items = list(results.items()) if isinstance(results, dict) else \
        [(r.label, r) for r in results]
    last = items[-1][1]
    tr = traces or last.traces
    other = other_traces if traces is not None else (other_traces or last.other_traces or None)
    n = sum(t.n for t in tr.values()) + sum(t.n for t in (other or {}).values())
    rows = []
    for lab, r in items:
        c2 = r.chi2_on(tr, other_traces=other)
        rows.append(dict(label=lab, k=r.k, n=n, chi2=c2, bic=c2 + r.k * np.log(n)))
    best = min(x["bic"] for x in rows)
    for x in rows:
        x["dbic"] = x["bic"] - best
    return sorted(rows, key=lambda x: x["bic"])


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    if isinstance(o, float) and not np.isfinite(o):
        return None if np.isnan(o) else (1e308 if o > 0 else -1e308)
    return o
