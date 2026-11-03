"""
:class:`WindAnalysis` -- the high-level interface: one object that holds the
configuration, the cube, the traces and the fit results, with one method per
step of the analysis.

Typical session
---------------
>>> from diskwind import WindAnalysis
>>> wa = WindAnalysis("elias20.toml")                       # doctest: +SKIP
>>> wa.check_noise().summary()          # 1. noise channels line-free?      # doctest: +SKIP
>>> wa.check_disk_rotation().summary()  # 2. disk rotation sense             # doctest: +SKIP
>>> pl = wa.fit("power-law", n_boot=0)  # 3. quick first fit                 # doctest: +SKIP
>>> bpl = wa.fit("broken-power-law", guide=pl, start_from=pl, iterations=2)  # doctest: +SKIP
>>> wa.compare(pl, bpl)                 # 4. chi2 / BIC on the same traces   # doctest: +SKIP
>>> m = wa.mass(bpl)                    # 5. mass and mass-loss rate         # doctest: +SKIP
>>> print(wa.derived(bpl, m).summary()) # 6. Vp, r_L, travel time, ...       # doctest: +SKIP
>>> wa.plot_cuts(bpl, others=[pl])      # figures                            # doctest: +SKIP

Or everything at once, as set in the configuration: ``wa.run()`` (also
``diskwind run elias20.toml`` on the command line).
"""
from __future__ import annotations

import os

import numpy as np

from . import plotting
from .config import WindConfig
from .constants import FWHM2SIG
from .cube import WindCube
from .fitting import FitResult, GlobalFitter, compare_fits
from .mass import MassAnalysis
from .model import ConicalWindModel
from .ringfit import RingFitter
from .trace import Trace, TraceExtractor

__all__ = ["WindAnalysis"]


class WindAnalysis:
    """A complete wind analysis.

    Parameters
    ----------
    config : WindConfig or str
        Configuration object or file (TOML / YAML / JSON / Python).
    cube : WindCube, optional
        An already loaded cube (default: loaded on first use).
    **options
        Configuration options to change (bare or dotted names), e.g.
        ``WindAnalysis("elias20.toml", n_boot=0, **{"output.dir": "test"})``.

    Attributes
    ----------
    config : WindConfig
    results : dict
        Fit results by label, in the order they were made.
    traces : dict
        Traces of the last extraction, ``{y0: Trace}``.
    """

    def __init__(self, config, cube=None, **options):
        cfg = config if isinstance(config, WindConfig) else WindConfig.load(config)
        if options:
            cfg = cfg.copy().set(**options)
        self.config = cfg
        self._cube = cube
        self.results = {}
        self.traces = {}

    # ---------------------------------------------------------------- data
    @property
    def cube(self):
        """The :class:`WindCube` (loaded on first access)."""
        if self._cube is None:
            self._cube = WindCube(self.config)
        return self._cube

    @property
    def geometry(self):
        return self.cube.geometry

    @property
    def sig_x(self):
        """Position uncertainty of trace points (arcsec)."""
        return self.config.trace.sigx_beam * self.cube.beam_x * FWHM2SIG

    def _log(self, msg):
        if self.config.output.verbose:
            print(f"[diskwind] {msg}")

    def pv(self, y0, **kw):
        """Transverse PV diagram (:meth:`WindCube.pv`)."""
        return self.cube.pv(y0, **kw)

    def check_noise(self, plot=True):
        """Noise-channel check (:func:`diskwind.diagnostics.check_noise`); saves
        ``noise_check.txt`` and the figure."""
        rep = self.cube.check_noise()
        with open(self.config.outpath("noise_check", "txt"), "w") as fh:
            fh.write(rep.summary() + "\n")
        if plot:
            rep.plot(save=self.config.outpath("noise_check"))
        self._log(rep.summary().splitlines()[-1])
        return rep

    def check_disk_rotation(self, y_offsets=(0.0, -0.15, -0.3), plot=True):
        """Disk rotation-sense check (:func:`diskwind.diagnostics.check_disk_rotation`)."""
        chk = self.cube.check_disk_rotation(y_offsets)
        if plot:
            chk.plot(save=self.config.outpath("disk_rotation_check"))
        self._log(chk.summary())
        return chk

    # -------------------------------------------------------------- traces
    def extract_traces(self, y_offsets=None, guide=None, **trace_options):
        """Extract the traces of the cuts (:class:`TraceExtractor`).

        Parameters
        ----------
        y_offsets : array-like, optional
            Default: ``cuts.y_offsets_arcsec``.
        guide : FitResult, str or None
            Guided extraction around a fit (or the label of a stored result).
        **trace_options
            ``[trace]`` options for this extraction.

        Returns
        -------
        dict
            ``{y0: Trace}`` (also stored in ``self.traces``).
        """
        g = self._resolve(guide)
        self.traces = TraceExtractor(self.cube).extract_all(y_offsets, guide=g, **trace_options)
        n = sum(t.n for t in self.traces.values())
        self._log(f"traces: {len(self.traces)} cuts, {n} points"
                  + (f" (guided by {g.label})" if g is not None else ""))
        return self.traces

    def _resolve(self, r):
        if r is None or isinstance(r, (FitResult, ConicalWindModel)):
            return r
        if r in ("last", "previous"):
            return list(self.results.values())[-1] if self.results else None
        return self.results[r]

    # ----------------------------------------------------------------- fit
    def model(self, velocity_law=None, radial=None, rotation=None, fixed=None, bounds=None,
              **values):
        """A :class:`ConicalWindModel` as set by ``[fit]``, with optional changes."""
        f = self.config.fit
        m = ConicalWindModel(velocity_law or f.velocity_law, radial=radial or f.radial,
                             rotation=rotation or f.rotation, zref_au=f.zref_au,
                             bpl_delta=f.bpl_delta)
        for k, v in ({**(f.fixed or {}), **(fixed or {})}).items():
            if k in m.names:
                m.fix(**{k: v})
        for k, b in ({**(f.bounds or {}), **(bounds or {})}).items():
            if k in m.names:
                m.set_bounds(**{k: b})
        return m.set(**values)

    def fit(self, velocity_law=None, radial=None, rotation=None, fixed=None, bounds=None,
            guide=None, start_from=None, iterations=None, n_boot=None, y_offsets=None,
            init=None, label=None, traces=None, save=True, **fit_options):
        """Global fit of a conical wind to all cuts.

        Parameters
        ----------
        velocity_law, radial, rotation : str, optional
            Model choices (default: ``[fit]``); see :class:`ConicalWindModel`.
        fixed : dict, optional
            Extra fixed parameters ``{name: value}``.
        bounds : dict, optional
            Extra bounds ``{name: (lo, hi)}``.
        guide : FitResult or str, optional
            Extract the traces around this fit (``'last'``: the previous result).
        start_from : FitResult or str, optional
            Start from this (simpler) fit, e.g. a power law for a broken power law.
        iterations : int, optional
            With a guide: re-extract the traces around the new model and refit this
            many times in total (default ``fit.iterations``).
        n_boot : int, optional
            Bootstrap resamplings (default ``fit.n_boot``; done on the last iteration).
        y_offsets : array-like, optional
            Cuts to use.
        init : dict, optional
            Starting values.
        label : str, optional
            Name of the result (default: the velocity law, plus the radial mode if not
            'cone' and '-jz' for a power-law j; a number is appended if that name is
            already used, so earlier fits are never overwritten).
        traces : dict, optional
            Use these traces instead of extracting them.
        save : bool
            Write ``<label>.json`` and ``<label>.txt`` to the output folder.
        **fit_options
            Other ``[fit]`` options for this fit, e.g. ``loss='soft_l1'``.

        Returns
        -------
        FitResult

        Examples
        --------
        Start simple, then let the velocity law break:

        >>> pl = wa.fit("power-law", n_boot=0)                       # doctest: +SKIP
        >>> bpl = wa.fit("broken-power-law", guide=pl, start_from=pl)  # doctest: +SKIP

        Test a free radial velocity, or hold the axis tilt fixed:

        >>> free_vr = wa.fit(radial="proportional", n_boot=0)        # doctest: +SKIP
        >>> no_tilt = wa.fit(fixed={"xslope": 0.0, "xoff_au": 0.0})  # doctest: +SKIP
        """
        from copy import copy
        fc = copy(self.config.fit)
        for k, v in fit_options.items():
            if not hasattr(fc, k):
                raise ValueError(f"unknown fit option {k!r}")
            setattr(fc, k, v)
        model = self.model(velocity_law, radial, rotation, fixed, bounds)
        guide, start_from = self._resolve(guide), self._resolve(start_from)
        n_boot = fc.n_boot if n_boot is None else n_boot
        iterations = fc.iterations if iterations is None else iterations
        if label is None:
            base = (model.law.name + ("" if model.radial == "cone" else f"-{model.radial}")
                    + ("" if model.rotation == "constant" else "-jz"))
            label, k = base, 2
            while label in self.results:          # never overwrite an earlier fit
                label, k = f"{base}-{k}", k + 1
        if traces is None:
            traces = self.extract_traces(y_offsets, guide=guide)
        if model.radial != "cone" and self.geometry.edge_on:
            self._log("edge-on view: Vz is unconstrained unless the flow is tied to the cone")
        n_it = max(1, int(iterations)) if guide is not None else 1
        res = None
        for it in range(n_it):
            last = it == n_it - 1
            fitter = GlobalFitter.from_config(model, self.geometry, traces, self.sig_x, fc,
                                              verbose=self.config.output.verbose)
            init_it = init if res is None else res.params
            res = fitter.fit(init=init_it, start_from=start_from, n_boot=n_boot if last else 0)
            if not last:
                self._log(f"guided iteration {it + 1}/{n_it}: re-extracting traces")
                traces = self.extract_traces(y_offsets, guide=res)
        res.label = label
        self.results[label] = res
        self.traces = traces
        if save:
            res.save(self.config.outpath(label, "json"))
            with open(self.config.outpath(label, "txt"), "w") as fh:
                fh.write(res.summary() + "\n")
        self._log(repr(res))
        return res

    def compare(self, *results, traces=None, show=True):
        """Compare fits on the same traces (default: those of the last one).

        Parameters
        ----------
        *results : FitResult or str
            Results or labels (default: all stored results).

        Returns
        -------
        list of dict
            ``label, k, n, chi2, bic, dbic`` sorted by BIC.
        """
        rs = [self._resolve(r) for r in results] or list(self.results.values())
        rows = compare_fits(rs, traces=traces)
        if show:
            print(f"{'model':28s} {'k':>3s} {'chi2':>10s} {'BIC':>10s} {'dBIC':>8s}")
            for r in rows:
                print(f"{r['label']:28s} {r['k']:3d} {r['chi2']:10.1f} {r['bic']:10.1f} "
                      f"{r['dbic']:8.1f}")
            print("(trace points oversample the beam: chi2/BIC differences are overestimated)")
        return rows

    def load_result(self, path_or_label):
        """Load a saved :class:`FitResult` (path, or label in the output folder)."""
        p = path_or_label if os.path.exists(str(path_or_label)) else \
            self.config.outpath(path_or_label, "json")
        r = FitResult.load(p)
        self.results[r.label] = r
        return r

    # ------------------------------------------------------------ others
    def ring_fits(self, y_offsets=None, prior_from=None, fix_from=None, guide=None):
        """Independent per-cut ring fits (:class:`RingFitter`), for full rings."""
        rs = RingFitter(self.cube).fit_all(y_offsets, prior_from=self._resolve(prior_from),
                                           fix_from=self._resolve(fix_from),
                                           guide=self._resolve(guide))
        rs.save_csv(self.config.outpath("ring_fits", "csv"))
        return rs

    def mass(self, result="last", plot=True, **mass_options):
        """Brightness, mass and mass-loss rate (:class:`MassAnalysis`)."""
        r = self._resolve(result) if result is not None else None
        m = MassAnalysis(self.cube, velocity_model=r, **mass_options).run()
        m.save(self.config.outpath("mass_slices", "txt"))
        if plot:
            m.plot(cube=self.cube, save=self.config.outpath("mass_slices"))
        self._log(m.summary().splitlines()[2])
        return m

    def derived(self, result="last", mass=None):
        """Derived quantities (:class:`~diskwind.derived.WindProperties`)."""
        return self._resolve(result).derived(mass=mass)

    # ------------------------------------------------------------ figures
    def _save(self, fig, name, save):
        if save:
            fig.savefig(self.config.outpath(name), dpi=self.config.output.dpi)
        return fig

    def plot_pv_grid(self, result=None, save=True, **kw):
        return self._save(plotting.plot_pv_grid(self.cube, result=self._resolve(result), **kw),
                          "pv_grid", save)

    def plot_gallery(self, result="last", save=True, **kw):
        return self._save(plotting.plot_gallery(self.cube, result=self._resolve(result), **kw),
                          "pv_gallery", save)

    def plot_slit_layout(self, save=True, **kw):
        return self._save(plotting.plot_slit_layout(self.cube, **kw), "slit_layout", save)

    def plot_cuts(self, result="last", others=(), save=True, **kw):
        r = self._resolve(result)
        fig = plotting.plot_cuts(self.cube, r, [self._resolve(o) for o in others], **kw)
        return self._save(fig, f"{r.label}_cuts", save)

    def plot_cut(self, y0, result="last", others=(), save=True, **kw):
        r = self._resolve(result)
        fig = plotting.plot_cut(self.cube, r, y0, [self._resolve(o) for o in others], **kw)
        return self._save(fig, f"{r.label}_cut_y{y0:.2f}", save)

    def plot_profiles(self, result="last", others=(), save=True, **kw):
        r = self._resolve(result)
        fig = plotting.plot_profiles(r, [self._resolve(o) for o in others], **kw)
        return self._save(fig, f"{r.label}_profiles", save)

    def plot_bootstrap(self, result="last", save=True):
        r = self._resolve(result)
        return self._save(plotting.plot_bootstrap(r), f"{r.label}_bootstrap", save)

    # --------------------------------------------------------------- run
    def run(self, first_law="power-law", n_boot=None, figures=True):
        """The standard analysis, driven by the configuration.

        1. noise and disk-rotation checks;
        2. an unguided fit with ``first_law`` (no bootstrap);
        3. the configured law (``fit.velocity_law``) guided by and started from it,
           ``fit.iterations`` times, with the bootstrap;
        4. model comparison, mass, derived quantities, figures.

        Returns
        -------
        dict
            ``first``, ``final``, ``comparison``, ``mass``, ``derived``.
        """
        cfg = self.config
        self.check_noise()
        self.check_disk_rotation()
        first = self.fit(first_law, n_boot=0, label=f"{first_law}-first")
        final = self.fit(guide=first, start_from=first, n_boot=n_boot)
        comp = self.compare(first, final)
        m = self.mass(final)
        props = self.derived(final, m)
        with open(cfg.outpath("derived", "txt"), "w") as fh:
            fh.write(props.summary() + "\n")
        print(props.summary())
        if figures:
            self.plot_cuts(final, others=[first])
            self.plot_profiles(final, others=[first])
            self.plot_gallery(final)
            self.plot_pv_grid(final)
            if final.boot is not None and len(final.boot):
                self.plot_bootstrap(final)
        self._log(f"done: results in {cfg.output_dir}")
        return dict(first=first, final=final, comparison=comp, mass=m, derived=props)

    def __repr__(self):
        return (f"WindAnalysis({self.config.source.name!r}, results={list(self.results)}, "
                f"cube={'loaded' if self._cube is not None else 'not loaded'})")
