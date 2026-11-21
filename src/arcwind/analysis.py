"""
:class:`WindAnalysis` -- the high-level interface: one object that holds the
configuration, the cube, the traces and the fit results, with one method per
step of the analysis.

Typical session
---------------
>>> from arcwind import WindAnalysis
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
``arcwind run elias20.toml`` on the command line).
"""
from __future__ import annotations

import os

from . import plotting
from .config import WindConfig
from .constants import FWHM2SIG
from .cube import WindCube
from .fitting import FitResult, GlobalFitter, compare_fits
from .mass import MassAnalysis
from .model import ConicalWindModel
from .ringfit import RingFitter
from .trace import TraceExtractor

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
        self._cubes = {}
        if cube is not None:
            self._cubes[getattr(cube, "lobe", "main")] = cube
        self.results = {}
        self.traces = {}
        self.counter_traces = {}

    # ---------------------------------------------------------------- data
    @property
    def default_lobe(self):
        """The lobe analyzed by default: 'main', 'counter' or 'both' (``source.lobes``)."""
        return self.config.source.lobes

    def cube_for(self, lobe="main"):
        """The :class:`WindCube` in the frame of ``lobe`` ('main' or 'counter'),
        loaded on first access."""
        if lobe not in self._cubes:
            self._cubes[lobe] = WindCube(self.config, lobe=lobe)
        return self._cubes[lobe]

    @property
    def cube(self):
        """The cube of the default lobe (the main lobe when both are analyzed)."""
        return self.cube_for("counter" if self.default_lobe == "counter" else "main")

    @property
    def counter_cube(self):
        """The cube in the frame of the counter lobe."""
        return self.cube_for("counter")

    @property
    def geometry(self):
        return self.cube.geometry

    def _lobe(self, lobe):
        lobe = lobe or self.default_lobe
        if lobe not in ("main", "counter", "both"):
            raise ValueError("lobe must be 'main', 'counter' or 'both'")
        return lobe

    @staticmethod
    def _in_frame(r, lobe):
        """A result (or None) in the frame of a single lobe."""
        return r.view(lobe) if isinstance(r, FitResult) else r

    @property
    def sig_x(self):
        """Position uncertainty of trace points (arcsec)."""
        return self.config.trace.sigx_beam * self.cube.beam_x * FWHM2SIG

    def _log(self, msg):
        if self.config.output.verbose:
            print(f"[arcwind] {msg}")

    def pv(self, y0, **kw):
        """Transverse PV diagram (:meth:`WindCube.pv`)."""
        return self.cube.pv(y0, **kw)

    def check_noise(self, plot=True):
        """Noise-channel check (:func:`arcwind.diagnostics.check_noise`); saves
        ``noise_check.txt`` and the figure."""
        rep = self.cube.check_noise()
        with open(self.config.outpath("noise_check", "txt"), "w") as fh:
            fh.write(rep.summary() + "\n")
        if plot:
            rep.plot(save=self.config.outpath("noise_check"))
        self._log(rep.summary().splitlines()[-1])
        return rep

    def check_disk_rotation(self, y_offsets=(0.0, -0.15, -0.3), plot=True):
        """Disk rotation-sense check (:func:`arcwind.diagnostics.check_disk_rotation`)."""
        chk = self.cube.check_disk_rotation(y_offsets)
        if plot:
            chk.plot(save=self.config.outpath("disk_rotation_check"))
        self._log(chk.summary())
        return chk

    # -------------------------------------------------------------- traces
    def extract_traces(self, y_offsets=None, guide=None, lobe="main", **trace_options):
        """Extract the traces of the cuts of one lobe (:class:`TraceExtractor`).

        Parameters
        ----------
        y_offsets : array-like, optional
            Default: ``cuts.y_offsets_arcsec``.
        guide : FitResult, str or None
            Guided extraction around a fit (or the label of a stored result); a fit
            of the other lobe or of both lobes is mirrored into this lobe's frame.
        lobe : {'main', 'counter'}
            Which lobe (each in its own frame).
        **trace_options
            ``[trace]`` options for this extraction.

        Returns
        -------
        dict
            ``{y0: Trace}`` (also stored in ``self.traces`` or ``self.counter_traces``).
        """
        g = self._in_frame(self._resolve(guide), lobe)
        cube = self.cube_for(lobe)
        if isinstance(g, ConicalWindModel):            # a bare model: evaluate it here
            model, geom = g, cube.geometry
            g = lambda y0: model.curve(y0, geom)       # noqa: E731
            g.label = "a model"
        tr = TraceExtractor(cube).extract_all(y_offsets, guide=g, **trace_options)
        if lobe == "counter":
            self.counter_traces = tr
        else:
            self.traces = tr
        n = sum(t.n for t in tr.values())
        self._log(f"traces ({lobe} lobe): {len(tr)} cuts, {n} points"
                  + (f" (guided by {g.label})" if g is not None else ""))
        return tr

    def _resolve(self, r):
        """A result from a FitResult, a model, a stored label or 'last'."""
        if r is None or isinstance(r, (FitResult, ConicalWindModel)):
            return r
        if r in ("last", "previous"):
            if not self.results:
                raise ValueError("no fit has been made yet")
            return list(self.results.values())[-1]
        if r not in self.results:
            raise KeyError(f"no result labeled {r!r}; stored results: {list(self.results)} "
                           "(load saved ones with load_result)")
        return self.results[r]

    # ----------------------------------------------------------------- fit
    def model(self, velocity_law=None, radial=None, rotation=None, fixed=None, bounds=None,
              **values):
        """A :class:`ConicalWindModel` as set by ``[fit]``, with optional changes.

        ``fixed`` and ``bounds`` given here must name parameters of the model (an
        error says which exist); those in the configuration are applied only where
        they exist, since one configuration serves several models.
        """
        cfg = self.config
        if velocity_law or radial or rotation:
            cfg = cfg.copy()
            cfg.fit.velocity_law = velocity_law or cfg.fit.velocity_law
            cfg.fit.radial = radial or cfg.fit.radial
            cfg.fit.rotation = rotation or cfg.fit.rotation
        m = ConicalWindModel.from_config(cfg)
        m.fix(**(fixed or {}))
        m.set_bounds(**(bounds or {}))
        return m.set(**values)

    def fit(self, velocity_law=None, radial=None, rotation=None, fixed=None, bounds=None,
            guide=None, start_from=None, iterations=None, n_boot=None, y_offsets=None,
            init=None, label=None, traces=None, save=True, lobe=None, **fit_options):
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
            Use these traces instead of extracting them (for ``lobe='both'``: a pair
            ``(main_traces, counter_traces)``).
        lobe : {'main', 'counter', 'both'}, optional
            Which lobe to fit (default ``source.lobes``). ``'both'`` fits one
            mirror-symmetric wind to the cuts of both lobes; the parameters refer to
            the main lobe's frame (``res.view('counter')`` for the other).
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

        Both lobes of a bipolar wind, jointly or one at a time:

        >>> both = wa.fit(lobe="both")                                # doctest: +SKIP
        >>> north, south = wa.fit(lobe="main"), wa.fit(lobe="counter")  # doctest: +SKIP

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
        lobe = self._lobe(lobe)
        frame = "counter" if lobe == "counter" else "main"
        guide = self._resolve(guide)
        start_from = self._in_frame(self._resolve(start_from), frame)
        n_boot = fc.n_boot if n_boot is None else n_boot
        iterations = fc.iterations if iterations is None else iterations
        if label is None:
            base = (model.law.name + ("" if model.radial == "cone" else f"-{model.radial}")
                    + ("" if model.rotation == "constant" else "-jz")
                    + ("" if lobe == "main" else f"-{lobe}"))
            label, k = base, 2
            while label in self.results:          # never overwrite an earlier fit
                label, k = f"{base}-{k}", k + 1

        def extract(g):
            if lobe == "both":
                return (self.extract_traces(y_offsets, guide=g, lobe="main"),
                        self.extract_traces(y_offsets, guide=g, lobe="counter"))
            return self.extract_traces(y_offsets, guide=g, lobe=frame), None

        if traces is None:
            traces, ctraces = extract(guide)
        elif lobe == "both":
            traces, ctraces = traces
        else:
            ctraces = None
        geom = self.cube_for(frame).geometry
        if model.radial != "cone" and geom.edge_on:
            self._log("edge-on view: Vz is unconstrained unless the flow is tied to the cone")
        n_it = max(1, int(iterations)) if guide is not None else 1
        res = None
        for it in range(n_it):
            last = it == n_it - 1
            fitter = GlobalFitter.from_config(model, geom, traces, self.sig_x, fc,
                                              verbose=self.config.output.verbose,
                                              counter_traces=ctraces)
            init_it = init if res is None else res.params
            res = fitter.fit(init=init_it, start_from=start_from, n_boot=n_boot if last else 0)
            res.frame = frame
            if not last:
                self._log(f"guided iteration {it + 1}/{n_it}: re-extracting traces")
                traces, ctraces = extract(res)
        res.label = label
        self.results[label] = res
        if frame == "counter":
            self.counter_traces = traces
        else:
            self.traces = traces
            if ctraces:
                self.counter_traces = ctraces
        if save:
            res.save(self.config.outpath(label, "json"))
            with open(self.config.outpath(label, "txt"), "w") as fh:
                fh.write(res.summary() + "\n")
        self._log(repr(res))
        return res

    def compare(self, *results, traces=None, other_traces=None, show=True):
        """Compare fits on the same traces (default: those of the last one).

        Parameters
        ----------
        *results : FitResult or str
            Results or labels (default: all stored results).
        traces, other_traces : dict, optional
            Common traces of the main lobe (and of the counter lobe).

        Returns
        -------
        list of dict
            ``label, k, n, chi2, bic, dbic`` sorted by BIC.
        """
        rs = [self._resolve(r) for r in results] or list(self.results.values())
        rows = compare_fits(rs, traces=traces, other_traces=other_traces)
        if show:
            print(f"{'model':28s} {'k':>3s} {'chi2':>10s} {'BIC':>10s} {'dBIC':>8s}")
            for r in rows:
                print(f"{r['label']:28s} {r['k']:3d} {r['chi2']:10.1f} {r['bic']:10.1f} "
                      f"{r['dbic']:8.1f}")
            print("(trace points are ~2 per beam: chi2/BIC differences are overestimated "
                  "by about that factor)")
        return rows

    def load_result(self, path_or_label):
        """Load a saved :class:`FitResult` (path, or label in the output folder)."""
        p = path_or_label if os.path.exists(str(path_or_label)) else \
            self.config.outpath(path_or_label, "json")
        r = FitResult.load(p)
        self.results[r.label] = r
        return r

    # ------------------------------------------------------------ others
    def _single(self, lobe):
        """A single lobe for per-lobe steps: 'both' -> 'main'."""
        lobe = self._lobe(lobe)
        return "main" if lobe == "both" else lobe

    def _sfx(self, lobe):
        return "" if lobe == "main" else f"_{lobe}"

    def ring_fits(self, y_offsets=None, prior_from=None, fix_from=None, guide=None, lobe=None):
        """Independent per-cut ring fits (:class:`RingFitter`) of one lobe, for full
        rings (``lobe`` default: the main lobe)."""
        lobe = self._single(lobe)
        v = lambda r: self._in_frame(self._resolve(r), lobe)
        rs = RingFitter(self.cube_for(lobe)).fit_all(y_offsets, prior_from=v(prior_from),
                                                     fix_from=v(fix_from), guide=v(guide))
        rs.save_csv(self.config.outpath("ring_fits" + self._sfx(lobe), "csv"))
        return rs

    def mass(self, result="last", plot=True, lobe=None, **mass_options):
        """Brightness, mass and mass-loss rate of one lobe (:class:`MassAnalysis`).

        With two lobes, call it once per lobe (``lobe='main'``, ``lobe='counter'``)."""
        lobe = self._single(lobe)
        r = self._in_frame(self._resolve(result), lobe) if result is not None else None
        cube = self.cube_for(lobe)
        m = MassAnalysis(cube, velocity_model=r, **mass_options).run()
        m.save(self.config.outpath("mass_slices" + self._sfx(lobe), "txt"))
        if plot:
            m.plot(cube=cube, save=self.config.outpath("mass_slices" + self._sfx(lobe)))
        self._log(f"{lobe} lobe: " + m.summary().splitlines()[2])
        return m

    def derived(self, result="last", mass=None):
        """Derived quantities (:class:`~arcwind.derived.WindProperties`)."""
        return self._resolve(result).derived(mass=mass)

    # ------------------------------------------------------------ figures
    # Every figure method takes ``lobe`` ('main' or 'counter'; default: the main
    # lobe, or the counter lobe if that is the only one analyzed). A two-lobe fit
    # is shown in the frame of the requested lobe.
    def _save(self, fig, name, save):
        """``save``: True (to ``output.dir``), False, or a file path."""
        if save:
            path = save if isinstance(save, str) else self.config.outpath(name)
            fig.savefig(path, dpi=self.config.output.dpi, bbox_inches="tight", pad_inches=0.05)
        return fig

    def _r(self, result, lobe):
        return self._in_frame(self._resolve(result), lobe)

    def plot_pv_grid(self, result=None, save=True, lobe=None, **kw):
        lobe = self._single(lobe)
        fig = plotting.plot_pv_grid(self.cube_for(lobe), result=self._r(result, lobe), **kw)
        return self._save(fig, "pv_grid" + self._sfx(lobe), save)

    def plot_gallery(self, result="last", save=True, lobe=None, **kw):
        lobe = self._single(lobe)
        fig = plotting.plot_gallery(self.cube_for(lobe), result=self._r(result, lobe), **kw)
        return self._save(fig, "pv_gallery" + self._sfx(lobe), save)

    def plot_slit_layout(self, save=True, lobe=None, **kw):
        lobe = self._single(lobe)
        return self._save(plotting.plot_slit_layout(self.cube_for(lobe), **kw),
                          "slit_layout" + self._sfx(lobe), save)

    def plot_cuts(self, result="last", others=(), save=True, lobe=None, **kw):
        lobe = self._single(lobe)
        r = self._r(result, lobe)
        fig = plotting.plot_cuts(self.cube_for(lobe), r, [self._r(o, lobe) for o in others],
                                 **kw)
        return self._save(fig, f"{r.label}_cuts" + self._sfx(lobe), save)

    def plot_cut(self, y0, result="last", others=(), save=True, lobe=None, **kw):
        lobe = self._single(lobe)
        r = self._r(result, lobe)
        fig = plotting.plot_cut(self.cube_for(lobe), r, y0, [self._r(o, lobe) for o in others],
                                **kw)
        return self._save(fig, f"{r.label}_cut_y{y0:.2f}" + self._sfx(lobe), save)

    def plot_profiles(self, result="last", others=(), save=True, **kw):
        r = self._resolve(result)
        fig = plotting.plot_profiles(r, [self._resolve(o) for o in others], **kw)
        return self._save(fig, f"{r.label}_profiles", save)

    def plot_angular_momentum(self, result="last", others=(), save=True, **kw):
        """j(z) and V_phi(z) with bootstrap bands (:func:`plotting.plot_angular_momentum`)."""
        r = self._resolve(result)
        fig = plotting.plot_angular_momentum(r, [self._resolve(o) for o in others], **kw)
        return self._save(fig, f"{r.label}_angular_momentum", save)

    def rotation_per_cut(self, result="last", nsigma=2.0, save=True, **kw):
        """j measured in each cut with the geometry and velocity law of ``result``
        held fixed (:func:`arcwind.rotation.rotation_per_cut`); writes
        ``<label>_rotation_per_cut.txt``.

        >>> rot = wa.rotation_per_cut(bpl)          # doctest: +SKIP
        >>> print(rot.summary())                    # doctest: +SKIP
        """
        from .rotation import rotation_per_cut
        r = self._resolve(result)
        prof = rotation_per_cut(r, nsigma=nsigma, **kw)
        if save:
            prof.save(self.config.outpath(f"{r.label}_rotation_per_cut", "txt"))
        return prof

    def plot_rotation(self, profile, result=None, save=True, **kw):
        """Figure of a :class:`~arcwind.rotation.RotationProfile`
        (``<label>_rotation_per_cut``); ``result`` adds the global fit's j(z)."""
        r = self._resolve(result) if result is not None else None
        fig = plotting.plot_rotation(profile, r, **kw)
        return self._save(fig, f"{profile.label}_rotation_per_cut", save)

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
        if self.default_lobe == "both":
            self.mass(final, lobe="counter")
        props = self.derived(final, m)
        with open(cfg.outpath("derived", "txt"), "w") as fh:
            fh.write(props.summary() + "\n")
        print(props.summary())
        if figures:
            lobes = ["main", "counter"] if self.default_lobe == "both" else [None]
            for lb in lobes:
                self.plot_cuts(final, others=[first], lobe=lb)
                self.plot_gallery(final, lobe=lb)
                self.plot_pv_grid(final, lobe=lb)
            self.plot_profiles(final, others=[first])
            if final.boot is not None and len(final.boot):
                self.plot_bootstrap(final)
        self._log(f"done: results in {cfg.output_dir}")
        return dict(first=first, final=final, comparison=comp, mass=m, derived=props)

    def __repr__(self):
        return (f"WindAnalysis({self.config.source.name!r}, results={list(self.results)}, "
                f"lobes={self.default_lobe!r}, cubes loaded={sorted(self._cubes)})")
