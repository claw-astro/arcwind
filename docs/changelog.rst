Changelog
=========

Unreleased
----------

* Contributing guide (``CONTRIBUTING.md``), support statement (``SUPPORT.md``),
  governance (``GOVERNANCE.md``) and Code of Conduct, summarized in a new docs
  page, :doc:`contributing`; issue and pull-request templates.
* Acknowledgment of AI assistance, with the authors' confirmation of review.
* Python 3.10 or later is required (Python 3.8 and 3.9 have reached end of life and
  can no longer be installed on current CI machines). Tests run on 3.10 to 3.13.

0.3.0
-----

* The package is now called **arcwind** (Analysis of Resolved Conical Winds In
  Nascent Disks): package, import name, command and documentation. The code and
  results are unchanged. Default output folder: ``arcwind_output``.
* Logo (helix), in ``docs/_static/logo/``: README banner, docs sidebar logos for
  light and dark themes, favicon, app icons and the script that draws them.

0.2.7
-----

* A standard place for the data of each source: inputs in
  ``examples/data/<source>/`` and results in ``examples/output/<source>/``
  (``examples/data/README.md``). ``configs/elias20.toml`` now reads the Elias 2-20
  files from ``data/elias20/`` and writes to ``output/elias20/``; ``examples/output/``
  is ignored by git. A missing cube now gives an error that says where it was
  looked for.
* Example 02 names its adopted fit ``broken-power-law`` whatever the rotation law,
  and examples 08-10 also find a fit saved as ``broken-power-law-jz``.
* Code check: duplicate ``Trace.subset`` removed (one method takes boolean or index
  arrays); unused imports and variables removed from the examples.
* Elias 2-20 tutorial: the data and output folders; the paragraph of per-cut
  rotation values removed.

0.2.6
-----

New

* Rotation measured cut by cut (``wa.rotation_per_cut``, ``arcwind.rotation``):
  the geometry and velocity law of a global fit are held fixed and one j is fitted
  to each cut. Detections (``|j| >= 2`` sigma, well-fitted cuts only), the median j of
  all well-fitted cuts with a bootstrap range, one j fitted jointly to the detected
  cuts and an upper limit from the others; detections that are not adjacent are
  flagged (``<label>_rotation_per_cut.txt``). Figure: ``wa.plot_rotation``, with
  the median shown as a line and band.
* Example 02 fits a constant j to the same traces as the adopted j(z) fit
  (``CONSTANT_J``) and measures the rotation per cut; example 11 adds the
  rotation per cut to ``compare_ymax.txt``. Legends name a constant-j fit
  "broken power law, constant j".
* Fix: a fixed parameter is no longer clipped to the default bounds of the cuts
  being fitted.
* ``plot_angular_momentum`` / ``WindAnalysis.plot_angular_momentum``: j(z) and
  V_phi(z) with bootstrap bands (16-84% and 2.5-97.5%) and the Keplerian values at
  the cone radius, for scale (``<label>_angular_momentum``).
* Legends of the cut figures name the velocity law only ("power law", "broken
  power law"; a cut range in the label is kept, e.g. ``y0 <= 2.00"``).
* Trace points are drawn as black filled circles with error bars in the cut
  figures (``plot_cut``, ``plot_cuts``, ``plot_pv``).
* New example 11: the Elias 2-20 fit repeated with the cuts restricted to
  y0 <= YMAX (default 2.0" and 2.4"), each from scratch, compared parameter by
  parameter (``compare_ymax.txt``, ``compare_ymax_profiles``,
  ``compare_ymax_angular_momentum``), with j and V_phi at the lowest and highest
  cuts and the launch radius at the lowest cut.

Elias 2-20 analysis (example 02, ``configs/elias20.toml``)

* The spatial window follows the inner arch
  (``mask.x_cone = [0.19, 0.10, 0.50, -0.14]``); a wider window mixed in a fainter
  outer arc.
* Cuts from y = 0.45" to 2.40" every 0.15"; the top-cut overrides (1.95-2.40")
  no longer replace the window. Top-cut S/N thresholds 3.0, ring coverage
  threshold 0.25.
* j may change with height, j(z) = j (z/100 au)^qj (``fit.rotation = "power-law"``,
  adopted); a constant j is fitted to the same traces for comparison.
* Every fit starts from scratch (an unguided single power law, then the broken
  power law guided by it); ``YMAX`` restricts the cuts. The single power law used
  for the comparison is fitted to the same traces as the broken power law.
* One figure per fitted cut (``<label>_cut_y<y0>``) and the j(z) figure.
* Elias 2-20 tutorial: the two-shell structure, the solutions of the fit, j(z)
  versus a constant j, and the rotation measured cut by cut.

0.2.5
-----

* Faster first fits: multistart fits (the first, unguided fit and the j
  restarts) solve every starting point coarsely, then the ``fit.n_refine = 4``
  best to full precision, both steps in parallel on ``fit.n_jobs`` processes.
  The first Elias 2-20 power-law fit takes ~20 s instead of ~200 s (2 cores)
  and reaches a lower chi2.
* Slit layout: the maps are shown in mJy/beam km/s (moment 0) and mJy/beam
  (continuum), still converted from the units of each file (``BUNIT``).
  ``plot.map_min`` and ``plot.continuum_min`` are in these units again (Elias
  2-20 example: 0.064 and 0.11).

0.2.4
-----

* Faster fits and error estimates (bootstrap refits about 20x faster for Elias 2-20,
  before parallelism):

  - trace points are now spaced by half the beam along the cut
    (``trace.dx_arcsec`` auto; ``0`` restores every pixel). The cut selection
    (``trace.min_points``) still counts every pixel. Fewer, nearly independent
    points also make chi2 and BIC differences more honest;
  - the bootstrap refits run in parallel on all cores (``fit.n_jobs``; falls back
    to serial refits if the parallel run fails, e.g. on Windows);
  - the refits start at the best fit with a looser tolerance
    (``fit.resample_tol = 1e-6``), and the model rings use 181 points
    (``fit.nphi``);
  - new ``fit.error_method = "jackknife"``: one refit per left-out cut instead of
    ``n_boot`` bootstrap draws.

* Slit layout: images are read in their own units (``BUNIT``, e.g.
  'Jy/beam m/s') and shown in Jy/beam km/s (moment 0) and Jy/beam (continuum).
  Earlier versions ignored 'm/s' and showed moment-0 maps from such files 1000x
  too bright in mJy. ``plot.map_min`` and ``plot.continuum_min`` are now in
  these units (Elias 2-20 example: 6.4e-5 and 1.1e-4).

* Plot style: every figure uses matplotlib's 'classic' style with 14 pt
  sans-serif text and 2 pt axes (``arcwind.plotting.set_style``; the settings
  are ``arcwind.plotting._style.STYLE`` and ``RC``). Font sizes of labels,
  legends and titles scale with it; tick labels are thinned so they do not
  overlap. The schematic examples (08-10) use the same settings.
* Tutorials reordered: the synthetic wind, the HH 30-like wind, your own
  velocity law, per-cut ring fits, and the Elias 2-20 analysis last, as the
  real application. The ring-fit tutorial and example 05 say which data they
  use (the synthetic cube of the first tutorial).
* ``WindAnalysis`` plot methods accept a file path as ``save``.

0.2.3
-----

* HH 30-like mock: the docs now state that the real CO outflow of HH 30 is
  monopolar; the mock's counter lobe is a mirror image added for the example.
  The HH 30-like tutorial has a new section that fits a single lobe, like the
  real source (``MockWind.hh30_like(lobes="main")``).
* Ring schematic (example 10), panel (c): drawn relative to the fitted wind
  axis, with the cone extended below the disk to its apex and the axis offset
  and tilt listed below the panel.
* Elias 2-20 example configuration: trace S/N thresholds lowered by 1
  (``trace.snr_min`` 4.0 -> 3.0; per-cut overrides at y >= 2.00" 3.5 -> 2.5) so
  more points, and the faint top cuts, are kept.
* ``WindCube.center_radec`` now includes ``center_offset_arcsec`` (the offset
  was always applied to the cuts; only the reported center missed it).

0.2.2
-----

* Tutorials: many more figures made by the package (noise and disk-rotation
  checks, PV grid with masks, single cut, gallery, bootstrap, profiles, mass,
  ring fits, velocity laws, HH 30-like galleries and profiles, Elias 2-20
  profiles and schematics), each with a short explanation.
* ``plot_profiles``: readable labels, a visible flow-direction panel (minimum
  span, dashed cone angle), bands on R(z), and a note when the flow is tied to
  the cone.
* ``plot_bootstrap``: readable parameter labels and fewer ticks.
* Panel grids (``plot_pv_grid``, ``plot_gallery``, ``plot_cuts``): no
  overlapping tick labels between panels; the shared labels no longer touch
  the ticks for single-row grids. Saved figures are trimmed to their content.
* ``model_definitions`` schematic: fixed overlapping labels in panel (d).

0.2.1
-----

* Review pass: starting values of a free radial velocity are correct for any
  velocity law (not only the built-in ones); model comparison can include the
  counter lobe (``compare(..., other_traces=...)``); a warning when both lobes
  are analyzed with an automatic Vsys or disk sense.
* Docs: simpler limitations and citing; the bipolar-wind notes moved to the end of
  the workflow page.

0.2.0
-----

New

* **Bipolar winds.** ``source.lobes = "main" | "counter" | "both"``;
  ``wa.fit(lobe="both")`` fits one mirror-symmetric wind to the cuts of both
  lobes, each extracted in its own frame. ``FitResult.view("counter")`` gives a
  result in the other lobe's frame; every per-lobe step and figure takes
  ``lobe=``.
* **HH 30-like synthetic cube.** ``MockWind.hh30_like()`` (and
  ``arcwind mock --preset hh30``): nearly edge-on, both lobes, the parameters of
  Louvet et al. (2018). Example 07 and a tutorial.
* **Disk center from the continuum.** ``arcwind center continuum.fits`` and
  ``source.center_from_continuum``.
* **Slit-layout figure.** Continuum from ``source.continuum_fits``, cmasher color
  maps (``plot.map_cmap``), a line label (``plot.line_label``) and map floors.

Fixed

* Blank channels (all zero or NaN, e.g. at the edge of a cube) are no longer used
  for the noise; they made the rms too low.
* The Gaussian reference level in the noise-check figure (one-sided, 3.2e-5).
* An f-string that did not compile on Python < 3.12; Python 3.8 is supported.
* Separate fits of the counter lobe store their traces with that lobe.
* A bare model can guide the trace extraction; unknown fixed or bounded
  parameters and unknown result labels now raise clear errors.

Changed

* The schematic figures are numbered examples (08-10).
* American spelling throughout.

0.1.0
-----

* First release: object-oriented rewrite of the original Elias 2-20 scripts.
* Configuration files (TOML/YAML/JSON) with per-cut overrides and validation.
* Velocity-law registry (constant, power law, broken power law, user laws);
  radial modes cone / proportional / independent; constant or power-law j.
* ``WindAnalysis`` pipeline, ``arcwind`` command line, synthetic cubes.
