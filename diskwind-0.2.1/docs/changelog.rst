Changelog
=========

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
  ``diskwind mock --preset hh30``): nearly edge-on, both lobes, the parameters of
  Louvet et al. (2018). Example 07 and a tutorial.
* **Disk center from the continuum.** ``diskwind center continuum.fits`` and
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
* ``WindAnalysis`` pipeline, ``diskwind`` command line, synthetic cubes.
