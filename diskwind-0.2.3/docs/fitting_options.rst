Fitting options: what to turn on and off
========================================

Each option can be set in three places, from most permanent to most local:

1. the configuration file (``[fit]``, ``[trace]``, ``[mask]``...);
2. the configuration object: ``wa.config.fit.loss = "soft_l1"``;
3. the arguments of one call: ``wa.fit(loss="soft_l1")``.

The model
---------

.. list-table::
   :header-rows: 1
   :widths: 24 30 46

   * - Choice
     - How
     - Use it when
   * - Velocity law
     - ``fit.velocity_law`` / ``wa.fit("broken-power-law")``
     - ``"constant"``: one speed. ``"power-law"``: acceleration with one index.
       ``"broken-power-law"``: two regimes with a break height. Or
       :doc:`your own <tutorials/custom_law>`.
   * - Radial velocity
     - ``fit.radial`` / ``wa.fit(radial=...)``
     - ``"cone"`` (default): flow along the cone wall -- needed with half rings.
       ``"proportional"``: free Vr, same height dependence -- full rings.
       ``"independent"``: Vr with its own law (flow angle varies with height).
   * - Rotation
     - ``fit.rotation``
     - ``"constant"`` j (angular momentum conserved), or ``"power-law"`` to test
       whether j varies with height (check ``qj`` against mock tests).
   * - Fixed parameters
     - ``fit.fixed = {xslope = 0.0}`` / ``wa.fit(fixed={...})``
     - Hold e.g. the axis tilt ``xslope``, offset ``xoff_au`` or ``j`` fixed,
       e.g. to test how much rotation depends on them.
   * - Bounds
     - ``fit.bounds = {Vz = [0, 30]}`` / ``wa.fit(bounds={...})``
     - Exclude unphysical regions (defaults are wide).
   * - Starting point
     - ``wa.fit(init={...})``, ``start_from=result``
     - Start a complex law from a simpler fit (recommended).

Parameter names
   ``R0_au, tan_theta, xoff_au, xslope``; the law's (``Vz, qv`` or
   ``Vz_b, z_b, q1, q2``); radial (``Vr`` or ``Vr_b, r1, r2``); ``j, qj``.
   ``ConicalWindModel(...).names`` lists them for any choice.

Lobes
-----

.. list-table::
   :header-rows: 1
   :widths: 24 30 46

   * - Choice
     - How
     - Use it when
   * - One lobe
     - ``source.lobes = "main"`` (default) or ``"counter"``; ``wa.fit(lobe=...)``
     - Only one lobe is detected, or to test the lobes one at a time.
   * - Both lobes, one wind
     - ``source.lobes = "both"`` / ``wa.fit(lobe="both")``
     - A bipolar wind: one mirror-symmetric model fitted to the cuts of both lobes
       (:ref:`two-lobes`).
   * - A fit in the other lobe's frame
     - ``result.view("counter")``
     - Plotting the counter lobe, or comparing separate fits of the two lobes.

The statistic and the optimizer
-------------------------------

.. list-table::
   :header-rows: 1
   :widths: 24 30 46

   * - Option
     - Default
     - Effect
   * - ``fit.loss``
     - ``"linear"``
     - ``"soft_l1"``, ``"huber"``, ``"cauchy"``: down-weight stray trace points.
   * - ``fit.equal_cut_weight``
     - ``false``
     - Every cut counts equally (faint, sparsely traced cuts are not swamped).
   * - ``trace.sigx_beam``
     - ``0.5``
     - Position error of trace points, in beam dispersions along X.
   * - ``trace.err_inflate``
     - ``2``
     - Velocity errors = this x the formal centroid error (Louvet et al. 2018).
   * - ``fit.j_multistart``
     - ``true``
     - Restart with j = +-j_s, +-3 j_s (avoids the wrong rotation sign).
   * - ``fit.max_nfev``, ``fit.nphi``
     - ``3000``, ``361``
     - Optimizer effort; azimuths per model curve.

Traces
------

.. list-table::
   :header-rows: 1
   :widths: 24 30 46

   * - Option
     - Default
     - Effect
   * - ``trace.mode``
     - ``"wings"``
     - ``"single"`` when only one half of the rings is seen.
   * - ``trace.snr_min``
     - ``5``
     - Minimum peak S/N of a trace point.
   * - guided extraction
     - ``wa.fit(guide=result)``
     - Only channels within ``trace.guide_dv_kms`` of the model are used.
   * - ``trace.peak_choice``
     - ``"brightest"``
     - ``"nearest"``: in guided extraction, trace the peak closest to the model.
   * - ``fit.iterations``
     - ``1``
     - With a guide: re-extract around the new model and refit.
   * - ``trace.xsmooth_arcsec``, ``trace.vsmooth_chan``
     - off
     - Smooth faint cuts (the noise is smoothed identically).
   * - ``trace.drop_clipped``
     - ``true``
     - Drop points whose line is cut by a mask edge (set ``false`` if the arch
       legs run into a mask edge and are needed).

Uncertainties
-------------

.. list-table::
   :header-rows: 1
   :widths: 24 30 46

   * - Option
     - Default
     - Effect
   * - ``fit.n_boot``
     - ``100``
     - Bootstrap over cuts; ``0`` for fast exploratory fits.
   * - ``fit.scale_errors_by_chi2``
     - ``true``
     - Covariance errors x sqrt(chi2_red) when > 1.
   * - ``source.vsys_err_kms``
     - ``0.1``
     - Added to the velocity amplitude.
   * - ``result.errors(kind)``
     - ``"max"``
     - ``"cov"``, ``"boot"`` or the larger of the two.

Recipes
-------

Fast exploration:

.. code-block:: python

   wa = WindAnalysis("src.toml", n_boot=0)

Robustness checks for a paper:

.. code-block:: python

   base = wa.fit(...)
   tests = {
       "no tilt":  wa.fit(fixed={"xslope": 0.0, "xoff_au": 0.0}, guide=base),
       "robust":   wa.fit(loss="soft_l1", guide=base),
       "j(z)":     wa.fit(rotation="power-law", guide=base),
       "q2 = 1":   wa.fit("broken-power-law", fixed={"q2": 1.0}, guide=base, start_from=base),
   }
   wa.compare(base, *tests.values(), traces=base.traces)

Test your settings on a synthetic cube with your geometry:

.. code-block:: python

   from diskwind.mock import MockWind, make_mock_cube
   m = MockWind(inclination_deg=49, approaching=False, theta_deg=20, vp=10, j=-40)
   make_mock_cube("mock.fits", m)
   cfg = wa.config.copy(); cfg.source.fits_file = "mock.fits"
