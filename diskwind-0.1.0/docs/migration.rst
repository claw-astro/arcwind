Migrating from the original scripts
===================================

The original analysis was a set of scripts sharing ``windkit.py`` and a Python
configuration with upper-case keys. diskwind reads those configurations
directly and can convert them:

.. code-block:: bash

   diskwind convert elias20_config.py        # -> elias20_config.toml

or ``WindConfig.load("elias20_config.py")``. The mapping of every key is in
``diskwind.config.LEGACY_KEYS``; for example ``TRACE_SNR_MIN`` ->
``trace.snr_min``, ``TRACE_X_CONE`` -> ``mask.x_cone``,
``TRACE_SLIT_OVERRIDES`` -> ``cuts.overrides``. ``VELOCITY_MODEL`` and
``GLOBAL_FLOW_ALONG_CONE`` become ``fit.velocity_law`` and ``fit.radial``.

Scripts and their replacements
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Script
     - diskwind
   * - ``fit_wind_global.py``
     - ``wa.fit("power-law")`` (``--guided``: ``guide=``; ``--iterations``)
   * - ``fit_wind_global_bpl.py --mode cone``
     - ``wa.fit("broken-power-law", guide=pl, start_from=pl)``
   * - ``fit_wind_global_bpl.py --mode separate``
     - ``wa.fit("broken-power-law", radial="independent", start_from=bpl)``
   * - ``fit_wind_rings.py``
     - ``wa.ring_fits()`` (``--prior-from-global``: ``prior_from=``)
   * - ``make_fig10_mass.py``
     - ``wa.mass(result)``
   * - ``compare_velocity_models.py``, ``compare_hh30.py``
     - ``wa.compare(...)``, ``result.derived(mass).summary()``,
       ``.compare_with("HH 30")``
   * - ``check_noise.py``
     - ``wa.check_noise()``
   * - ``make_disk_rotation_check.py``
     - ``wa.check_disk_rotation()``
   * - ``make_fig7_pv_grid.py``
     - ``wa.plot_pv_grid()``, ``wa.plot_slit_layout()``
   * - ``make_fig7_slit_layout_fancy.py``
     - ``wa.plot_slit_layout(mom0=..., continuum=...)``
   * - ``make_fig8_gallery.py``
     - ``wa.plot_gallery()``
   * - ``make_global_profiles.py``
     - ``wa.plot_profiles()``
   * - ``make_mock_cube.py``
     - ``diskwind.mock.make_mock_cube`` / ``diskwind mock``
   * - schematic figures
     - ``examples/paper_elias20/``

Numerical agreement
-------------------

On the synthetic cube of the original validation, diskwind reproduces the
original single-power-law and broken-power-law fits exactly (same traces, same
parameters, chi2 and BIC) and the same masses and mass-loss rates.
