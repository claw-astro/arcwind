Outputs
=======

Everything is written to ``output.dir`` (created if needed; a relative path is
relative to the configuration file). The examples keep inputs and results per source:
``examples/data/<source>/`` for the cubes and ``examples/output/<source>/`` for the
results (e.g. ``data/elias20/`` and ``output/elias20/``); see
``examples/data/README.md``.

``<label>.json``
   The full fit result: model, parameters, errors, bootstrap samples, traces,
   geometry. ``FitResult.load(path)`` or ``wa.load_result(label)`` restores it.
``<label>.txt``
   The text summary (``result.summary()``).
``noise_check.txt``, ``mass_slices.txt``, ``ring_fits.csv``, ``derived.txt``
   Tables of the corresponding steps.
Per-lobe outputs
   With two lobes, files of the counter lobe end in ``_counter``
   (``mass_slices_counter.txt``, ``<label>_cuts_counter.pdf`` ...), and a two-lobe
   result stores the traces of both lobes.
Figures
   ``<label>_cuts``, ``<label>_cut_y<y0>``, ``<label>_profiles``, ``<label>_angular_momentum``, ``<label>_rotation_per_cut`` (and ``.txt``),
   ``<label>_bootstrap``, ``pv_grid``, ``pv_gallery``, ``slit_layout``,
   ``mass_slices``, ``noise_check``, ``disk_rotation_check`` (format
   ``output.fig_format``).

Using a saved result
--------------------

.. code-block:: python

   from diskwind import FitResult
   r = FitResult.load("diskwind_output/broken-power-law.json")
   r.params, r.errors(), r.theta_deg
   r.profiles(z)                     # with bands from the bootstrap
   r.curve(1.0)                      # model curve of a cut
   r.table("latex")                  # rows for a LaTeX table
