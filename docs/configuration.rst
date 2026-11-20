Configuration
=============

A configuration is a :class:`~arcwind.WindConfig`: sections of options with
defaults. Files can be TOML (recommended), YAML, JSON or Python.

.. code-block:: bash

   arcwind init mysource.toml      # annotated template
   arcwind options                 # every option, default and meaning
   arcwind options mask            # one section

In Python:

.. code-block:: python

   from arcwind import WindConfig
   cfg = WindConfig.load("mysource.toml")
   cfg.trace.snr_min = 4                    # attribute access
   cfg.set(snr_min=4, **{"fit.n_boot": 0})  # bare or dotted names
   cfg.save("mysource_v2.toml")
   print(WindConfig.describe("fit"))

Rules
-----

* ``"auto"`` means unset: arcwind decides (TOML has no null).
* Option names are checked: a typo raises an error with a suggestion.
* Relative paths (``fits_file``, ``output.dir``) are relative to the file.
* ``cuts.y_offsets_arcsec`` is a list or ``{start, stop, step}``.

The disk center
---------------

The origin of the wind frame is the disk center, best measured on the continuum.
Three ways to set it in ``[source]``:

* ``center_ra``/``center_dec`` (absolute), plus an optional
  ``center_offset_arcsec``;
* ``center_offset_arcsec`` alone: an offset from the reference pixel of the cube,
  e.g. a published center relative to the phase center (Elias 2-20 uses the
  ellipse fit of Huang et al. 2018);
* ``center_from_continuum = true`` with ``continuum_fits``: an elliptical Gaussian
  is fitted to the continuum (keep ``center_offset_arcsec = [0, 0]``).

``arcwind center continuum.fits`` runs the same fit and prints the lines to paste.

Per-cut overrides
-----------------

Any option can differ for one cut, e.g. for the faint top of a wind:

.. code-block:: toml

   [cuts.overrides."2.00"]
   x_range_arcsec = [-1.8, 1.8]      # replaces the x_cone window for this cut
   vrel_range_kms = [4.8, 12.1]
   slit_width_arcsec = 0.10
   snr_min = 3.5
   xsmooth_arcsec = 0.10
   vsmooth_chan = 2
   guide_dv_kms = [1.0, 2.5]         # asymmetric guide window

``cfg.for_cut(2.0)`` returns the configuration that cut sees.

Python configuration files
--------------------------

A ``.py`` file may define ``CONFIG = {"source": {...}, ...}``; Python is the
only format that can hold ``mask.disk_height_func`` (a function
``r_au -> H_au``).

Reference
---------

.. include:: config_reference.rst
