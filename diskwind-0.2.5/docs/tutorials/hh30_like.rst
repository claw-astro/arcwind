Tutorial: an HH 30-like bipolar wind
====================================

This tutorial fits both lobes of a nearly edge-on bipolar wind, using a synthetic
cube with the parameters that Louvet et al. (2018) measured for HH 30
(``examples/07_hh30_like_both_lobes.py``; about 5 minutes, no data needed).

The real CO outflow of HH 30 is monopolar: it is seen only above the north-east
face of the disk. The counter lobe of the mock is added, as a mirror image of the
main lobe, to show how two lobes are fitted.

The HH 30-like preset
---------------------

:meth:`MockWind.hh30_like() <diskwind.mock.MockWind.hh30_like>` builds the
disk and wind below. Any of these values can be changed with keyword arguments,
e.g. ``MockWind.hh30_like(j=0.0)``.

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Quantity
     - HH 30-like mock
     - Louvet et al. (2018)
   * - distance, stellar mass, Vsys
     - 140 pc, 0.45 Msun, 6.9 km/s
     - same
   * - wind axis
     - PA 31.2 deg, i = 89 deg (main lobe tilted away)
     - PA 31.2 deg, i = 91 deg
   * - cone half-opening angle, R0
     - 35 deg, 20 au
     - ~35 deg, R0 <= 22 au
   * - poloidal speed
     - 9.3 km/s, constant
     - 9.3 +/- 0.7 km/s, constant
   * - specific angular momentum
     - 38 au km/s, with the disk
     - 38 +/- 15 au km/s, with the disk
   * - lobes
     - both, mirror-symmetric
     - one: the CO outflow is monopolar (north-east face only)
   * - beam, channels
     - 0.26" x 0.19" (PA 29 deg), 0.3 km/s
     - same
   * - CO disk
     - Keplerian, r <= 150 au, z/r <= 0.2
     - r ~ 145 au (continuum)

An inclination of 91 deg is the same view as 89 deg with the lobe tilted the
other way; diskwind uses 0-90 deg and the lobe orientation
(``wind_approaching``). Close to edge-on, the tilt and the disk rotation sense
cannot be estimated reliably from the cube, so the mock configuration states
them, as one would from the jet or scattered-light images.

.. code-block:: python

   from diskwind import WindAnalysis
   from diskwind.mock import MockWind, make_mock_cube

   truth = MockWind.hh30_like()
   make_mock_cube("hh30_mock.fits", truth)
   cfg = truth.config("hh30_mock.fits")     # source.lobes = "both"
   wa = WindAnalysis(cfg)
   wa.plot_slit_layout()                    # cuts on both sides of the disk

.. figure:: ../_static/hh30_slit_layout.png
   :width: 70%

   The HH 30-like mock: integrated intensity, with the cuts of both lobes.

One lobe, like the real source
------------------------------

To match the real, monopolar CO outflow of HH 30, make the mock with the main
(north-east) lobe only and fit it as a single-lobe wind:

.. code-block:: python

   truth = MockWind.hh30_like(lobes="main")
   make_mock_cube("hh30_one_lobe.fits", truth)
   cfg = truth.config("hh30_one_lobe.fits")   # source.lobes = "main"
   wa = WindAnalysis(cfg)
   wa.plot_slit_layout()

   res = wa.fit("constant", n_boot=30)
   print(res.summary())
   wa.plot_gallery(res)

.. figure:: ../_static/hh30_one_lobe_slit_layout.png
   :width: 60%

   The one-lobe mock: the wind leaves only the north-east face of the disk, and
   the cuts are placed on that side only.

The fit recovers theta = 34.8 +/- 0.2 deg (true 35), a poloidal speed
Vz / cos(theta) = 9.3 km/s (9.3), j = +37.6 +/- 0.6 au km/s (38; the same
sense as the disk) and R0 = 19.4 +/- 1.4 au (20).

.. figure:: ../_static/hh30_one_lobe_gallery.png
   :width: 80%

   The one-lobe fit. At this inclination both halves of each ring are detected
   in a single lobe, so one lobe already constrains the wind well; the second
   lobe mainly tests whether the wind is symmetric.

For a real nearly edge-on source with one lobe, set ``vsys_kms``,
``wind_approaching`` and ``disk_rotation_sense`` in the configuration (the
mock does this for you), for the reasons given above.

Fitting both lobes
------------------

With ``lobe="both"`` (or ``source.lobes = "both"`` in the configuration), the
cuts of both lobes are extracted, each lobe in its own frame, and **one**
mirror-symmetric wind is fitted to all of them:

.. code-block:: python

   both = wa.fit("constant", lobe="both")
   print(both.summary())

On the mock (28 cuts, 14 per lobe) this recovers theta = 35.1 deg (true 35),
Vp = 9.25 km/s (true 9.3) and j = 38.6 au km/s (true 38). The parameters refer to the main
lobe's frame. ``both.view("counter")`` gives the same result in the counter
lobe's frame, which is what the figures of the counter lobe use:

.. code-block:: python

   wa.plot_cuts(both, lobe="main")
   wa.plot_cuts(both, lobe="counter")

.. figure:: ../_static/hh30_counter_cuts.png
   :width: 100%

   The counter lobe, in its own frame, with the two-lobe model (red) and the
   traces (black).

.. code-block:: python

   wa.plot_gallery(both, lobe="main")
   wa.plot_gallery(both, lobe="counter")

.. figure:: ../_static/hh30_gallery_main.png
   :width: 80%

   Main lobe. Seen edge-on, both halves of each ring are detected and the ring
   is centered on Vsys (cross): Vz has no component along the line of sight.
   The ring's velocity extent gives Vr and its width in X gives R(z).

.. figure:: ../_static/hh30_gallery_counter.png
   :width: 80%

   Counter lobe, in its own frame, with the same model mirrored. In the lowest cut
   the disk (the tilted band through Vsys) crosses the ring; it is masked out of
   the traces.

Is the speed constant? A power law fitted to the same cuts gives
``qv = 0.003``: the speed is constant, as in the mock. Its formal BIC is lower
nonetheless (Delta BIC ~ 12), because neighboring trace points (two per beam) are
not independent -- judge such comparisons by the size of the change, not by the
BIC alone:

.. code-block:: python

   pl = wa.fit("power-law", lobe="both", start_from=both)
   wa.compare(both, pl)
   wa.plot_profiles(both, others=[pl])

.. figure:: ../_static/hh30_profiles.png
   :width: 100%

   The constant-speed fit (solid, with bootstrap bands) and the power law
   (dashed) agree to within ~0.2 km/s at all heights. Edge-on, the
   ring-center velocity (second panel) is ~0, and the rotation velocity falls
   as 1/R (fourth panel, blue), as for a constant specific angular momentum.

Each lobe on its own
--------------------

Fitting the lobes separately tests whether the wind is really symmetric:

.. code-block:: python

   north = wa.fit("constant", lobe="main")
   south = wa.fit("constant", lobe="counter")
   south.view("main").params        # in the main lobe's sign conventions

On the mock, the main lobe gives theta = 34.9 deg, Vz = 7.62 km/s and
j = +40 au km/s, the counter lobe 35.3 deg, 7.51 km/s and +38 au km/s (true:
35 deg, 7.62 km/s, +38 au km/s).

In its own frame the counter lobe has ``xoff_au`` and ``j`` of the opposite sign
(:ref:`two-lobes`); ``view("main")`` converts them, so the two sets of
parameters can be compared directly.

Per-cut ring fits
-----------------

Because both halves of every ring are detected, each cut can also be fitted on
its own, as in Louvet et al. (2018); see :doc:`ring_fits`:

.. code-block:: python

   rings = wa.ring_fits(lobe="main")
   rings.cone                        # R0 and theta from R(z)

On the mock the ring radii give R0 = 18.7 +/- 1.6 au and theta = 34.9 +/- 0.2 deg.

Mass
----

Masses are computed per lobe:

.. code-block:: python

   m_main = wa.mass(both, lobe="main")
   m_counter = wa.mass(both, lobe="counter")

.. figure:: ../_static/hh30_mass.png
   :width: 85%

   ``mass_slices.png`` for the main lobe (the counter lobe is written to
   ``mass_slices_counter.png``).
