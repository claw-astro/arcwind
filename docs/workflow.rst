The analysis step by step
=========================

Every step is a method of :class:`~arcwind.WindAnalysis`; each also exists as a
class you can use on its own (:doc:`api/index`).

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Step
     - Method
     - Output (in ``output.dir``)
   * - Load and resample the cube
     - ``wa.cube`` (:class:`~arcwind.WindCube`)
     - --
   * - Noise channels line-free?
     - :meth:`~arcwind.WindAnalysis.check_noise`
     - ``noise_check.txt/.pdf``
   * - Disk rotation sense
     - :meth:`~arcwind.WindAnalysis.check_disk_rotation`
     - ``disk_rotation_check.pdf``
   * - Look at the PV diagrams and masks
     - ``plot_pv_grid(show_mask=True)``, ``plot_slit_layout()``
     - ``pv_grid.pdf``, ``slit_layout.pdf``
   * - Traces
     - :meth:`~arcwind.WindAnalysis.extract_traces`
     - (stored in the results)
   * - Global fit
     - :meth:`~arcwind.WindAnalysis.fit`
     - ``<label>.json``, ``<label>.txt``
   * - Model comparison
     - :meth:`~arcwind.WindAnalysis.compare`
     - printed table
   * - Per-cut ring fits (full rings)
     - :meth:`~arcwind.WindAnalysis.ring_fits`
     - ``ring_fits.csv``
   * - Mass, mass-loss rate
     - :meth:`~arcwind.WindAnalysis.mass`
     - ``mass_slices.txt/.pdf``
   * - Derived quantities
     - :meth:`~arcwind.WindAnalysis.derived`
     - (``WindProperties``)
   * - Figures
     - ``plot_cuts``, ``plot_cut``, ``plot_profiles``, ``plot_angular_momentum``, ``plot_rotation``, ``plot_gallery``, ``plot_bootstrap``
     - ``<label>_cuts.pdf`` ...

1. Set up and check
-------------------

.. code-block:: python

   from arcwind import WindAnalysis
   wa = WindAnalysis("mysource.toml")
   print(wa.cube)

The cube is resampled onto the wind frame, centered on the disk center (from the
configuration, or fitted to the continuum: see :doc:`configuration`). Unset values are estimated and printed:
Vsys from the line wings of the disk half without wind, the disk rotation sense,
and the lobe orientation. **Put the values you trust in the configuration** --
especially ``vsys_kms``, ``disk_rotation_sense`` and ``wind_approaching``.

``check_noise()`` compares the rms of the noise channels at both ends of the cube
cut by cut. Blank channels (all zero or NaN, as at the edges of some cubes) are
never used. If the line reaches one end, set ``noise.vel_ranges_kms`` to
line-free ranges.

2. Masks and traces
-------------------

Look at ``wa.plot_pv_grid(show_mask=True)``: cyan outlines are the pixels the
trace may not use. Adjust ``[mask]`` until the wind emission is inside and the
disk, cloud and unrelated emission outside:

* ``vrel_range_kms``, ``v_exclude_kms`` -- velocity windows;
* ``x_cone`` -- a spatial window that widens with height;
* ``kepler`` -- the Keplerian disk mask.

Choose ``trace.mode = "single"`` if only one half of each ring is seen, and
``trace.mode = "ring"`` if the rings are much taller in velocity than they are
wide (a narrow, fast cone): the default ``"wings"`` mode traces one velocity per
position and misses the steep sides of such rings, which ``"ring"`` traces
channel by channel (open squares in the cut plots). Faint
cuts can get their own settings with ``[cuts.overrides."y0"]``.

3. Fit, from simple to complex
------------------------------

.. code-block:: python

   first = wa.fit("power-law", n_boot=0)             # unguided, fast
   final = wa.fit("broken-power-law", guide=first,   # traces around the first model
                  start_from=first, iterations=2, n_boot=100)
   wa.compare(first, final)

Guided extraction keeps only emission near the previous model, which removes
contamination picked up by the first, unguided traces. ``iterations`` repeats
"re-extract around the new model, refit". See :doc:`fitting_options` for
everything that can be changed.

4. Check the fit
----------------

* ``wa.plot_cuts(final)``: the model must follow the trace in every cut;
* ``final.per_cut_d2``: mean squared distance per cut (large = poor fit there);
* ``wa.plot_bootstrap(final)``: are the parameters stable?
* compare models on the same traces with ``wa.compare``.

5. Physical quantities
----------------------

.. code-block:: python

   mass = wa.mass(final)               # slices along the axis, Mdot
   props = wa.derived(final, mass)
   props.at([100, 400])                # Vz, Vr, Vp, flow angle, dlnVp/dlnz
   props.launch_radius()               # j^2 / G M*
   props.travel_time()
   print(props.compare_with("HH 30"))

6. Everything at once
---------------------

``wa.run()`` or ``arcwind run mysource.toml`` performs the checks, a first
power-law fit, the configured law guided by it, the comparison, the mass and the
figures.

Bipolar winds
-------------

With ``source.lobes = "both"`` every step works per lobe: ``wa.counter_cube``,
``wa.extract_traces(lobe="counter")``, ``wa.mass(result, lobe="counter")``,
``wa.plot_cuts(result, lobe="counter")``, while ``wa.fit()`` fits one wind to
both lobes. The automatic Vsys and disk-sense estimates use the side of the disk
without wind, so for a bipolar wind set ``vsys_kms`` and ``disk_rotation_sense``
in the configuration. See :doc:`tutorials/hh30_like`.
