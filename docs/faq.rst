Troubleshooting and FAQ
=======================

The fit puts the cone on the wrong side / the curves look upside down
   Check ``source.wind_approaching`` (the lobe orientation ``s``) and
   ``source.wind_pa_deg``. A red-shifted lobe usually means ``false``.

The rotation sign is not robust
   Rotation comes from the velocity asymmetry between the two legs of each
   arch. Refit with the axis offset and tilt fixed
   (``fixed={"xslope": 0, "xoff_au": 0}``), look at the fraction of bootstrap
   samples keeping the sign (in ``summary()``), and test on a mock cube.

j(z) gives huge rotation somewhere, or qj sits at its bound
   A power-law j(z) fitted to a wind that rotates measurably only in a few cuts
   is extrapolated outside them. Compare with ``rotation="constant"`` on the same
   traces, and measure the rotation cut by cut (``wa.rotation_per_cut(result)``):
   quote the median, the detections and the upper limit, and keep j(z) as the
   fitting form for the geometry and velocities.

Vr comes out negative
   With half rings, "far half, Vr > 0" and "near half, Vr < 0" can look alike.
   ``fit.check_vr_sign`` refits with Vr >= 0 and reports both chi2; use
   ``radial = "cone"`` to tie the flow to the cone.

The model misses the fast emission at large heights
   Try ``"broken-power-law"`` started from the power-law fit; for faint cuts use
   per-cut overrides (wider slits, smoothing, lower ``snr_min``, an asymmetric
   ``guide_dv_kms``) and ``peak_choice = "nearest"``.

Some cuts are skipped
   They have fewer than ``trace.min_points`` points. Look at them with
   ``wa.pv(y0).plot(show_mask=True)``; widen the masks or lower ``snr_min``.

Only a few trace points per cut, although the rings are bright
   If the rings are much taller in velocity than they are wide (a narrow, fast
   cone), most of each ring is a steep wall. The default ``trace.mode = "wings"``
   fits one velocity per position, and along a wall that profile is broad or cut
   by a mask, so the point is dropped. Use ``trace.mode = "ring"``: the walls are
   then traced channel by channel (open squares in ``plot_cuts``). If the walls
   are patchy, check that ``mask.x_cone`` does not cut the outer edge of the
   walls (points whose profile reaches a mask edge are dropped), and raise
   ``trace.wall_max_fwhm_beam`` (default 4) if the walls are broad.

The noise check says CHECK
   A noise channel contains emission or one end of the cube is noisier: set
   ``noise.vel_ranges_kms`` to line-free channels. Blank channels (all zero or
   NaN) are excluded automatically; they would make the noise look too low.

How do I change the look of the figures?
   Every figure calls ``arcwind.plotting.set_style()`` (matplotlib's 'classic'
   style, 14 pt text, 2 pt axes). Edit ``STYLE`` and ``RC`` in
   ``arcwind/plotting/_style.py``, or restyle the returned figures yourself.

The slit-layout map is not in the cmasher colors
   Install cmasher (``pip install cmasher``); without it ``plot.map_cmap =
   "cmr.freeze"`` falls back to 'bone'.

chi2 / BIC differences are huge
   Trace points are spaced by half a beam (``trace.dx_arcsec``), so neighboring
   points are not independent and chi2 and BIC differences are overestimated by
   about a factor of 2. Judge models by the cuts (``plot_cuts``) and the
   stability of the parameters.

How long does a fit take?
   A few seconds for a guided fit or one started from an earlier fit; ~10-20 s
   for a first, unguided fit (several starting points, solved in parallel).
   The bootstrap refits ``n_boot`` times,
   starting from the best fit and in parallel on all cores (about 1 s per refit
   per core for Elias 2-20). ``fit.error_method = "jackknife"`` needs only one
   refit per cut. Use ``n_boot = 0`` while exploring.
