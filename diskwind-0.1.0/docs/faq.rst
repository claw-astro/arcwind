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

The noise check says CHECK
   A noise channel contains emission or one end of the cube is noisier: set
   ``noise.vel_ranges_kms`` to line-free channels.

chi2 / BIC differences are huge
   Trace points are sampled every pixel and oversample the beam, so chi2 and BIC
   differences are overestimated by about the number of points per beam. Judge
   models by the cuts (``plot_cuts``) and the stability of the parameters.

How long does a fit take?
   Seconds to a minute per fit for ~15 cuts; the bootstrap repeats the fit
   ``n_boot`` times. Use ``n_boot = 0`` while exploring.
