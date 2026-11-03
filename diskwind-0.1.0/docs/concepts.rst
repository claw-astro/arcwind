Concepts and conventions
========================

.. figure:: _static/model_definitions.png
   :width: 100%

   Definitions of the model (made with ``examples/paper_elias20/make_model_definitions.py``).
   (a) The cuts on the sky; (b) side view: a cut meets the near and far sides of
   the cone at different heights; (c) one ring and its velocity components;
   (d) the resulting curve in the PV diagram of a cut.

Frames
------

``X, Y``
   Sky offsets from the disk center. ``Y`` runs along the projected wind axis,
   towards the lobe (position angle ``wind_pa_deg``); ``X`` is perpendicular to it,
   towards ``wind_pa_deg - 90`` (the disk major axis).
``z, R, phi``
   Height along the wind axis, distance from the axis, and azimuth about ``+z``,
   counted counter-clockwise from ``+X`` seen from the tip of the lobe.
``i, s``
   The inclination of the axis to the line of sight (0-90 deg), and the lobe
   orientation: ``s = +1`` if the axis is tilted towards the observer (lobe
   blue-shifted on average), ``-1`` if away.
Disk rotation sense
   ``+1`` if the disk rotates counter-clockwise about ``+z`` (its ``+X`` side is
   red-shifted). A wind with ``j > 0`` then rotates with the disk.

The model
---------

The wind is a thin conical shell whose axis may be offset and tilted along X:

.. math::

   R(z) = R_0 + z\tan\theta, \qquad x_c(z) = x_0 + \xi z .

The gas has an axial velocity :math:`V_z(z)` (the :mod:`velocity law
<diskwind.laws>`), a radial velocity :math:`V_r(z)` and a rotation
:math:`V_\phi = j(z)/R`. A point of the ring at height :math:`z` and azimuth
:math:`\phi` is seen at

.. math::

   X = x_c + R\cos\phi, \qquad Y = z\sin i + sR\sin\phi\cos i,

with line-of-sight velocity

.. math::

   v - V_{\rm sys} = -sV_z\cos i + V_r\sin i\sin\phi + V_\phi\sin i\cos\phi .

Each ring is an ellipse on the sky, centered at :math:`Y = z\sin i`. A cut at
fixed :math:`Y = y_0` therefore meets the cone at heights that depend on
azimuth,

.. math::

   z(\phi) = \frac{y_0 - sR_0\sin\phi\cos i}{\sin i + s\tan\theta\sin\phi\cos i},

which requires :math:`\tan\theta < \tan i` (the fit allows up to
:math:`0.95\tan i`). Substituting :math:`z(\phi)` gives the **model curve** of
the cut in the PV diagram: an ellipse for an edge-on source (Louvet et al.
2018, their App. B), a distorted closed curve otherwise.

The fit
-------

In each PV diagram the locus of the emission (the **trace**) is measured by
Gaussian fits to the spectrum at every position (:mod:`diskwind.trace`). The
model is compared with the trace points of all cuts through their normalised
orthogonal distance to the model curve,

.. math::

   d_k^2 = \min_\phi\left[\left(\frac{x_k - X(\phi)}{\sigma_x}\right)^2 +
   \left(\frac{v_k - v(\phi)}{\sigma_k}\right)^2\right],

and :math:`\sum_k\rho(w_k^2d_k^2)` is minimised (:mod:`diskwind.fitting`).
Uncertainties come from the covariance matrix and a bootstrap over cuts.

Why tie the flow to the cone?
-----------------------------

When only half of each ring is detected, the top of each arch constrains only
:math:`V_z\cos i + V_r\sin i`; with :math:`V_r` free the fit can find
unphysical solutions (gas moving towards the axis, or sideways expansion with
:math:`V_z\approx0`). ``radial = "cone"`` (:math:`V_r = V_z\tan\theta`, as in
Louvet et al. 2018) removes this degeneracy. With full rings, ``radial =
"proportional"`` lets the data measure the flow direction.

Limitations
-----------

* The model fits the locus of the emission, not its intensity.
* :math:`R_0` is extrapolated below the lowest cut; it is not a launch radius.
* The wind axis is assumed to share the disk inclination.
* Trace points oversample the beam: chi2 and BIC differences are overestimated
  by roughly the number of points per beam.
* With partial rings, test the method on a synthetic cube with your geometry
  (:mod:`diskwind.mock`): biases can exceed the statistical errors.
