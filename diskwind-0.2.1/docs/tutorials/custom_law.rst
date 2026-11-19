Tutorial: your own velocity law
===============================

Velocity laws are small classes registered by name
(``examples/04_custom_velocity_law.py``):

.. code-block:: python

   import numpy as np
   from diskwind.laws import VelocityLaw, register_law

   @register_law
   class Saturating(VelocityLaw):
       """V(z) = V_inf (1 - exp(-z / z_acc))."""
       name = "saturating"
       param_names = ("V_inf", "z_acc")
       amplitude = "V_inf"                 # scales the velocity (Vsys error goes here)
       radial_map = {"V_inf": "Vr_inf"}     # names of an independent radial copy

       def evaluate(self, z, p):
           return p["V_inf"] * (1.0 - np.exp(-z / p["z_acc"]))

       def default_bounds(self, z_cuts):
           return {"V_inf": (0.0, 50.0), "z_acc": (1.0, 10 * np.max(z_cuts))}

       def guess(self, v_ref, z_cuts):
           return {"V_inf": 1.2 * v_ref, "z_acc": 0.5 * float(np.median(z_cuts))}

   res = wa.fit("saturating")

or ``velocity_law = "saturating"`` in the configuration (the module defining the
law must be imported first).

What each part does
-------------------

``evaluate(z, p)``
   The velocity at heights ``z`` (au); ``p`` maps parameter names to values
   (scalars, or arrays when evaluating bootstrap draws -- write it with numpy
   broadcasting).
``default_bounds(z_cuts)``
   Bounds given the central heights of the cuts.
``guess(v_ref, z_cuts)``
   Starting values from a typical speed of the traces.
``starts_from_profile(vfunc, z_cuts)`` (optional)
   Starting points that follow a previous fit's profile, used by
   ``wa.fit(..., start_from=previous)``.
``latex()`` (optional)
   The formula, shown by ``model.describe()``.
