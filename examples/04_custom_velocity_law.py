"""
04 -- Adding your own velocity law.

A law needs a name, parameter names, the amplitude parameter, the names of an
independent radial copy, and three methods. It is then available everywhere
(configuration files included) under its name.
"""
import numpy as np

from arcwind import WindAnalysis
from arcwind.laws import VelocityLaw, register_law


@register_law
class Saturating(VelocityLaw):
    """V(z) = V_inf (1 - exp(-z / z_acc)): accelerates, then saturates."""
    name = "saturating"
    param_names = ("V_inf", "z_acc")
    amplitude = "V_inf"
    radial_map = {"V_inf": "Vr_inf"}          # z_acc shared if radial = 'independent'

    def evaluate(self, z, p):
        return p["V_inf"] * (1.0 - np.exp(-z / p["z_acc"]))

    def default_bounds(self, z_cuts):
        return {"V_inf": (0.0, 50.0), "z_acc": (1.0, 10 * np.max(z_cuts))}

    def guess(self, v_ref, z_cuts):
        return {"V_inf": 1.2 * v_ref, "z_acc": 0.5 * float(np.median(z_cuts))}

    def latex(self):
        return r"V_z = V_\infty (1 - e^{-z/z_{\rm acc}})"


wa = WindAnalysis("mock.toml", **{"output.dir": "out_04_custom", "fit.n_boot": 0})
res = wa.fit("saturating")
print(res.summary())
