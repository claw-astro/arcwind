"""
Velocity laws: how the axial velocity of the wind varies with height.

A law is a small class with named parameters. Built in:

========================  =====================================  ==========================
name                      formula                                parameters
========================  =====================================  ==========================
``"constant"``            :math:`V_z = V_z`                      ``Vz``
``"power-law"``           :math:`V_z (z/z_{\\rm ref})^{q_v}`      ``Vz`` (at zref), ``qv``
``"broken-power-law"``    smoothly broken power law (below)      ``Vz_b``, ``z_b``, ``q1``, ``q2``
========================  =====================================  ==========================

The smoothly broken power law is

.. math::

   V_z = V_{z,b}\\left(\\frac{z}{z_b}\\right)^{q_1}
         \\left[\\frac{1 + (z/z_b)^{1/\\Delta}}{2}\\right]^{(q_2-q_1)\\Delta},

which behaves as :math:`z^{q_1}` below the break height :math:`z_b` and as
:math:`z^{q_2}` above it (:math:`\\Delta` = ``fit.bpl_delta``).

Add your own law by subclassing :class:`VelocityLaw` and registering it:

>>> import numpy as np
>>> from diskwind.laws import VelocityLaw, register_law
>>> @register_law
... class Exponential(VelocityLaw):
...     '''V = V_inf (1 - exp(-z / z_acc)): accelerates and saturates.'''
...     name = "exponential"
...     param_names = ("V_inf", "z_acc")
...     amplitude = "V_inf"
...     radial_map = {"V_inf": "Vr_inf"}
...     def evaluate(self, z, p):
...         return p["V_inf"] * (1 - np.exp(-z / p["z_acc"]))
...     def default_bounds(self, z_cuts):
...         return {"V_inf": (0.0, 50.0), "z_acc": (1.0, 5 * np.max(z_cuts))}
...     def guess(self, v_ref, z_cuts):
...         return {"V_inf": 1.5 * v_ref, "z_acc": float(np.median(z_cuts))}
>>> "exponential" in available_laws()
True

and then ``cfg.fit.velocity_law = "exponential"``.
"""
from __future__ import annotations

import numpy as np

__all__ = ["VelocityLaw", "ConstantLaw", "PowerLaw", "BrokenPowerLaw", "register_law",
           "get_law", "available_laws"]

_REGISTRY = {}


def register_law(cls):
    """Class decorator: make a :class:`VelocityLaw` available by its ``name``."""
    if not getattr(cls, "name", None):
        raise ValueError("a velocity law needs a 'name'")
    _REGISTRY[cls.name] = cls
    return cls


def available_laws():
    """Names of the registered velocity laws."""
    return sorted(_REGISTRY)


def get_law(name, **kwargs):
    """An instance of the law ``name`` (``kwargs``: e.g. ``zref_au``, ``delta``)."""
    if isinstance(name, VelocityLaw):
        return name
    if name not in _REGISTRY:
        raise ValueError(f"unknown velocity law {name!r}; available: {available_laws()}")
    return _REGISTRY[name](**kwargs)


class VelocityLaw:
    """Base class of velocity laws V(z).

    Subclasses define

    * ``name`` -- registry name;
    * ``param_names`` -- parameter names (unique within the model);
    * ``amplitude`` -- the parameter that scales the velocity (used for the Vsys
      error and for 'proportional' radial velocities);
    * ``radial_map`` -- names of the parameters of an independent radial law
      (``fit.radial = 'independent'``); parameters not in it are shared
      (e.g. the break height);
    * :meth:`evaluate`, :meth:`default_bounds`, :meth:`guess`.

    Parameters
    ----------
    zref_au : float
        Reference height of power laws (au).
    delta : float
        Smoothness of broken power laws.
    """
    name = None
    param_names = ()
    amplitude = None
    radial_map = {}

    def __init__(self, zref_au=100.0, delta=0.1):
        self.zref_au = float(zref_au)
        self.delta = float(delta)

    def evaluate(self, z, p):
        """Velocity at heights ``z`` (au) for parameters ``p`` (mapping name -> value)."""
        raise NotImplementedError

    def __call__(self, z, p):
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            return self.evaluate(np.asarray(z, float), p)

    def shape(self, z, p):
        """V(z) / amplitude (used for 'proportional' radial velocities)."""
        a = p[self.amplitude]
        return self(z, p) / np.where(a == 0, 1.0, a)

    def default_bounds(self, z_cuts):
        """Bounds ``{name: (lo, hi)}`` given the central heights of the cuts (au)."""
        raise NotImplementedError

    def guess(self, v_ref, z_cuts):
        """Starting values given a typical speed ``v_ref`` (km/s) of the traces."""
        raise NotImplementedError

    def starts_from_profile(self, vfunc, z_cuts):
        """Starting points that follow a previous velocity profile ``vfunc(z)``
        (e.g. a simpler fit); default: one :meth:`guess` at the median height."""
        zm = float(np.median(z_cuts))
        return [self.guess(float(vfunc(np.array([zm]))[0]), z_cuts)]

    def latex(self):
        """LaTeX formula of the law."""
        return ""

    def __repr__(self):
        return f"{type(self).__name__}({', '.join(self.param_names)})"


def _local_index(vfunc, z):
    z = np.asarray([0.9 * z, 1.1 * z], float)
    v = np.clip(vfunc(z), 1e-3, None)
    return float(np.log(v[1] / v[0]) / np.log(z[1] / z[0]))


@register_law
class ConstantLaw(VelocityLaw):
    """Constant axial velocity, ``Vz(z) = Vz``."""
    name = "constant"
    param_names = ("Vz",)
    amplitude = "Vz"
    radial_map = {"Vz": "Vr"}

    def evaluate(self, z, p):
        return p["Vz"] * np.ones_like(z)

    def default_bounds(self, z_cuts):
        return {"Vz": (0.0, 50.0)}

    def guess(self, v_ref, z_cuts):
        return {"Vz": float(v_ref)}

    def latex(self):
        return r"V_z = V_z"


@register_law
class PowerLaw(VelocityLaw):
    """``Vz(z) = Vz (z / zref)^qv``; ``Vz`` is the velocity at ``zref``."""
    name = "power-law"
    param_names = ("Vz", "qv")
    amplitude = "Vz"
    radial_map = {"Vz": "Vr", "qv": "qr"}

    def evaluate(self, z, p):
        return p["Vz"] * np.clip(z / self.zref_au, 1e-3, None) ** p["qv"]

    def default_bounds(self, z_cuts):
        return {"Vz": (0.0, 50.0), "qv": (-3.0, 3.0)}

    def guess(self, v_ref, z_cuts):
        return {"Vz": float(v_ref), "qv": 0.0}

    def starts_from_profile(self, vfunc, z_cuts):
        return [{"Vz": float(vfunc(np.array([self.zref_au]))[0]),
                 "qv": _local_index(vfunc, float(np.median(z_cuts)))}]

    def latex(self):
        return rf"V_z = V_z\,(z/{self.zref_au:g}\,{{\rm au}})^{{q_v}}"


@register_law
class BrokenPowerLaw(VelocityLaw):
    """Smoothly broken power law: ``~ z^q1`` below ``z_b``, ``~ z^q2`` above,
    ``= Vz_b`` at ``z_b`` (smoothness ``delta``)."""
    name = "broken-power-law"
    param_names = ("Vz_b", "z_b", "q1", "q2")
    amplitude = "Vz_b"
    radial_map = {"Vz_b": "Vr_b", "q1": "r1", "q2": "r2"}

    def evaluate(self, z, p):
        x = np.clip(z / np.maximum(p["z_b"], 1e-3), 1e-4, None)
        D = self.delta
        return p["Vz_b"] * x ** p["q1"] * ((1.0 + x ** (1.0 / D)) / 2.0) ** \
            ((p["q2"] - p["q1"]) * D)

    def default_bounds(self, z_cuts):
        z = np.asarray(z_cuts, float)
        return {"Vz_b": (0.0, 50.0), "z_b": (0.9 * z.min(), 1.1 * z.max()),
                "q1": (-2.0, 3.0), "q2": (-2.0, 3.0)}

    def guess(self, v_ref, z_cuts):
        return {"Vz_b": float(v_ref), "z_b": float(np.median(z_cuts)), "q1": 0.0, "q2": 0.0}

    def starts_from_profile(self, vfunc, z_cuts):
        """A grid of break heights (30, 50, 70, 85% quantiles of the cut heights) and
        upper indices (q2 = q1 + 0, 0.5, 1), each following ``vfunc`` at the break."""
        out = []
        for zb in np.quantile(np.asarray(z_cuts, float), [0.3, 0.5, 0.7, 0.85]):
            q = _local_index(vfunc, zb)
            for dq in (0.0, 0.5, 1.0):
                out.append({"Vz_b": float(vfunc(np.array([zb]))[0]), "z_b": float(zb),
                            "q1": q, "q2": q + dq})
        return out

    def latex(self):
        return (r"V_z = V_{z,b}(z/z_b)^{q_1}\left[\frac{1+(z/z_b)^{1/\Delta}}{2}\right]"
                rf"^{{(q_2-q_1)\Delta}},\ \Delta={self.delta:g}")
