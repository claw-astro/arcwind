"""
The conical wind model.

The wind is a thin, axisymmetric conical shell

.. math::

   R(z) = R_0 + z\\tan\\theta, \\qquad x_c(z) = x_0 + \\xi z,

whose gas moves with an axial velocity :math:`V_z(z)` (a :mod:`velocity law
<arcwind.laws>`), a radial velocity :math:`V_r` and a rotation
:math:`V_\\phi = j(z)/R`:

* ``radial = "cone"``: :math:`V_r = V_z\\tan\\theta` (flow along the cone wall);
* ``radial = "proportional"``: :math:`V_r = V_{r}\\,V_z(z)/V_z^{\\rm amp}`
  (free amplitude, same height dependence);
* ``radial = "independent"``: :math:`V_r` follows its own copy of the law;
* ``rotation = "constant"``: :math:`j` constant; ``"power-law"``:
  :math:`j\\,(z/z_{\\rm ref})^{q_j}`.

Parameter names: ``R0_au``, ``tan_theta``, ``xoff_au`` (:math:`x_0`), ``xslope``
(:math:`\\xi`), the law's parameters (e.g. ``Vz_b, z_b, q1, q2``), the radial
parameters (``Vr`` or e.g. ``Vr_b, r1, r2``), ``j`` and ``qj``.

Examples
--------
>>> from arcwind import ConicalWindModel, ViewGeometry
>>> m = ConicalWindModel("broken-power-law", radial="cone")
>>> m.names
['R0_au', 'tan_theta', 'xoff_au', 'xslope', 'Vz_b', 'z_b', 'q1', 'q2', 'j']
>>> m.set(R0_au=30, tan_theta=0.34, xoff_au=30, xslope=-0.08, Vz_b=9.75, z_b=491,
...       q1=0.21, q2=1.25, j=-39)
ConicalWindModel(broken-power-law, radial=cone, rotation=constant)
>>> Vz, Vr = m.velocities([100.0, 500.0])
>>> g = ViewGeometry(inclination_deg=49, s=-1, distance_pc=137.5, vsys=1.9)
>>> c = m.curve(1.0, g)              # model curve of the cut at y0 = 1"
>>> c.X.shape, bool(c.valid.any())
((361,), True)
>>> m.fix(xslope=0.0)                # hold a parameter fixed in fits
ConicalWindModel(broken-power-law, radial=cone, rotation=constant)
>>> m.free_names[:3]
['R0_au', 'tan_theta', 'xoff_au']
"""
from __future__ import annotations

import copy
from dataclasses import dataclass

import numpy as np

from .geometry import cone_cut_heights
from .laws import get_law

__all__ = ["ConicalWindModel", "ModelCurve", "GEOMETRY_PARAMS"]

GEOMETRY_PARAMS = ("R0_au", "tan_theta", "xoff_au", "xslope")
_DEFAULTS = dict(R0_au=5.0, tan_theta=0.5, xoff_au=0.0, xslope=0.0, j=0.0, qj=0.0)
_DESCRIPTIONS = dict(
    R0_au="cone radius at z = 0 (au)", tan_theta="tan of the cone half-opening angle",
    xoff_au="offset of the axis along X at z = 0 (au)",
    xslope="tilt of the axis along X, dx_c/dz", j="specific angular momentum R V_phi (au km/s)",
    qj="index of j(z) = j (z/zref)^qj", Vr="radial velocity amplitude (km/s)")


@dataclass
class ModelCurve:
    """The model in one PV diagram: the locus traced by the cone in the cut.

    Attributes
    ----------
    phi : ndarray
        Azimuth (rad).
    X : ndarray
        Position along the cut (arcsec).
    v : ndarray
        LSRK velocity (km/s).
    z, R : ndarray
        Height and radius (au) at which the cut meets the cone.
    valid : ndarray of bool
        Where the solution exists (z > 0, R > 0).
    """
    phi: np.ndarray
    X: np.ndarray
    v: np.ndarray
    z: np.ndarray
    R: np.ndarray
    valid: np.ndarray

    def __getitem__(self, key):
        return getattr(self, key)

    def masked(self):
        """(X, v) with NaN where invalid, ready for ``ax.plot``."""
        return np.where(self.valid, self.X, np.nan), np.where(self.valid, self.v, np.nan)


class ConicalWindModel:
    """A conical wind: geometry, velocity law, radial-velocity and rotation laws.

    Parameters
    ----------
    velocity_law : str or VelocityLaw
        Name of a registered law (``arcwind.laws.available_laws()``) or an instance.
    radial : {'cone', 'proportional', 'independent'}
        How the radial velocity is set (module docstring).
    rotation : {'constant', 'power-law'}
        Specific angular momentum law.
    zref_au : float
        Reference height of power laws.
    bpl_delta : float
        Smoothness of the broken power law.
    **values
        Initial parameter values.

    Attributes
    ----------
    values : dict
        Current parameter values.
    fixed : dict
        Parameters held fixed in fits, ``{name: value}``.
    bounds : dict
        User bounds, ``{name: (lo, hi)}`` (others: defaults at fit time).
    """

    def __init__(self, velocity_law="power-law", radial="cone", rotation="constant",
                 zref_au=100.0, bpl_delta=0.1, **values):
        if radial not in ("cone", "proportional", "independent"):
            raise ValueError("radial must be 'cone', 'proportional' or 'independent'")
        if rotation not in ("constant", "power-law"):
            raise ValueError("rotation must be 'constant' or 'power-law'")
        self.law = get_law(velocity_law, zref_au=zref_au, delta=bpl_delta)
        self.radial = radial
        self.rotation = rotation
        self.zref_au = float(zref_au)
        names = list(GEOMETRY_PARAMS) + list(self.law.param_names)
        if radial == "proportional":
            names.append("Vr")
        elif radial == "independent":
            names += [self.law.radial_map[k] for k in self.law.param_names
                      if k in self.law.radial_map]
        names.append("j")
        if rotation == "power-law":
            names.append("qj")
        self.names = names
        lg = self.law.guess(5.0, np.array([100.0]))
        self.values = {k: float(_DEFAULTS.get(k, lg.get(k, self._radial_default(k, lg))))
                       for k in names}
        self.fixed = {}
        self.bounds = {}
        self.set(**values)

    @classmethod
    def from_config(cls, cfg):
        """Model as set by ``cfg.fit`` (law, radial, rotation, fixed, bounds)."""
        f = cfg.fit
        m = cls(f.velocity_law, radial=f.radial, rotation=f.rotation, zref_au=f.zref_au,
                bpl_delta=f.bpl_delta)
        for k, v in (f.fixed or {}).items():
            if k in m.names:
                m.fix(**{k: v})
        for k, b in (f.bounds or {}).items():
            if k in m.names:
                m.set_bounds(**{k: b})
        return m

    def _radial_default(self, k, lg):
        inv = {v: kk for kk, v in self.law.radial_map.items()}
        if k in inv:
            return lg[inv[k]] * (0.5 if inv[k] == self.law.amplitude else 1.0)
        return 1.0

    # ------------------------------------------------------------ parameters
    def _check(self, k):
        if k not in self.names:
            raise KeyError(f"{k!r} is not a parameter of this model; parameters: {self.names}")

    def set(self, **values):
        """Set parameter values. Returns the model (so calls can be chained)."""
        for k, v in values.items():
            self._check(k)
            self.values[k] = float(v)
            if k in self.fixed:
                self.fixed[k] = float(v)
        return self

    def fix(self, *names, **values):
        """Hold parameters fixed: ``fix('xslope')`` at the current value or
        ``fix(xslope=0.0)`` at a given value."""
        for k in names:
            self._check(k)
            self.fixed[k] = self.values[k]
        for k, v in values.items():
            self._check(k)
            self.values[k] = float(v)
            self.fixed[k] = float(v)
        return self

    def free(self, *names):
        """Release fixed parameters."""
        for k in names:
            self.fixed.pop(k, None)
        return self

    def set_bounds(self, **bounds):
        """Bounds ``name=(lo, hi)`` used in fits."""
        for k, b in bounds.items():
            self._check(k)
            self.bounds[k] = (float(b[0]), float(b[1]))
        return self

    @property
    def free_names(self):
        """Names of the parameters that are fitted."""
        return [k for k in self.names if k not in self.fixed]

    def vector(self, names=None):
        """Parameter values as an array (in ``names`` order, default all)."""
        return np.array([self.values[k] for k in (names or self.names)], float)

    def copy(self):
        return copy.deepcopy(self)

    #: Parameters that change sign in the frame of the opposite lobe: the axis
    #: offset (X is reversed) and the specific angular momentum (the axis is).
    MIRROR_ODD = ("xoff_au", "j")

    def mirror_signs(self):
        """+1/-1 per parameter (``names`` order): the factor that maps the parameters
        of one lobe onto the frame of the opposite lobe."""
        return np.array([-1.0 if k in self.MIRROR_ODD else 1.0 for k in self.names])

    def mirrored(self):
        """The same wind described in the frame of the opposite lobe.

        A mirror-symmetric bipolar wind has the same cone and velocities in both
        lobes; in the counter lobe's own frame ``xoff_au`` and ``j`` change sign
        (``xslope``, a tilt, does not). Use with :meth:`ViewGeometry.mirrored`.
        """
        m = self.copy()
        for k in self.MIRROR_ODD:
            if k in m.values:
                m.values[k] = -m.values[k]
            if k in m.fixed:
                m.fixed[k] = -m.fixed[k]
            if k in m.bounds:
                lo, hi = m.bounds[k]
                m.bounds[k] = (-hi, -lo)
        return m

    def default_bounds(self, geom, z_cuts):
        """Bounds of every parameter for a geometry and cut heights (user bounds win)."""
        b = dict(R0_au=(0.0, 500.0), tan_theta=(0.0, geom.tan_theta_max),
                 xoff_au=(-500.0, 500.0), xslope=(-2.0, 2.0), j=(-3000.0, 3000.0),
                 qj=(-3.0, 3.0), Vr=(-50.0, 50.0))
        lb = self.law.default_bounds(z_cuts)
        b.update(lb)
        if self.radial == "independent":
            for k, rk in self.law.radial_map.items():
                b[rk] = lb[k]
        b.update(self.bounds)
        return {k: b[k] for k in self.names}

    def _p(self, p):
        """Mapping name -> value from None (current), dict, or array(s) in names order."""
        if p is None:
            return self.values
        if isinstance(p, dict):
            out = dict(self.values)
            out.update(p)
            return out
        p = np.asarray(p, float)
        return {k: p[..., i] for i, k in enumerate(self.names)}

    # -------------------------------------------------------------- physics
    def velocities(self, z, p=None):
        """Axial and radial velocity (km/s) at heights ``z`` (au).

        With ``p`` an array of draws of shape (n, k), returns arrays (n, len(z)).
        """
        P = self._p(p)
        z = np.asarray(z, float)
        if np.ndim(P[self.names[0]]) > 0:
            P = {k: np.asarray(v)[..., None] for k, v in P.items()}
        Vz = self.law(z, P)
        if self.radial == "cone":
            Vr = Vz * P["tan_theta"]
        elif self.radial == "proportional":
            Vr = P["Vr"] * self.law.shape(z, P)
        else:
            Q = dict(P)
            for k, rk in self.law.radial_map.items():
                Q[k] = P[rk]
            Vr = self.law(z, Q)
        return Vz, Vr

    def j_of_z(self, z, p=None):
        """Specific angular momentum (au km/s) at heights ``z``."""
        P = self._p(p)
        z = np.asarray(z, float)
        if np.ndim(P["j"]) > 0:
            P = {k: np.asarray(v)[..., None] for k, v in P.items()}
        if self.rotation == "constant":
            return P["j"] * np.ones_like(z)
        return P["j"] * np.clip(z / self.zref_au, 1e-3, None) ** P["qj"]

    def radius(self, z, p=None):
        """Cone radius R(z) (au)."""
        P = self._p(p)
        return P["R0_au"] + np.asarray(z, float) * P["tan_theta"]

    def profiles(self, z, p=None):
        """All model quantities along the flow.

        Returns
        -------
        dict of ndarray
            ``z, R, xc, Vz, Vr, Vp, Vphi, j, flow_angle_deg`` (au, km/s, deg).
        """
        P = self._p(p)
        z = np.asarray(z, float)
        Vz, Vr = self.velocities(z, p)
        jz = self.j_of_z(z, p)
        nd = np.ndim(P["R0_au"]) > 0
        g = (lambda k: np.asarray(P[k])[..., None]) if nd else (lambda k: P[k])
        R = g("R0_au") + z * g("tan_theta")
        return dict(z=z, R=R, xc=g("xoff_au") + g("xslope") * z, Vz=Vz, Vr=Vr,
                    Vp=np.hypot(Vz, Vr), Vphi=jz / R, j=jz,
                    flow_angle_deg=np.degrees(np.arctan2(Vr, Vz)))

    def curve(self, y0, geom, p=None, nphi=361):
        """Locus of the model in the PV diagram of the cut at ``y0`` (arcsec).

        For each azimuth the height where the cone meets the cut is solved
        exactly (:func:`arcwind.geometry.cone_cut_heights`).

        Returns
        -------
        ModelCurve
        """
        P = self._p(p)
        si, ci, s, d = geom.sini, geom.cosi, geom.s, geom.distance_pc
        phi = np.linspace(0.0, 2 * np.pi, nphi)
        sp, cp = np.sin(phi), np.cos(phi)
        R0, T = P["R0_au"], P["tan_theta"]
        z, denom = cone_cut_heights(y0 * d, R0, T, phi, geom)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            R = R0 + z * T
            Vz, Vr = self.velocities(np.where(z > 0, z, 1.0), P)
            X = (P["xoff_au"] + P["xslope"] * z + R * cp) / d
            v = geom.vsys - s * Vz * ci + Vr * si * sp + self.j_of_z(z, P) / R * si * cp
        valid = (denom > 0) & (z > 0) & (R > 0) & np.isfinite(v)
        return ModelCurve(phi=phi, X=X, v=v, z=z, R=R, valid=valid)

    def ring_at(self, y0, geom, p=None):
        """The model at the central height of the cut, as ring quantities.

        Returns
        -------
        dict
            ``zc, R, xc`` (au), ``vcent`` (km/s), ``A = Vr sin i``,
            ``B = Vphi sin i`` (km/s) and ``T = tan theta``.
        """
        P = self._p(p)
        zc = float(geom.zc_au(y0))
        T = P["tan_theta"]
        R = P["R0_au"] + T * zc
        Vz, Vr = self.velocities(np.array([zc]), P)
        return dict(zc=zc, R=R, xc=P["xoff_au"] + P["xslope"] * zc,
                    vcent=geom.vsys - geom.s * float(Vz[0]) * geom.cosi,
                    A=float(Vr[0]) * geom.sini,
                    B=float(self.j_of_z(zc, P)) / R * geom.sini, T=T)

    def starts(self, geom, z_cuts, v_ref, x_ref_au, vr_ref=None, n_theta=6, previous=None):
        """Starting points for a fit (list of dicts).

        Without ``previous``: a grid of ``n_theta`` cone openings between 0.1 and 0.9
        times the maximum, with velocities from ``v_ref``, ``R0 = 5 au``,
        ``x0 = x_ref_au`` and ``j = 0``. With ``previous`` (a model or result):
        its geometry and rotation and starting points of this law that follow
        its velocity profile.
        """
        out = []
        if previous is None:
            lg = self.law.guess(v_ref, z_cuts)
            for T in np.linspace(0.1, 0.9, n_theta) * geom.tan_theta_max:
                st = dict(self.values)
                st.update(lg)
                st.update(tan_theta=T, R0_au=5.0, xoff_au=x_ref_au, j=0.0)
                self._radial_start(st, None, lg, vr_ref=vr_ref)
                out.append(st)
            return out
        pm = previous.model if hasattr(previous, "model") else previous
        base = dict(self.values)
        for k in ("R0_au", "tan_theta", "xoff_au", "xslope", "j", "qj"):
            if k in pm.values and k in self.values:
                base[k] = pm.values[k]
        vz_prev = (lambda z: pm.velocities(z)[0])
        vr_prev = (lambda z: pm.velocities(z)[1])
        same = (pm.law.name == self.law.name)
        for lg in ([dict((k, pm.values[k]) for k in self.law.param_names)] if same
                   else self.law.starts_from_profile(vz_prev, z_cuts)):
            st = dict(base)
            st.update(lg)
            self._radial_start(st, vr_prev, lg, from_profile=True)
            out.append(st)
        return out

    def _radial_start(self, st, vr_func, lg, from_profile=False, vr_ref=None):
        """Starting values of the radial-velocity parameters (in place).

        From a previous profile ``vr_func``, the radial amplitude is chosen so that
        Vr matches it at the law's reference height; otherwise it is ``vr_ref``.
        """
        if self.radial == "cone":
            return
        z_ref = np.array([st.get("z_b", self.zref_au)])
        if from_profile:
            # amplitude = Vr(z_ref) / shape(z_ref), so that the law gives Vr(z_ref)
            amp = float(vr_func(z_ref)[0] / max(float(self.law.shape(z_ref, st)[0]), 1e-6))
        else:
            amp = float(0.5 * lg[self.law.amplitude] if vr_ref is None else vr_ref)
        if self.radial == "proportional":
            st["Vr"] = amp
        else:                                   # independent: a copy of the law
            for k, rk in self.law.radial_map.items():
                st[rk] = amp if k == self.law.amplitude else lg[k]

    # ---------------------------------------------------------------- misc
    def describe(self):
        """Text table of the parameters, their values and status."""
        lines = [f"{self!r}", f"  law: {self.law.latex()}"]
        for k in self.names:
            st = "fixed" if k in self.fixed else "free"
            b = self.bounds.get(k)
            lines.append(f"  {k:10s} = {self.values[k]:10.4g}  {st:5s}"
                         + (f"  bounds {b}" if b else "")
                         + f"  {_DESCRIPTIONS.get(k, '')}")
        return "\n".join(lines)

    def __repr__(self):
        return (f"ConicalWindModel({self.law.name}, radial={self.radial}, "
                f"rotation={self.rotation})")
