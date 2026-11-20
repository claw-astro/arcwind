"""
Derived physical quantities of a fitted wind, with uncertainties from the
bootstrap (or covariance) draws of the fit.

Examples
--------
>>> props = result.derived(mass=mass_result)                   # doctest: +SKIP
>>> props.at([64, 446])["Vp"]           # median, 16th, 84th percentiles  # doctest: +SKIP
>>> props.launch_radius()               # r_L = j^2 / G M*                 # doctest: +SKIP
>>> props.travel_time(64, 446)          # yr                               # doctest: +SKIP
>>> print(props.summary())                                                  # doctest: +SKIP
>>> print(props.compare_with("HH 30"))                                     # doctest: +SKIP
"""
from __future__ import annotations

import numpy as np

from .constants import AU_KM, GM_SUN_AU_KMS2, MSUN_G, YR_S

__all__ = ["WindProperties", "REFERENCE_WINDS"]
_trapz = getattr(np, "trapezoid", None) or np.trapz

#: Published winds to compare with (values as given in the papers).
REFERENCE_WINDS = {
    "HH 30": dict(ref="Louvet et al. 2018, A&A 618, A120", distance_pc=140.0, mstar_msun=0.45,
                  inclination_deg=91.0, theta_deg=(35.0, 1.0), R0_au_max=22.0,
                  Vp=(9.3, 0.7), j=(38.0, 15.0), z_range_au=(15, 250), mass_min=1.7e-5,
                  mdot_min=9e-8, alpha=0.8,
                  notes="Vp constant with z, measured at z <= 250 au; CO detected to ~700 au"),
}


def _pct(a, axis=0):
    md, lo, hi = np.nanpercentile(np.asarray(a, float), [50, 16, 84], axis=axis)
    return md, lo, hi


class WindProperties:
    """Physical quantities derived from a :class:`~diskwind.fitting.FitResult`.

    Parameters
    ----------
    result : FitResult
    mass : MassResult, optional
        Slice masses for the mass-loss, momentum and energy rates.
    mstar_msun : float, optional
        Stellar mass (default: from the fit geometry).
    """

    def __init__(self, result, mass=None, mstar_msun=None):
        self.result = result
        self.model = result.model
        self.mass = mass
        self.mstar = mstar_msun or result.geometry.mstar_msun
        self.draws = result.draws()
        zc = [result.geometry.zc_au(y) for y in result.traces] or [100.0]
        self.z_cuts = np.array(zc, float)

    def at(self, z):
        """Quantities at heights ``z`` (au): ``(median, p16, p84)`` arrays for
        ``Vz, Vr, Vp, flow_angle_deg, R, Vphi, index`` (index = dlnVp/dlnz)."""
        z = np.atleast_1d(np.asarray(z, float))
        P = self.model.profiles(z, self.draws)
        zz = np.concatenate([0.999 * z, 1.001 * z])
        Q = self.model.profiles(zz, self.draws)["Vp"]
        n = len(z)
        idx = (np.log(Q[:, n:]) - np.log(Q[:, :n])) / (np.log(1.001) - np.log(0.999))
        out = {k: _pct(P[k]) for k in ("Vz", "Vr", "Vp", "flow_angle_deg", "R", "Vphi")}
        out["index"] = _pct(idx)
        return out

    def launch_radius(self):
        """r_L = j^2 / (G M*) (au), the launch radius if j is conserved from the disk
        without torque: (median, p16, p84)."""
        if not self.mstar:
            raise ValueError("needs the stellar mass")
        j = self.draws[:, self.model.names.index("j")]
        return _pct(j ** 2 / (GM_SUN_AU_KMS2 * self.mstar))

    def escape_ratio(self, z=None):
        """Vp / v_esc at the lowest cut (or at ``z``), with v_esc at r = (R^2 + z^2)^1/2."""
        z = float(self.z_cuts.min() if z is None else z)
        P = self.model.profiles(np.array([z]), self.draws)
        r = np.hypot(P["R"][:, 0], z)
        return _pct(P["Vp"][:, 0] / np.sqrt(2 * GM_SUN_AU_KMS2 * self.mstar / r))

    def travel_time(self, z1=None, z2=None, n=400):
        """Time (yr) to travel from ``z1`` to ``z2`` (au), integral of dz / Vz."""
        z1 = float(self.z_cuts.min() if z1 is None else z1)
        z2 = float(self.z_cuts.max() if z2 is None else z2)
        z = np.linspace(z1, z2, n)
        Vz = self.model.velocities(z, self.draws)[0]
        t = _trapz(1.0 / np.clip(Vz, 1e-3, None), z, axis=1) * AU_KM / YR_S
        return _pct(t)

    def dynamical_time(self, z=None):
        """z / Vz(z) (yr) at ``z`` (default: the highest cut)."""
        z = float(self.z_cuts.max() if z is None else z)
        Vz = self.model.velocities(np.array([z]), self.draws)[0][:, 0]
        return _pct(z * AU_KM / np.clip(Vz, 1e-3, None) / YR_S)

    def rates(self):
        """From the slice masses: per-slice Mdot (Msun/yr), momentum rate Mdot Vp
        (Msun km/s / yr) and kinetic luminosity Mdot Vp^2 / 2 (Lsun); medians over
        slices with emission. Lower limits like the masses."""
        if self.mass is None or "mdot" not in self.mass.slices:
            raise ValueError("needs a MassResult computed with a velocity model")
        S = self.mass.slices
        has = S["mass"] > 0
        Vp = self.model.profiles(S["z"])["Vp"]
        mdot = S["mdot"]
        pdot = mdot * Vp
        lkin = 0.5 * mdot * MSUN_G / YR_S * (Vp * 1e5) ** 2 / 3.828e33
        return dict(z=S["z"], mdot=mdot, pdot=pdot, lkin_lsun=lkin,
                    median_mdot=float(np.median(mdot[has])) if has.any() else np.nan,
                    median_pdot=float(np.median(pdot[has])) if has.any() else np.nan,
                    median_lkin_lsun=float(np.median(lkin[has])) if has.any() else np.nan)

    def summary(self):
        """Text summary of the derived quantities."""
        f = lambda t, fm="{:.2f}": (fm + " (+" + fm + "/-" + fm + ")").format(
            t[0], t[2] - t[0], t[0] - t[1])
        zlo, zhi = self.z_cuts.min(), self.z_cuts.max()
        A = self.at([zlo, zhi])
        L = [f"# derived from {self.result!r}",
             f"Vp at z = {zlo:.0f} au: {f(tuple(a[0] for a in A['Vp']))} km/s; at "
             f"{zhi:.0f} au: {f(tuple(a[1] for a in A['Vp']))} km/s",
             f"flow angle atan(Vr/Vz) at {zlo:.0f} au: "
             f"{f(tuple(a[0] for a in A['flow_angle_deg']), '{:.1f}')} deg",
             f"dlnVp/dlnz at {zlo:.0f} / {zhi:.0f} au: "
             f"{f(tuple(a[0] for a in A['index']))} / {f(tuple(a[1] for a in A['index']))}",
             f"travel time {zlo:.0f} -> {zhi:.0f} au: {f(self.travel_time(), '{:.0f}')} yr; "
             f"z/Vz at {zhi:.0f} au: {f(self.dynamical_time(), '{:.0f}')} yr"]
        if self.mstar:
            L.append(f"r_L = j^2/GM* = {f(self.launch_radius(), '{:.2f}')} au;  "
                     f"Vp/v_esc at {zlo:.0f} au = {f(self.escape_ratio())}")
        if self.mass is not None and "mdot" in self.mass.slices:
            r = self.rates()
            L.append(f"median Mdot = {r['median_mdot']:.3g} Msun/yr, momentum rate = "
                     f"{r['median_pdot']:.3g} Msun km/s/yr, L_kin = {r['median_lkin_lsun']:.3g} "
                     "Lsun (lower limits)")
        return "\n".join(L)

    def compare_with(self, name="HH 30"):
        """Side-by-side text table with a published wind in :data:`REFERENCE_WINDS`."""
        h = REFERENCE_WINDS[name]
        th, eth = self.result.theta_deg
        j = self.result.params["j"]
        ej = self.result.errors()["j"]
        zlo, zhi = self.z_cuts.min(), self.z_cuts.max()
        A = self.at([zlo, zhi])
        rows = [("quantity", "this source", name),
                ("cone half-angle (deg)", f"{th:.1f} +/- {eth:.1f}",
                 f"{h['theta_deg'][0]:.0f} +/- {h['theta_deg'][1]:.0f}"),
                ("Vp (km/s)", f"{A['Vp'][0][0]:.1f} - {A['Vp'][0][1]:.1f}",
                 f"{h['Vp'][0]:.1f} +/- {h['Vp'][1]:.1f}"),
                ("heights (au)", f"{zlo:.0f} - {zhi:.0f}",
                 f"{h['z_range_au'][0]} - {h['z_range_au'][1]}"),
                ("|j| (au km/s)", f"{abs(j):.0f} +/- {ej:.0f}",
                 f"{h['j'][0]:.0f} +/- {h['j'][1]:.0f}")]
        if self.mstar:
            rl = self.launch_radius()
            rows.append(("r_L = j^2/GM* (au)", f"{rl[0]:.1f}",
                         f"{h['j'][0] ** 2 / (GM_SUN_AU_KMS2 * h['mstar_msun']):.1f}"))
        if self.mass is not None:
            rows.append(("mass (Msun)", f">= {self.mass.total_mass:.2g}", f">= {h['mass_min']:.2g}"))
            if np.isfinite(self.mass.median_mdot):
                rows.append(("Mdot (Msun/yr)", f">= {self.mass.median_mdot:.2g}",
                             f">= {h['mdot_min']:.2g}"))
        w = [max(len(r[i]) for r in rows) for i in range(3)]
        out = [f"{a:{w[0]}s}  {b:{w[1]}s}  {c:{w[2]}s}" for a, b, c in rows]
        out.append(f"# {name}: {h['ref']}; {h['notes']}")
        return "\n".join(out)
