"""
Brightness, mass and mass-loss rate of the wind in slices along its axis
(Louvet et al. 2018, Sect. 4.2).

For every pixel of the resampled cube the emission is integrated over the
velocity window of the trace and inside its spatial window, converted to
:math:`W=\\int T_{\\rm mb}\\,dV` (K km/s, Rayleigh-Jeans), and kept where
:math:`W > n\\sigma_W`. For optically thin emission in LTE

.. math::

   N_{\\rm up} = \\frac{8\\pi k\\nu^2}{hc^3A_{ul}}W, \\qquad
   N_{\\rm CO} = \\frac{N_{\\rm up}\\,Q(T_{\\rm ex})}{g_{\\rm up}
   e^{-E_{\\rm up}/kT_{\\rm ex}}}, \\qquad
   M = \\sum \\frac{N_{\\rm CO}}{X_{\\rm CO}}\\,\\mu m_{\\rm H_2}\\,\\Omega_{\\rm pix}d^2 .

The mass-loss rate of a slice of thickness :math:`\\Delta z` is
:math:`\\dot M = M V_z(z)/\\Delta z`. Masses and rates are lower limits when part
of the emission is undetected (half rings, absorption, filtering, opacity).

Examples
--------
>>> from arcwind import MassAnalysis
>>> ma = MassAnalysis(cube, velocity_model=result)            # doctest: +SKIP
>>> m = ma.run()                                               # doctest: +SKIP
>>> m.total_mass, m.total_mass_err, m.alpha                    # doctest: +SKIP
>>> m.slices["mdot"]                    # per slice, Msun/yr   # doctest: +SKIP
>>> print(f"{column_density_co(38.7, 30.0) / 1e-4:.3g}")   # HH 30: N_H2 = 2.25e20 cm^-2
2.25e+20
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .constants import ARCSEC_PER_RAD, AU_CM, C_CGS, H_P, K_B, M_H, MSUN_G, PC_CM, YR_S
from .masks import x_window

__all__ = ["MassAnalysis", "MassResult", "column_density_up", "column_density_co",
           "partition_function"]

_DEF = dict(nu=230.538e9, a_ul=6.910e-7, g_up=5.0, e_up=16.596, b_rot=57.635968e9)


def column_density_up(W_kkms, nu=_DEF["nu"], a_ul=_DEF["a_ul"]):
    """Upper-level column density (cm^-2) from ``W`` (K km/s), optically thin."""
    return 8 * np.pi * K_B * nu ** 2 / (H_P * C_CGS ** 3 * a_ul) * np.asarray(W_kkms) * 1e5


def partition_function(tex, b_rot=_DEF["b_rot"]):
    """Rigid-rotor partition function of a linear molecule at ``tex`` (K)."""
    J = np.arange(0, 200)
    return float(np.sum((2 * J + 1) * np.exp(-H_P * b_rot * J * (J + 1) / (K_B * tex))))


def column_density_co(W_kkms, tex, nu=_DEF["nu"], a_ul=_DEF["a_ul"], g_up=_DEF["g_up"],
                      e_up=_DEF["e_up"], b_rot=_DEF["b_rot"]):
    """Total CO column density (cm^-2) from ``W`` (K km/s) at ``tex`` (K), LTE, thin."""
    return (column_density_up(W_kkms, nu, a_ul) * partition_function(tex, b_rot)
            / (g_up * np.exp(-e_up / tex)))


@dataclass
class MassResult:
    """Slices along the axis and totals.

    Attributes
    ----------
    slices : dict of ndarray
        ``y`` (arcsec), ``z`` (au), ``W`` (mean K km/s), ``W_err`` (noise),
        ``W_err_tot`` (+ calibration), ``npix``, ``mass``, ``mass_err``,
        ``mass_tlo``, ``mass_thi`` (T_ex range), and with a velocity model
        ``Vz``, ``mdot``, ``mdot_err``, ``mdot_tlo``, ``mdot_thi`` (Msun, Msun/yr).
    total_mass, total_mass_err : float
    median_mdot, median_mdot_err : float
    alpha, alpha_err : float
        Brightness index, W ~ z^-alpha.
    W_map, sigW_map : ndarray
        Integrated-intensity map and its noise on the (y, x) grid.
    edges : ndarray
        Slice edges (arcsec).
    """
    slices: dict
    total_mass: float
    total_mass_err: float
    total_mass_trange: tuple
    alpha: float
    alpha_err: float
    alpha_fit: tuple = None
    median_mdot: float = np.nan
    median_mdot_err: float = np.nan
    median_mdot_trange: tuple = (np.nan, np.nan)
    W_map: np.ndarray = None
    sigW_map: np.ndarray = None
    edges: np.ndarray = None
    settings: dict = field(default_factory=dict)

    def summary(self):
        s = self.settings
        L = [f"# T_ex = {s.get('tex')} K, X_CO = {s.get('xco'):g}, mu = {s.get('mu')}, "
             f"slices dz = {s.get('dz')} au, threshold {s.get('nsigma')} sigma",
             "# LOWER LIMITS if part of the emission is undetected",
             f"total mass = {self.total_mass:.3g} +/- {self.total_mass_err:.2g} Msun "
             f"(T_ex range {self.total_mass_trange[0]:.3g}-{self.total_mass_trange[1]:.3g})",
             f"brightness ~ z^-alpha, alpha = {self.alpha:.2f} +/- {self.alpha_err:.2f}"]
        if np.isfinite(self.median_mdot):
            L.append(f"median mass-loss rate = {self.median_mdot:.3g} +/- "
                     f"{self.median_mdot_err:.2g} Msun/yr (T_ex range "
                     f"{self.median_mdot_trange[0]:.3g}-{self.median_mdot_trange[1]:.3g})")
        return "\n".join(L)

    def save(self, path):
        """Write the slice table (text)."""
        keys = [k for k in ("y", "z", "W", "W_err", "W_err_tot", "npix", "mass", "mass_err",
                            "mass_tlo", "mass_thi", "Vz", "mdot", "mdot_err", "mdot_tlo",
                            "mdot_thi") if k in self.slices]
        with open(path, "w") as fh:
            fh.write("\n".join("# " + x for x in self.summary().splitlines()) + "\n")
            fh.write("# " + " ".join(keys) + "\n")
            for i in range(len(self.slices["y"])):
                fh.write(" ".join(f"{self.slices[k][i]:.5g}" for k in keys) + "\n")
        return path

    def plot(self, cube=None, **kw):
        """Figure of the map, brightness, mass and Mdot (see :func:`plotting.plot_mass`)."""
        from .plotting import plot_mass
        return plot_mass(self, cube=cube, **kw)


class MassAnalysis:
    """Mass and mass-loss rate of the wind (settings: ``[mass]``).

    Parameters
    ----------
    cube : WindCube
    velocity_model : FitResult or ConicalWindModel, optional
        Provides Vz(z) for the mass-loss rate (and its bootstrap draws for errors).
    **options
        Any ``[mass]`` option for this analysis, e.g. ``tex_k=40``.
    """

    def __init__(self, cube, velocity_model=None, **options):
        from copy import copy
        self.cube = cube
        self.cfg = cube.config
        self.opt = copy(self.cfg.mass)
        for k, v in options.items():
            if not hasattr(self.opt, k):
                raise ValueError(f"unknown mass option {k!r}")
            setattr(self.opt, k, v)
        self.vmodel = velocity_model

    def ncol_factor(self, tex=None):
        """N_H2 (cm^-2) per K km/s."""
        o = self.opt
        return column_density_co(1.0, tex or o.tex_k, o.nu_hz, o.a_ul, o.g_up, o.e_up_k,
                                 o.b_rot_hz) / o.xco

    def moment_map(self, y_range=None):
        """Integrated intensity W (K km/s) over the trace windows, and its noise."""
        c, cfg = self.cube, self.cfg
        to_K = c.k_per_jy / c.ufac
        rms_K = c.rms * to_K
        vrel = c.v - c.vsys
        W = np.full(c.cube.shape[1:], np.nan)
        sig = np.full(c.cube.shape[1:], np.nan)
        lo_y, hi_y = y_range if y_range is not None else (-np.inf, np.inf)
        for iy, yy in enumerate(c.y):
            if yy < lo_y or yy > hi_y:
                continue
            cc = cfg.for_cut(yy)
            vlo, vhi = cc.mask.vrel_range_kms or (-np.inf, np.inf)
            ch = (vrel >= vlo) & (vrel <= vhi)
            for a, b in cc.mask.v_exclude_kms or []:
                ch &= ~((vrel >= a) & (vrel <= b))
            xlo, xhi = x_window(cc, yy)
            xs = (c.x >= xlo) & (c.x <= xhi)
            spec = c.cube[ch][:, iy, :] * to_K
            W[iy, xs] = (np.nansum(spec, axis=0) * c.dv)[xs]
            sig[iy, xs] = rms_K * c.dv * np.sqrt(ch.sum())
        return W, sig

    def run(self):
        """Compute the slices and totals. Returns :class:`MassResult`."""
        c, o = self.cube, self.opt
        g = c.geometry
        rng = np.random.default_rng(self.cfg.fit.seed)
        ys = self.cfg.y_offsets
        if o.y_range_arcsec is not None:
            y0, y1 = o.y_range_arcsec
        else:
            y0, y1 = max(ys.min() - 0.5 * c.slit_width, 0.0), ys.max()
        dy = o.dz_au * g.sini / g.distance_pc
        edges = np.arange(y0, y1 + 1e-9, dy)
        W, sigW = self.moment_map((y0 - dy, y1 + dy))
        pix_sr = (c.step / ARCSEC_PER_RAD) ** 2
        ncol = self.ncol_factor()
        m_pix = pix_sr * (g.distance_pc * PC_CM) ** 2 * o.mu * 2 * M_H / MSUN_G
        nbeam = 1.1331 * c.bmaj * c.bmin / c.step ** 2
        rows = []
        for a, b in zip(edges[:-1], edges[1:]):
            sel = (c.y >= a) & (c.y < b)
            Wy, Sy = W[sel], sigW[sel]
            good = np.isfinite(Wy) & (Wy > o.nsigma * Sy)
            npx = int(good.sum())
            if npx == 0:
                rows.append([0.5 * (a + b), np.nan, np.nan, 0, 0.0, 0.0])
                continue
            s2 = np.mean(Sy[good] ** 2)
            rows.append([0.5 * (a + b), float(np.mean(Wy[good])),
                         float(np.sqrt(s2 * min(npx, nbeam) / npx)), npx,
                         float(np.sum(Wy[good]) * ncol * m_pix),
                         float(np.sqrt(s2 * npx * min(npx, nbeam)) * ncol * m_pix)])
        yc, mW, eW, npx, mass, e_noise = np.array(rows).T
        zc = yc / g.sini * g.distance_pc
        cal = o.flux_cal
        e_mass = np.hypot(e_noise, cal * mass)
        tg = np.linspace(min(o.tex_range_k), max(o.tex_range_k), 200)
        tfac = np.array([self.ncol_factor(t) for t in tg]) / ncol
        S = dict(y=yc, z=zc, W=mW, W_err=eW, W_err_tot=np.hypot(eW, cal * mW), npix=npx,
                 mass=mass, mass_err=e_mass, mass_tlo=mass * tfac.min(),
                 mass_thi=mass * tfac.max())
        nmc = o.nmc
        Md = (mass[None] + e_noise[None] * rng.standard_normal((nmc, len(mass))))
        Md *= 1.0 + cal * rng.standard_normal((nmc, 1))
        out = dict(slices=S, total_mass=float(np.nansum(mass)),
                   total_mass_err=float(np.std(Md.sum(axis=1))),
                   total_mass_trange=(float(np.nansum(mass) * tfac.min()),
                                      float(np.nansum(mass) * tfac.max())))
        if self.vmodel is not None:
            vm = self.vmodel
            model = vm.model if hasattr(vm, "model") else vm
            Vz = model.velocities(zc)[0]
            draws = vm.draws() if hasattr(vm, "draws") else None
            mdot = mass * Vz * 1e5 * YR_S / (o.dz_au * AU_CM)
            if draws is not None:
                Vd = model.velocities(zc, draws)[0]
                eVz = np.nanstd(Vd, axis=0)
            else:
                Vd, eVz = Vz[None], np.zeros_like(Vz)
            with np.errstate(invalid="ignore", divide="ignore"):
                e_mdot = mdot * np.sqrt((e_mass / np.where(mass > 0, mass, np.nan)) ** 2
                                        + (eVz / Vz) ** 2)
            S.update(Vz=Vz, mdot=mdot, mdot_err=np.where(mass > 0, e_mdot, 0.0),
                     mdot_tlo=mdot * tfac.min(), mdot_thi=mdot * tfac.max())
            has = mass > 0
            if has.any():
                Vs = Vd[rng.integers(0, len(Vd), nmc)]
                md = Md * Vs * 1e5 * YR_S / (o.dz_au * AU_CM)
                med = float(np.median(mdot[has]))
                out.update(median_mdot=med,
                           median_mdot_err=float(np.std(np.median(md[:, has], axis=1))),
                           median_mdot_trange=(med * tfac.min(), med * tfac.max()))
        ok = np.isfinite(mW) & (mW > 0)
        alpha = e_alpha = np.nan
        fitpar = None
        if ok.sum() >= 3:
            if ok.sum() >= 4:
                cf, cc = np.polyfit(np.log(zc[ok]), np.log(mW[ok]), 1, cov=True)
                e_alpha = float(np.sqrt(cc[0, 0]))
                fitpar = (cf, cc)
            else:
                cf = np.polyfit(np.log(zc[ok]), np.log(mW[ok]), 1)
                fitpar = (cf, None)
            alpha = -float(cf[0])
        return MassResult(alpha=alpha, alpha_err=e_alpha, alpha_fit=fitpar, W_map=W,
                          sigW_map=sigW, edges=edges,
                          settings=dict(tex=o.tex_k, xco=o.xco, mu=o.mu, dz=o.dz_au,
                                        nsigma=o.nsigma, flux_cal=cal,
                                        tex_range=tuple(o.tex_range_k)), **out)
