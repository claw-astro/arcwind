"""Mass figure and diagnostic figures."""
from __future__ import annotations

import numpy as np

from ._style import finish, plt

__all__ = ["plot_mass", "plot_noise_check", "plot_disk_rotation_check"]


def plot_mass(mass, cube=None, save=None):
    """Integrated-intensity map with the slices (left); mean brightness with the
    power-law fit, slice mass and Mdot vs height (right).

    Parameters
    ----------
    mass : MassResult
    cube : WindCube, optional
        Needed for the map panel.

    Returns
    -------
    Figure
    """
    p = plt()
    S = mass.slices
    has_mdot = "mdot" in S
    npan = 3 if has_mdot else 2
    fig = p.figure(figsize=(12 if cube is not None else 7, 3.2 * npan))
    gs = fig.add_gridspec(npan, 2 if cube is not None else 1,
                          width_ratios=[1, 1.4] if cube is not None else None)
    col = 1 if cube is not None else 0
    z = S["z"]
    if cube is not None:
        from ..masks import x_window
        ax0 = fig.add_subplot(gs[:, 0])
        ax0.set_facecolor("k")
        W = mass.W_map
        st = cube.step
        ext = [cube.x[0] - st / 2, cube.x[-1] + st / 2, cube.y[0] - st / 2, cube.y[-1] + st / 2]
        im = ax0.imshow(W, origin="lower", extent=ext, cmap="magma", vmin=0,
                        vmax=np.nanpercentile(W, 99.5), aspect="auto")
        p.colorbar(im, ax=ax0, fraction=0.05, pad=0.1, location="bottom",
                   label=r"$\int T_{\rm mb}\,dV$ (K km/s)")
        ax0.contour(cube.x, cube.y, np.ma.masked_invalid(W / mass.sigW_map),
                    levels=[mass.settings["nsigma"]], colors="w", linewidths=0.7)
        for e in mass.edges:
            ax0.axhline(e, color="0.6", lw=0.4, ls=(0, (4, 3)))
        e0, e1 = mass.edges[0], mass.edges[-1]
        yy = cube.y[(cube.y >= e0) & (cube.y <= e1)]
        xw = np.array([x_window(cube.config.for_cut(v), v) for v in yy])
        for k in (0, 1):
            ax0.plot(xw[:, k], yy, color="0.75", lw=0.8, ls=":")
        ax0.set_ylim(e0 - 0.1, e1 + 0.2)
        ax0.set_xlabel("X (arcsec)")
        ax0.set_ylabel("Y along the wind axis (arcsec)")
    ax1 = fig.add_subplot(gs[0, col])
    ax1.errorbar(z, S["W"], yerr=S["W_err_tot"], fmt="o", ms=4, color="k", capsize=1.5)
    if np.isfinite(mass.alpha) and mass.alpha_fit is not None:
        cf = mass.alpha_fit[0]
        ok = np.isfinite(S["W"]) & (S["W"] > 0)
        zz = np.linspace(z[ok].min(), z[ok].max(), 50)
        ax1.plot(zz, np.exp(cf[1]) * zz ** cf[0], "g-",
                 label=rf"$\propto z^{{-{mass.alpha:.2f}\pm{mass.alpha_err:.2f}}}$")
        ax1.legend(fontsize=7)
    ax1.set_ylabel(r"$\langle\int T_{\rm mb}dV\rangle$ (K km/s)")
    ax2 = fig.add_subplot(gs[1, col], sharex=ax1)
    ax2.fill_between(z, S["mass_tlo"], S["mass_thi"], color="0.8", lw=0, label="T$_{ex}$ range")
    ax2.errorbar(z, S["mass"], yerr=S["mass_err"], fmt="ko-", ms=4, capsize=1.5)
    ax2.set_ylabel(r"slice mass ($M_\odot$)")
    ax2.legend(fontsize=7)
    last = ax2
    if has_mdot:
        ax3 = fig.add_subplot(gs[2, col], sharex=ax1)
        ax3.fill_between(z, S["mdot_tlo"], S["mdot_thi"], color="0.8", lw=0)
        ax3.errorbar(z, S["mdot"], yerr=S["mdot_err"], fmt="bo-", ms=4, capsize=1.5)
        ax3.set_ylabel(r"$\dot M$ ($M_\odot$ yr$^{-1}$)")
        last = ax3
    last.set_xlabel("height z (au)")
    fig.tight_layout()
    return finish(fig, save)


def plot_noise_check(report, save=None):
    """Per-channel rms and fraction of pixels > 4 sigma; PV rms per cut from the
    current, blue-end and red-end channels."""
    p = plt()
    fig, ax = p.subplots(1, 3, figsize=(15, 4.2))
    ax[0].plot(report.v, report.rms_ratio, "k.-", lw=0.8)
    ax[0].set_ylabel("channel rms / median")
    ax[1].semilogy(report.v, np.maximum(report.frac_4sigma, 1e-6), "k.-", lw=0.8)
    ax[1].axhline(6.3e-5, color="grey", ls=":")
    ax[1].set_ylabel(r"fraction of pixels $>4\sigma$")
    for a in ax[:2]:
        a.set_xlabel(r"$V_{\rm LSRK}$ (km s$^{-1}$)")
        for v in report.noise_v:
            a.axvline(v, color="tab:blue", lw=0.6, alpha=0.6)
        if report.wind_window is not None:
            a.axvspan(*report.wind_window, color="tab:red", alpha=0.12, label="wind window")
    if report.wind_window is not None:
        ax[0].legend(fontsize=8)
    ys = np.array(list(report.cuts))
    c = np.array(list(report.cuts.values()))
    ax[2].plot(ys, c[:, 0], "ko-", label="current")
    ax[2].plot(ys, c[:, 1], "s--", color="tab:blue", label="blue end only")
    ax[2].plot(ys, c[:, 2], "^--", color="tab:red", label="red end only")
    ax[2].set_xlabel(r"$y_0$ (arcsec)")
    ax[2].set_ylabel(f"PV rms ({report.units})")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    return finish(fig, save)


def plot_disk_rotation_check(check, save=None):
    """Cuts parallel to the major axis with Keplerian curves for both senses."""
    from matplotlib.colors import PowerNorm
    p = plt()
    n = len(check.pvs)
    fig, axs = p.subplots(1, n, figsize=(4.6 * n, 4.2), squeeze=False)
    for ax, pv, cc in zip(axs[0], check.pvs, check.curves):
        ax.imshow(np.clip(pv.data, 0, None), origin="lower", aspect="auto", extent=pv.extent,
                  cmap="Greys", norm=PowerNorm(0.5, vmin=0, vmax=np.nanmax(pv.data)))
        for (sense, h), v in cc.items():
            ax.plot(check.x, v, ls="-" if sense > 0 else "--",
                    color="tab:red" if sense > 0 else "tab:blue", lw=1.0 if h == 0 else 0.6,
                    alpha=1 if h == 0 else 0.6,
                    label=f"disk_rotation_sense = {sense:+d}" if h == 0 else None)
        ax.axhline(check.vsys, color="0.5", lw=0.5, ls=":")
        ax.set_ylim(check.vsys - 9, check.vsys + 9)
        ax.set_title(f"cut at y = {pv.y0:+.2f}\"", fontsize=9)
        ax.set_xlabel("X (arcsec)")
    axs[0, 0].set_ylabel("V (km/s)")
    axs[0, 0].legend(fontsize=7, loc="upper left")
    fig.suptitle(check.summary(), fontsize=8)
    fig.tight_layout()
    return finish(fig, save)
