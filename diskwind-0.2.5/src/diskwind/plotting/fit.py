"""Figures of fit results: all cuts, one cut, bootstrap, profiles, ring fits."""
from __future__ import annotations

import numpy as np

from matplotlib.ticker import MaxNLocator

from ._style import cut_label, draw_curve, finish, grid, plt, tidy_grid_ticks

__all__ = ["plot_cuts", "plot_cut", "plot_bootstrap", "plot_profiles", "plot_ring_fits"]

_COLORS = ["tab:red", "tab:blue", "tab:green", "tab:purple", "tab:orange"]


def _draw_cut(ax, cube, y0, results, trace, n_boot, show_mask):
    pv = cube.pv(y0, with_mask=show_mask)
    ax.imshow(pv.data, origin="lower", aspect="auto", extent=pv.extent, cmap="Greys",
              vmin=0, vmax=0.8 * np.nanmax(pv.data))
    if show_mask and pv.excluded is not None:
        ax.contour(pv.x, pv.v, pv.excluded.astype(float), levels=[0.5], colors="c",
                   linewidths=0.4)
    main = results[0]
    if n_boot and main.boot is not None:
        for b in main.boot[:n_boot]:
            draw_curve(ax, main.curve(y0, b), color="gold", lw=0.4, alpha=0.4)
    if trace is not None and trace.n:
        ax.errorbar(trace.x, trace.v, yerr=trace.err, fmt="+", ms=3, color="k", lw=0.5)
    for k, r in enumerate(results):
        draw_curve(ax, r.curve(y0), color=_COLORS[k % len(_COLORS)], lw=1.3 if k == 0 else 0.9,
                   ls="-" if k == 0 else "--", label=r.label)


def plot_cuts(cube, result, others=(), n_boot_curves=30, show_mask=True, ncols=4, save=None):
    """Every fitted cut with its trace and the model (and other models dashed).

    Parameters
    ----------
    cube : WindCube
    result : FitResult
        Main model (red); its traces are shown; ``n_boot_curves`` bootstrap curves (gold).
    others : sequence of FitResult
        Models to compare (dashed).

    Returns
    -------
    Figure
    """
    res = [result] + list(others)
    ys = list(result.traces) or list(cube.config.y_offsets)
    fig, axs = grid(len(ys), ncols)
    for ax, y0 in zip(axs.flat, ys):
        _draw_cut(ax, cube, y0, res, result.traces.get(y0), n_boot_curves, show_mask)
        ax.text(0.03, 0.95, cut_label(y0, cube.zc_au(y0)), transform=ax.transAxes, va="top",
                fontsize="x-small")
    pc = cube.config.plot
    axs.flat[0].set_xlim(*(pc.xlim or (cube.x[0], cube.x[-1])))
    if pc.vlim:
        axs.flat[0].set_ylim(*pc.vlim)
    axs.flat[0].legend(fontsize="xx-small", loc="lower left")
    tidy_grid_ticks(axs)
    fig.supxlabel("X (arcsec)")
    fig.supylabel(r"$V_{\rm LSR}$ (km s$^{-1}$)")
    fig.suptitle(f"{cube.config.source.name}: {result.label} (red)"
                 + ("; yellow: bootstrap" if result.boot is not None else ""), fontsize="small")
    fig.subplots_adjust(wspace=0.04, hspace=0.04)
    return finish(fig, save, cube.config.output.dpi)


def plot_cut(cube, result, y0, others=(), n_boot_curves=60, save=None):
    """Zoom on one cut: PV, trace, model(s), bootstrap curves.

    Returns
    -------
    Figure
    """
    p = plt()
    fig, ax = p.subplots(figsize=(6, 4.5))
    tr = result.traces.get(float(y0))
    if tr is None:
        from ..trace import TraceExtractor
        tr = TraceExtractor(cube).extract(y0, guide=result.curve(y0))
    _draw_cut(ax, cube, y0, [result] + list(others), tr, n_boot_curves, True)
    zc = cube.zc_au(y0)
    Vz, Vr = result.velocities([zc])
    ax.text(0.02, 0.02, f"at z$_c$ = {zc:.0f} au: V$_z$ = {Vz[0]:.2f}, V$_r$ = {Vr[0]:.2f} km/s",
            transform=ax.transAxes, fontsize="x-small", bbox=dict(fc="w", alpha=0.7, lw=0))
    pc = cube.config.plot
    ax.set_xlim(*(pc.xlim or (cube.x[0], cube.x[-1])))
    if pc.vlim:
        ax.set_ylim(*pc.vlim)
    ax.legend(fontsize="x-small", loc="lower right")
    ax.set_title(f"{cube.config.source.name}  {cut_label(y0, zc)}", fontsize="small")
    ax.set_xlabel("distance to the wind axis (arcsec)")
    ax.set_ylabel(r"velocity (km s$^{-1}$)")
    fig.tight_layout()
    return finish(fig, save, cube.config.output.dpi)


#: Readable axis labels for common parameter names (others are shown as named).
PARAM_LABELS = {
    "R0_au": r"$R_0$ (au)", "xoff_au": r"axis offset $x_0$ (au)", "xslope": r"axis tilt $\xi$",
    "Vz": r"$V_z$ (km s$^{-1}$)", "Vz_b": r"$V_{z,b}$ (km s$^{-1}$)", "Vr": r"$V_r$ (km s$^{-1}$)",
    "j": r"$j$ (au km s$^{-1}$)", "qv": r"$q_v$", "q1": r"$q_1$", "q2": r"$q_2$",
    "z_b": r"$z_b$ (au)", "qj": r"$q_j$",
}


def plot_bootstrap(result, save=None):
    """Histograms of the bootstrap samples of the free parameters.

    Returns
    -------
    Figure
    """
    if result.boot is None or not len(result.boot):
        raise ValueError("this result has no bootstrap samples (fit with n_boot > 0)")
    p = plt()
    show = result.free
    fig, axs = p.subplots(1, len(show), figsize=(2.3 * len(show), 2.3))
    for ax, k in zip(np.atleast_1d(axs), show):
        vals = result.boot[:, result.names.index(k)]
        best = result.params[k]
        if k == "tan_theta":
            vals, best, k = np.degrees(np.arctan(vals)), np.degrees(np.arctan(best)), "θ (deg)"
        ax.hist(vals, bins=20, color="0.7")
        ax.axvline(best, color="r")
        ax.set_xlabel(PARAM_LABELS.get(k, k))
        ax.set_yticks([])
        ax.xaxis.set_major_locator(MaxNLocator(3))
        ax.tick_params(axis="x", labelsize="small")
    meth = result.settings.get("error_method", "bootstrap")
    fig.suptitle(f"{meth} over cuts (N = {len(result.boot)})", fontsize="small")
    fig.tight_layout()
    return finish(fig, save)


def plot_profiles(result, others=(), z=None, save=None):
    """The fitted wind as a function of height, with bootstrap (16-84%) bands.

    Panels: (1) de-projected velocities Vz, Vr, Vp; (2) what the data measure
    directly, the ring-center velocity Vcent - Vsys and Vr sin i; (3) the
    direction of the flow, atan(Vr/Vz), with the cone half-angle; (4) the radius
    R(z) and the rotation velocity V_phi. Thin gray lines: heights of the cuts;
    dotted line: break height of a broken power law. Other results are dashed.

    Returns
    -------
    Figure
    """
    p = plt()
    g = result.geometry
    zc = np.array([g.zc_au(y) for y in result.traces]) if result.traces else np.array([50, 500])
    if z is None:
        z = np.linspace(max(1.0, 0.3 * zc.min()), 1.1 * zc.max(), 200)
    fig, axs = p.subplots(1, 4, figsize=(17, 3.8))
    names = {"Vz": r"$V_z$", "Vr": r"$V_r$", "Vp": r"$V_p$"}
    for k, r in enumerate([result] + list(others)):
        P = r.profiles(z)
        main = k == 0
        ls, lw = ("-", 1.6) if main else ("--", 1.0)
        tag = "" if main else f" ({r.label})"
        for q, col in (("Vz", "k"), ("Vr", "tab:green"), ("Vp", "tab:purple")):
            if main:
                axs[0].fill_between(z, P[q + "_lo"], P[q + "_hi"], color=col, alpha=0.2, lw=0)
            axs[0].plot(z, P[q], color=col, ls=ls, lw=lw,
                        label=names[q] + tag if (main or q == "Vp") else None)
        axs[1].plot(z, -g.s * P["Vz"] * g.cosi, "k", ls=ls, lw=lw,
                    label=r"$V_{\rm cent}-V_{\rm sys}$" + tag if main else None)
        axs[1].plot(z, P["Vr"] * g.sini, color="tab:green", ls=ls, lw=lw,
                    label=r"$V_r\sin i$" + tag if main else None)
        if main:
            axs[2].fill_between(z, P["flow_angle_deg_lo"], P["flow_angle_deg_hi"], color="k",
                                alpha=0.2, lw=0)
            axs[3].fill_between(z, P["R_lo"], P["R_hi"], color="k", alpha=0.2, lw=0)
        axs[2].plot(z, P["flow_angle_deg"], "k", ls=ls, lw=lw, label="flow direction" + tag)
        axs[3].plot(z, P["R"], "k", ls=ls, lw=lw, label="R" + tag)
        if main:
            ax2 = axs[3].twinx()
            ax2.fill_between(z, P["Vphi_lo"], P["Vphi_hi"], color="tab:blue", alpha=0.2, lw=0)
            ax2.plot(z, P["Vphi"], color="tab:blue", lw=lw)
            ax2.set_ylabel(r"$V_\phi$ (km s$^{-1}$)", color="tab:blue")
            ax2.tick_params(axis="y", colors="tab:blue")
    th = result.theta_deg[0]
    axs[2].axhline(th, color="tab:orange", lw=1.2, ls="--", label="cone half-angle")
    lo, hi = axs[2].get_ylim()               # at least +-5 deg, so flat lines stay readable
    mid = 0.5 * (lo + hi)
    axs[2].set_ylim(min(lo, mid - 5), max(hi, mid + 5))
    if result.model.radial == "cone":
        axs[2].text(0.03, 0.04, "flow tied to the cone", transform=axs[2].transAxes, fontsize="x-small",
                    color="0.35")
    if "z_b" in result.params:
        for ax in axs:
            ax.axvline(result.params["z_b"], color="0.5", ls=":", lw=0.8)
    for ax, t, yl in zip(axs, ("de-projected velocities", "measured (projected)",
                                "flow direction", r"radius and rotation"),
                         (r"km s$^{-1}$", r"km s$^{-1}$", "deg", "R (au)")):
        ax.set_xlabel("height z (au)")
        ax.set_title(t, fontsize="small")
        ax.set_ylabel(yl)
        ax.legend(fontsize="x-small", loc="best", framealpha=0.8)
        for zz in zc:
            ax.axvline(zz, color="0.88", lw=0.5, zorder=0)
    fig.tight_layout()
    return finish(fig, save)


def plot_ring_fits(ringset, save=None):
    """Per-cut ring fits vs height: R(z) with the cone, Vcent - Vsys and Vr sin i,
    Vphi sin i, and j (Louvet+2018 Figs. 12 and 14).

    Returns
    -------
    Figure
    """
    p = plt()
    t = ringset.table()
    g = t["good"].astype(bool)
    fig, axs = p.subplots(1, 4, figsize=(16, 3.6))
    z = t["z_au"]
    for ax, (k, lab) in zip(axs, (("R_au", "R (au)"), ("vcent_rel", r"V$_{cent}$-V$_{sys}$"),
                                 ("Vphi_sini", r"V$_\phi \sin i$"), ("j", "j (au km/s)"))):
        ax.errorbar(z[g], t[k][g], yerr=t[k + "_err"][g], fmt="ko", ms=3)
        ax.errorbar(z[~g], t[k][~g], yerr=t[k + "_err"][~g], fmt="o", mfc="none", color="0.6",
                    ms=3)
        ax.set_xlabel("z (au)")
        ax.set_ylabel(lab)
    axs[1].errorbar(z[g], t["Vr_sini"][g], yerr=t["Vr_sini_err"][g], fmt="gs", ms=3,
                    label=r"V$_r \sin i$")
    axs[1].legend(fontsize="x-small", loc="center right")
    if ringset.cone:
        c = ringset.cone
        zz = np.linspace(0, z.max() * 1.05, 50)
        axs[0].plot(zz, c["R0"] + c["tan"] * zz, "r-",
                    label=f"θ = {c['theta']:.1f} ± {c['theta_err']:.1f}°")
        axs[0].legend(fontsize="x-small")
    fig.tight_layout()
    return finish(fig, save)
