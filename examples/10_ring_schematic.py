#!/usr/bin/env python
"""
10 -- Schematic of the fitted wind: rings, velocities and the observer's view.

    python 10_ring_schematic.py [--config configs/elias20.toml] [--nring 8]
           [--zmax 450] [--azim 35] [--elev 18]

  (a) 3-D view: rings colored by V_phi, poloidal velocity arrows, rotation of the
      wind and of the disk, the cone;
  (b) the model as the observer sees it (projected at the disk inclination);
  (c) meridional cut: cone walls, rings at the cut heights, poloidal velocities.

The rings are evaluated from the saved fit <output.dir>/broken-power-law.json
(run 02_elias20_analysis.py first) at the heights of the fitted cuts, or at
--nring evenly spaced heights. Writes <output.dir>/ring_schematic.<fig_format>.
"""
import argparse
import os

import numpy as np

import _schematic_helpers as wk


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="configs/elias20.toml")
    ap.add_argument("--nring", type=int, default=None)
    ap.add_argument("--zmax", type=float, default=None)
    ap.add_argument("--azim", type=float, default=35.0)
    ap.add_argument("--elev", type=float, default=18.0)
    ap.add_argument("--model", choices=wk.VELOCITY_MODELS, default=None,
                    help="velocity law (default: VELOCITY_MODEL in the config)")
    args = ap.parse_args()
    cfg = wk.load_config(args.config)
    plt = wk.setup_matplotlib()
    from matplotlib import cm, colors

    gpath = wk.result_path(cfg)
    if not os.path.exists(gpath):
        raise SystemExit(f"[schematic] {gpath} not found: run the fit first (diskwind fit)")
    G = wk.geometry_meta(cfg)      # geometry metadata (s, i, d, cuts)
    VM = wk.VelocityModel(cfg, args.model or cfg.VELOCITY_MODEL)
    gp = dict(zip(VM.names, VM.p))
    ge = dict(zip(VM.names, np.nanstd(VM.draws, axis=0)))
    gp.setdefault("qj", 0.0)
    ge.setdefault("qj", 0.0)

    # geometry stored with the saved fit
    inc_deg, dist = G["inc_deg"], G["dist_pc"]
    s, dsense = G["s"], G["disk_sense"]
    inc = np.radians(inc_deg)
    si, ci = np.sin(inc), np.cos(inc)

    R0, T, xo, xsl = gp["R0_au"], gp["tan_theta"], gp["xoff_au"], gp["xslope"]
    j, qj = gp["j"], gp["qj"]
    zref = wk.GLOBAL_ZREF_AU
    Vz0, Vr0 = (float(v[0]) for v in VM.vel(np.array([zref])))

    ys = G.get("y_offsets") or list(np.atleast_1d(cfg.FIT_Y_OFFSETS_ARCSEC))
    z_slits = np.array(ys) / si * dist
    zmax = args.zmax or 1.05 * z_slits.max()
    if args.nring:
        z = np.linspace(zmax / args.nring, zmax, args.nring)
    else:
        z = z_slits[z_slits <= zmax]
    R = R0 + T * z
    xc = xo + xsl * z
    Vz, Vr = VM.vel(z)
    Vphi = VM.jz(z) / R

    rdisk = min(cfg.DISK_ROUT_AU, 1.2 * R.max() + 40)
    vlim = max(np.abs(Vphi).max(), 0.1)
    norm = colors.TwoSlopeNorm(vmin=-vlim, vcenter=0.0, vmax=vlim)
    cmap = cm.coolwarm
    phi = np.linspace(0, 2 * np.pi, 181)
    zz = np.linspace(0, zmax, 60)

    fig = plt.figure(figsize=(19, 7.2))

    # ------------------------------------------------ (a) 3D view
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    PP, RR = np.meshgrid(phi, np.linspace(0, rdisk, 2))
    ax.plot_surface(RR * np.cos(PP), RR * np.sin(PP), 0 * RR, color="0.75", alpha=0.35,
                    linewidth=0, shade=False)
    ax.plot(rdisk * np.cos(phi), rdisk * np.sin(phi), 0, color="0.4", lw=0.8)
    Pc, Zc = np.meshgrid(phi, zz)
    Rc = R0 + Zc * T
    ax.plot_surface(xo + xsl * Zc + Rc * np.cos(Pc), Rc * np.sin(Pc), Zc, color="tab:orange",
                    alpha=0.08, linewidth=0, shade=False)
    ax.plot(xo + xsl * zz, 0 * zz, zz, color="k", lw=1)
    dz = np.median(np.diff(z)) if len(z) > 1 else 0.2 * zmax
    vscale = 0.9 * dz / max(np.hypot(Vz, Vr).max(), 1e-3)
    for zi, Ri, xi, vz, vr, vp in zip(z, R, xc, Vz, Vr, Vphi):
        col = cmap(norm(vp))
        ax.plot(xi + Ri * np.cos(phi), Ri * np.sin(phi), np.full_like(phi, zi), color=col, lw=2)
        for p0 in (0, np.pi / 2, np.pi, 1.5 * np.pi):
            ax.quiver(xi + Ri * np.cos(p0), Ri * np.sin(p0), zi,
                      vr * np.cos(p0) * vscale, vr * np.sin(p0) * vscale, vz * vscale,
                      color="k", lw=0.6, arrow_length_ratio=0.25)
        p0 = 0.75 * np.pi
        sg = np.sign(vp) if vp != 0 else 1
        ax.quiver(xi + Ri * np.cos(p0), Ri * np.sin(p0), zi,
                  -sg * np.sin(p0) * 0.3 * Ri, sg * np.cos(p0) * 0.3 * Ri, 0,
                  color=col, lw=1.2, arrow_length_ratio=0.35)
    if dsense:
        p0 = 0.25 * np.pi
        ax.quiver(0.9 * rdisk * np.cos(p0), 0.9 * rdisk * np.sin(p0), 0,
                  -dsense * np.sin(p0) * 0.3 * rdisk, dsense * np.cos(p0) * 0.3 * rdisk, 0,
                  color="0.3", lw=1.5, arrow_length_ratio=0.3)
        ax.text(0.9 * rdisk * np.cos(p0), 0.9 * rdisk * np.sin(p0), 0, " disk", color="0.3",
                fontsize=7)
    lim = max(rdisk, R.max() + np.abs(xc).max())
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(0, zmax * 1.05)
    ax.set_box_aspect((1, 1, zmax * 1.05 / (2 * lim)))
    ax.view_init(elev=args.elev, azim=args.azim)
    ax.set_xlabel("x (au)")
    ax.set_ylabel("y (au)")
    ax.set_zlabel("z (au)")
    from matplotlib.ticker import MaxNLocator
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.set_major_locator(MaxNLocator(3))
    ax.tick_params(labelsize="small")
    ax.set_title("(a) model rings (color = V$_\\phi$ = j/R), velocity arrows", fontsize="small")
    cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.55, pad=0.12)
    cb.set_label(r"$V_\phi$ (km/s; > 0: counter-clockwise about +z)", fontsize="small")
    cb.ax.tick_params(labelsize="small")

    # ------------------------------------------------ (b) observer's view
    ax = fig.add_subplot(1, 3, 2)
    Xd, Yd = rdisk * np.cos(phi), s * rdisk * np.sin(phi) * ci
    ax.fill(Xd, Yd, color="0.85", zorder=0)
    ax.plot(Xd, Yd, color="0.5", lw=0.8)
    for zi, Ri, xi, vp in zip(z, R, xc, Vphi):
        X = xi + Ri * np.cos(phi)
        Y = s * Ri * np.sin(phi) * ci + zi * si
        near = np.sin(phi) < 0  # depth toward observer W = -R sin(phi) sin i > 0
        col = cmap(norm(vp))
        ax.plot(np.where(near, X, np.nan), np.where(near, Y, np.nan), color=col, lw=1.8)
        ax.plot(np.where(~near, X, np.nan), np.where(~near, Y, np.nan), color=col, lw=0.9, ls="--")
    ax.plot(xo + xsl * zz, zz * si, color="k", lw=0.8)
    for sg in (+1, -1):  # projected cone edges (limb at phi = 0, pi)
        ax.plot(xo + xsl * zz + sg * (R0 + T * zz), zz * si, color="tab:orange", lw=0.8, alpha=0.7)
    ax.set_aspect("equal")
    ax.set_xlabel("x (au, along disk major axis)")
    ax.set_ylabel("sky offset along wind axis (au)")
    ax.set_title(f"(b) as seen by the observer (i = {inc_deg}°, lobe "
                 f"{'approaching' if s > 0 else 'receding'})\nsolid = near side, dashed = far side",
                 fontsize="small")

    # ------------------------------------------------ (c) meridional cut
    # Drawn relative to the fitted wind axis, x - x_c(z), so the cone is symmetric
    # (the axis offset x0 and tilt xi are given in the text below the panel).
    ax = fig.add_subplot(1, 3, 3)
    Rz = R0 + T * zz
    ax.fill_betweenx(zz, -Rz, Rz, color="tab:orange", alpha=0.08)
    for sg in (+1, -1):
        ax.plot(sg * Rz, zz, color="tab:orange", lw=1.4)
    z_apex = -R0 / T if T > 0 else -np.inf            # where the walls would meet
    zlo = max(z_apex, -0.25 * zmax)
    zb = np.linspace(zlo, 0, 50)
    for sg in (+1, -1):                                 # extension below the disk
        ax.plot(sg * (R0 + T * zb), zb, color="tab:orange", lw=0.9, ls="--")
    ax.plot([0, 0], [zlo, zz.max()], color="tab:orange", lw=0.6, ls=":")
    if np.isfinite(z_apex) and z_apex >= -0.25 * zmax:
        ax.plot(0, z_apex, "x", color="tab:orange", ms=5)
        ax.text(0.04 * zmax, z_apex, f"apex, z = {z_apex:.0f} au", fontsize="x-small", va="center",
                color="tab:orange")
    ax.plot(np.r_[R, -R], np.r_[z, z], "o", color="k", ms=3.5)
    for zi, Ri, vz, vr in zip(z, R, Vz, Vr):
        for sg in (+1, -1):
            ax.annotate("", xy=(sg * (Ri + vr * vscale), zi + vz * vscale),
                        xytext=(sg * Ri, zi),
                        arrowprops=dict(arrowstyle="->", color="tab:blue", lw=0.9))
    lim = 1.15 * R.max()
    ax.fill_between([-lim, lim], -0.02 * zmax, 0, color="0.8", zorder=0)
    ax.annotate("", xy=(-R0, 0.03 * zmax), xytext=(R0, 0.03 * zmax),
                arrowprops=dict(arrowstyle="<->", color="0.3", lw=0.8))
    ax.text(0, 0.05 * zmax, f"2R$_0$ = {2 * R0:.0f} au", fontsize="x-small", ha="center", va="bottom",
            color="0.3", bbox=dict(fc="w", ec="none", pad=0.5))
    vref = 5.0
    x0r = -0.5 * vref * vscale                            # inside the cone, at the top
    ax.annotate("", xy=(x0r + vref * vscale, 0.97 * zmax), xytext=(x0r, 0.97 * zmax),
                arrowprops=dict(arrowstyle="->", color="tab:blue", lw=0.9))
    ax.text(x0r, 0.92 * zmax, f"{vref:g} km/s", color="tab:blue", fontsize="x-small", va="top")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(zlo - 0.03 * zmax, zmax * 1.05)
    ax.set_aspect("equal")
    ax.set_xlabel(r"distance from the fitted wind axis, $x - x_c(z)$ (au)")
    ax.set_ylabel("height z above the disk (au)")
    th = np.degrees(np.arctan(T))
    eth = np.degrees(ge["tan_theta"] / (1 + T ** 2))
    if VM.model == "single":
        vtxt = (f"V$_z$ = {Vz0:.2f} ± {ge['Vz']:.2f},  V$_r$ = {Vr0:.2f} ± {ge['Vr']:.2f} km/s"
                f"  (at z = {zref:.0f} au; ~ z$^{{{gp['qv']:.2f}}}$)\n")
    else:
        vtxt = (f"V$_z$ = {gp['Vz_b']:.2f} ± {ge['Vz_b']:.2f} km/s at z$_b$ = {gp['z_b']:.0f} au;"
                f"  ~ z$^{{{gp['q1']:.2f}}}$ below, z$^{{{gp['q2']:.2f}}}$ above\n")
    txt = (f"θ = {th:.1f} ± {eth:.1f}°   R$_0$ = {R0:.0f} ± {ge['R0_au']:.0f} au\n" + vtxt +
           f"j = R V$_\\phi$ = {j:.0f} ± {ge['j']:.0f} au km/s at z = {zref:.0f} au"
           + (f", ~ z$^{{{qj:.2f} \\pm {ge['qj']:.2f}}}$" if qj != 0 or ge["qj"] != 0 else "")
           + f"\naxis: x$_c$(z) = x$_0$ + ξz, x$_0$ = {xo:.0f} ± {ge['xoff_au']:.0f} au,"
             f" ξ = {xsl:.3f} ± {ge['xslope']:.3f}")
    ax.text(0.5, -0.2, txt, transform=ax.transAxes, fontsize="x-small", ha="center", va="top")
    ax.set_title("(c) meridional cut: global cone, rings and poloidal velocity", fontsize="small")

    rot = ("" if not dsense else
           f"; wind rotates {'WITH' if j * dsense > 0 else 'AGAINST'} the disk (disk sense {dsense:+d})")
    fig.suptitle(f"{cfg.SOURCE_NAME}: ring model from the global fit "
                 f"({VM.name}){rot}", fontsize="medium", y=1.0)
    fig.tight_layout()
    from diskwind.plotting._style import thin_ticks
    thin_ticks(fig, 4)
    out = wk.figpath(cfg, "ring_schematic")
    fig.savefig(out, dpi=cfg.FIG_DPI)
    print(f"[schematic] wrote {out}")


if __name__ == "__main__":
    main()
