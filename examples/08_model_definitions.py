#!/usr/bin/env python
"""
08 -- Schematic of the model definitions (four panels).

    python 08_model_definitions.py [--config configs/elias20.toml] [--y0 0.7]
           [--illustrative]

  (a) the sky as observed: disk and wind position angles, the axes X and Y,
      projected rings and one transverse cut;
  (b) side view: inclination i, lobe orientation s, R0, theta, z, R and the two
      heights at which one cut meets the cone;
  (c) one ring seen from the tip of the lobe: phi, R and the velocity components;
  (d) the resulting curve in the PV diagram of the cut, with the fitted distance d_k.

The geometry comes from the configuration and the cone and velocities from the
saved fit <output.dir>/broken-power-law.json (run 02_elias20_analysis.py first),
or illustrative values with --illustrative.
Writes <output.dir>/model_definitions.<fig_format>.
"""
import argparse
import os

import numpy as np

import _schematic_helpers as wk

ZREF = wk.GLOBAL_ZREF_AU
RED, GRAY, BLUE, PURPLE, ORANGE = "#c0392b", "#95a5a6", "#2471a3", "#8e44ad", "#00897b"


def arrow(ax, p0, p1, color="k", lw=1.3):
    ax.annotate("", xy=p1, xytext=p0, arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                                                      mutation_scale=9, shrinkA=0, shrinkB=0))


def arc(ax, c, r, a0, a1, **kw):
    t = np.linspace(min(a0, a1), max(a0, a1), 60)
    ax.plot(c[0] + r * np.cos(t), c[1] + r * np.sin(t), color="k", lw=0.8, **kw)
    return 0.5 * (a0 + a1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="configs/elias20.toml")
    ap.add_argument("--y0", type=float, default=0.7, help="cut shown in (a) and (c), arcsec")
    ap.add_argument("--illustrative", action="store_true",
                    help="ignore the saved fit and use illustrative cone values")
    args = ap.parse_args()
    cfg = wk.load_config(args.config)
    plt = wk.setup_matplotlib()
    # plot style (classic, 14 pt, thick axes)
    plt.style.use("classic")
    plt.rcParams["font.size"] = 14
    plt.rcParams["axes.linewidth"] = 2.0
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.weight"] = "500"
    plt.rcParams["mathtext.bf"] = "serif:normal"

    inc = np.radians(float(cfg.INCLINATION_DEG))
    si, ci = np.sin(inc), np.cos(inc)
    dist = float(cfg.DISTANCE_PC)
    s = {True: 1, False: -1}.get(cfg.WIND_APPROACHING, -1)
    dsense = cfg.DISK_ROTATION_SENSE or 1
    vsys = cfg.VSYS_KMS if cfg.VSYS_KMS is not None else 0.0
    P = dict(R0_au=15.0, tan_theta=np.tan(np.radians(30)), xoff_au=0.0, xslope=0.0,
             Vz=5.0, Vr=5.0 * np.tan(np.radians(30)), j=40.0, qv=0.3, qj=0.0)
    src = "illustrative"
    gpath = wk.result_path(cfg)
    VM = None
    if os.path.exists(gpath) and not args.illustrative:
        G = wk.geometry_meta(cfg)
        try:
            VM = wk.VelocityModel(cfg, cfg.VELOCITY_MODEL)
        except FileNotFoundError:
            VM = wk.VelocityModel(cfg, "single")
        P = dict(zip(VM.names, VM.p))
        s, dsense = int(G.get("s", s)), int(G.get("disk_sense", dsense))
        vsys = float(G.get("vsys", vsys))
        src = VM.name
    T = P["tan_theta"]
    y0 = args.y0
    Ycut = y0 * dist

    def Rz(z):
        return P["R0_au"] + z * T

    def cut_z(phi):
        sp = np.sin(phi)
        den = si + s * T * sp * ci
        z = (Ycut - s * sp * ci * P["R0_au"]) / den
        return np.where((den > 0) & (z > 0), z, np.nan)

    fig, axs = plt.subplots(2, 2, figsize=(11.85, 10.8),
                            gridspec_kw=dict(hspace=0.32, wspace=0.28))
    axs = axs.ravel()
    titles = {}

    # ---------------- (a) the sky as observed ----------------
    ax = axs[0]
    pa_w, pa_d = float(cfg.WIND_PA_DEG), float(cfg.DISK_PA_DEG)
    Yh = np.array([np.sin(np.radians(pa_w)), np.cos(np.radians(pa_w))])         # (E, N)
    Xh = np.array([np.sin(np.radians(pa_w - 90)), np.cos(np.radians(pa_w - 90))])

    def to_sky(Xa, Ya):           # arcsec in the wind frame -> (dE, dN)
        return Xa * Xh[0] + Ya * Yh[0], Xa * Xh[1] + Ya * Yh[1]
    ph = np.linspace(0, 2 * np.pi, 721)
    nr = np.sin(ph) < 0           # near half of a ring
    Rd = 0.9                      # disk radius (arcsec, illustrative)
    e, n = to_sky(Rd * np.cos(ph), s * Rd * np.sin(ph) * ci)
    ax.fill(e, n, color=GRAY, alpha=0.3, lw=0)
    ax.plot(e, n, color="0.5", lw=0.8)
    for zr in np.array([0.35, 0.7, 1.0, 1.2]) * y0 / si * dist:
        Rr = Rz(zr)
        Xa = (P["xoff_au"] + P["xslope"] * zr + Rr * np.cos(ph)) / dist
        Ya = (zr * si + s * Rr * np.sin(ph) * ci) / dist
        e, n = to_sky(Xa, Ya)
        ax.plot(np.where(nr, e, np.nan), np.where(nr, n, np.nan), color=RED, lw=1.0)
        ax.plot(np.where(~nr, e, np.nan), np.where(~nr, n, np.nan), color=RED, lw=0.8,
                ls=(0, (3, 2)), alpha=0.75)
    xh = 1.6
    w = cfg.SLIT_WIDTH_ARCSEC or 0.1
    e, n = to_sky(np.array([-xh, xh, xh, -xh]), np.array([y0 - w / 2] * 2 + [y0 + w / 2] * 2))
    ax.fill(e, n, color=BLUE, alpha=0.25, lw=0)
    e, n = to_sky(np.array([-xh, xh]), np.array([y0, y0]))
    ax.plot(e, n, color=BLUE, lw=1.0)
    rot = (pa_w % 180) - (180 if pa_w % 180 > 90 else 0)
    e, n = to_sky(np.array([-0.85 * xh]), np.array([y0 - w / 2 - 0.3]))
    ax.text(e[0], n[0], r"cut at $Y=y_0$", color=BLUE, fontsize=12, ha="center", va="center",
            rotation=rot, rotation_mode="anchor")
    L = 1.45
    e, n = to_sky(np.array([0, 0]), np.array([0, L]))
    arrow(ax, (e[0], n[0]), (e[1], n[1]), color=RED, lw=1.6)
    ax.text(e[1] + 0.12, n[1] - 0.12, r"$+Y$ (wind)", color=RED, fontsize=12.2, va="top", ha="center",
            bbox=dict(fc="w", ec="none", alpha=0.85, pad=0.6))
    e, n = to_sky(np.array([0, 1.2]), np.array([0, 0]))
    arrow(ax, (e[0], n[0]), (e[1], n[1]), color="#27ae60", lw=1.6)
    ax.text(e[1], n[1], r"  $+X$", color="#27ae60", fontsize=12.2, va="center")
    pdr = np.radians(pa_d)
    ax.plot([-1.5 * np.sin(pdr), 1.5 * np.sin(pdr)], [-1.5 * np.cos(pdr), 1.5 * np.cos(pdr)],
            ls="--", color="0.5", lw=0.8)
    ax.plot([0, 0], [0, 1.6], color="0.5", lw=0.6, ls=":")
    for pa, r, rl, col, lab in ((pa_d % 180, 0.3, 1.05, "0.35", r"PA$_d$"),
                                (pa_w, 0.5, 1.15, RED, r"PA$_w$")):
        tt_ = np.radians(np.linspace(0, pa, 80))
        ax.plot(r * np.sin(tt_), r * np.cos(tt_), color=col, lw=0.9)
        tm = np.radians(pa / 2)
        ax.plot([r * np.sin(tm), (rl - 0.22) * np.sin(tm)], [r * np.cos(tm), (rl - 0.22) * np.cos(tm)],
                color=col, lw=0.5)
        ax.text(rl * np.sin(tm), rl * np.cos(tm), lab, color=col, fontsize=13, ha="center",
                va="center")
    ax.plot(0, 0, "k+", ms=7)
    cxy = (1.0, 1.1)
    arrow(ax, cxy, (cxy[0], cxy[1] + 0.45), lw=1)
    arrow(ax, cxy, (cxy[0] + 0.45, cxy[1]), lw=1)
    ax.text(cxy[0], cxy[1] + 0.5, "N", ha="center", fontsize=10.8)
    ax.text(cxy[0] + 0.52, cxy[1], "E", va="center", ha="right", fontsize=10.8)
    ax.set_xlim(1.8, -1.8)
    ax.set_ylim(-1.8, 1.8)
    ax.set_xticks([1, 0, -1])
    ax.set_yticks([-1, 0, 1])
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\Delta\alpha\,\cos\delta$ (arcsec)")
    ax.set_ylabel(r"$\Delta\delta$ (arcsec)")
    titles[0] = "(a) the sky, as observed"

    # ---------------- (b) side view ----------------
    ax = axs[1]
    zhat = np.array([si, s * ci])        # (sky Y, toward observer)
    yw = np.array([s * ci, -si])         # phi = 90 deg direction
    z90, z270 = float(cut_z(np.array([np.pi / 2]))[0]), float(cut_z(np.array([1.5 * np.pi]))[0])
    Zmax = 1.12 * np.nanmax([z90, z270, Ycut / si])
    zz = np.linspace(0, Zmax, 40)
    walls = []
    for sg in (1, -1):
        p = np.outer(zz, zhat) + sg * np.outer(Rz(zz), yw)
        walls.append(p)
        ax.plot(p[:, 0], p[:, 1], color=RED, lw=1.6)
    ax.plot([0, Zmax * zhat[0]], [0, Zmax * zhat[1]], color=RED, lw=0.8, ls="-.")
    ax.text(*(Zmax * zhat), r"  $z$", color=RED, fontsize=14.9, va="center")
    allp = np.vstack(walls)
    Dm = 0.4 * Zmax
    ax.plot([-Dm * yw[0], Dm * yw[0]], [-Dm * yw[1], Dm * yw[1]], color=GRAY, lw=4, alpha=0.6,
            solid_capstyle="butt")
    lo_end = -Dm * yw if (-yw)[1] < 0 else Dm * yw
    dang = np.degrees(np.arctan2(yw[1], yw[0]))
    dang = dang - 180 if dang > 90 else (dang + 180 if dang < -90 else dang)
    dpos = 0.5 * lo_end
    nrm_ = np.array([-yw[1], yw[0]])
    nrm_ = nrm_ if nrm_[1] > 0 else -nrm_
    dpos = dpos + 0.09 * Zmax * nrm_
    ax.text(dpos[0], dpos[1], "disk", color="0.45",
            fontsize=11.5, ha="center", va="center", rotation=dang, rotation_mode="anchor")
    # theta: measured from the axis, between +z and a line parallel to the near wall
    za = 0.42 * Zmax
    A = za * zhat                                    # point on the axis
    wall = zhat - T * yw                             # direction of the near wall
    wall = wall / np.hypot(*wall)
    lw_ = 0.50 * Zmax
    ax.plot([A[0], A[0] + lw_ * wall[0]], [A[1], A[1] + lw_ * wall[1]], color=RED, lw=0.9,
            ls=(0, (4, 2)))
    rth = 0.42 * Zmax
    tm = arc(ax, A, rth, np.arctan2(zhat[1], zhat[0]), np.arctan2(wall[1], wall[0]))
    ax.text(A[0] + 1.18 * rth * np.cos(tm), A[1] + 1.18 * rth * np.sin(tm), r"$\theta$",
            fontsize=14.9, ha="center", va="center")
    # R0: from the axis to the base of the near wall, along the midplane
    b0 = -P["R0_au"] * yw
    off = -0.09 * Zmax * zhat                        # offset from the disk band, outside the cone
    ax.annotate("", xy=b0 + off, xytext=off,
                arrowprops=dict(arrowstyle="<->,head_length=0.25,head_width=0.15", color="k", lw=0.8,
                                shrinkA=0, shrinkB=0))
    for e_ in (np.zeros(2), b0):                     # short ticks back to the disk
        ax.plot([e_[0], e_[0] + 1.1 * off[0]], [e_[1], e_[1] + 1.1 * off[1]], color="k", lw=0.5)
    ax.text(*(0.5 * b0 + 1.9 * off), r"$R_0$", fontsize=13.5, ha="center", va="center")
    # one ring: radius R at height z1
    z1 = 0.6 * Zmax
    c1 = z1 * zhat
    p90 = c1 + Rz(z1) * yw
    ax.plot([c1[0], p90[0]], [c1[1], p90[1]], color="k", lw=1.0)
    ax.text(*(0.5 * (c1 + p90)), r" $R(z)$", fontsize=13.5, va="top")
    # inclination: between +z and the line of sight through the star
    los = np.array([0.0, np.sign(zhat[1]) or 1.0])
    ri = 0.22 * Zmax
    ax.axvline(0, color="k", lw=0.8, ls="--", zorder=1)   # line of sight through the star
    tm = arc(ax, (0, 0), ri, np.arctan2(zhat[1], zhat[0]), np.arctan2(los[1], los[0]))
    ax.text(1.22 * ri * np.cos(tm), 1.22 * ri * np.sin(tm), r"$i$", fontsize=14.9, ha="center",
            va="center")

    # observer
    top = max(allp[:, 1].max(), 0) + 0.12 * Zmax
    xo = allp[:, 0].max() + 0.12 * Zmax
    arrow(ax, (xo, top - 0.24 * Zmax), (xo, top + 0.08 * Zmax), lw=1.6)
    ax.text(xo - 0.03 * Zmax, top - 0.08 * Zmax, "to\nobserver", fontsize=11.5, ha="right",
            va="center")
    # the cut and its two intersections
    ax.axvline(Ycut, color=BLUE, lw=1.2)
    ax.text(Ycut, top + 0.1 * Zmax, r" $Y=y_0$", color=BLUE, fontsize=11.5, va="top")
    ax.text(-0.02 * Zmax, allp[:, 1].min() + 0.02 * Zmax, "line of sight to the star",
            fontsize=10.1, ha="right", va="bottom", rotation=90)
    for zc_, sg, lab in ((z90, 1, "far"), (z270, -1, "near")):
        if np.isfinite(zc_):
            pt = zc_ * zhat + sg * Rz(zc_) * yw
            ax.plot(*pt, "o", color=BLUE, ms=5, zorder=5)
            ax.text(pt[0], pt[1], f"  {lab}", color=BLUE, fontsize=11.5,
                    va="center", bbox=dict(fc="w", ec="none", alpha=0.8, pad=0.4))
    ax.set_xlim(allp[:, 0].min() - 0.28 * Zmax, allp[:, 0].max() + 0.3 * Zmax)
    ax.set_ylim(allp[:, 1].min() - 0.07 * Zmax, top + 0.12 * Zmax)
    ax.set_aspect("equal")
    ax.set_xlabel(r"sky offset along the projected axis, $Y$ (au)")
    ax.set_ylabel("distance toward the observer (au)")
    titles[1] = (f"(b) side view (lobe {'approaching' if s > 0 else 'receding'}, $s$ = {s:+d})"
                 .replace("-", "\u2212"))

    # ---------------- (c) ring from the lobe tip ----------------
    ax = axs[2]
    xc0 = 0.32                       # ring center offset from the axis (exaggerated)
    t = np.linspace(0, 2 * np.pi, 361)
    ax.plot(xc0 + np.cos(t), np.sin(t), color=RED, lw=1.8)
    ax.fill_between(xc0 + np.cos(t[180:]), np.sin(t[180:]), 0, color=RED, alpha=0.07, lw=0)
    ax.plot([-1.45, 1.72], [0, 0], color="0.6", lw=0.6)
    ax.plot([0, 0], [-1.35, 1.4], color="0.6", lw=0.6)
    arrow(ax, (1.47, 0), (1.72, 0), color="0.4", lw=0.9)
    arrow(ax, (0, 1.15), (0, 1.4), color="0.4", lw=0.9)
    ax.text(1.75, 0, r"$+X$", fontsize=12.2, va="center")
    ax.text(0.05, 1.42, r"$+y_w$", fontsize=12.2, va="bottom", ha="left")
    ax.plot(0, 0, "o", mfc="w", mec="k", ms=11)
    ax.plot(0, 0, "k.", ms=4)
    ax.annotate(r"axis: $+z$ and $V_z$" "\nout of the page", xy=(-0.06, 0.06),
                xytext=(-0.62, 0.62), fontsize=10.1, ha="right", va="bottom",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.3", shrinkA=0, shrinkB=4))
    ax.plot(xc0, 0, "+", color=RED, ms=9, mew=1.4)
    # axis offset of the ring center
    ax.plot([xc0, xc0], [0, -0.32], color="k", lw=0.6, ls=":")
    ax.annotate("", xy=(xc0, -0.27), xytext=(0, -0.27),
                arrowprops=dict(arrowstyle="<->", lw=0.8, shrinkA=0, shrinkB=0))
    ax.text(-0.05, -0.36, r"$x_c=x_0+\xi z$", fontsize=11.5, ha="left", va="top")
    # radius and azimuth from the ring center
    ph0 = np.radians(35)
    pr = np.array([xc0 + np.cos(ph0), np.sin(ph0)])
    ax.plot([xc0, pr[0]], [0, pr[1]], color="k", lw=1.0)
    ax.text(xc0 + 0.45 * np.cos(ph0) - 0.02, 0.45 * np.sin(ph0) + 0.04, r"$R(z)$", fontsize=13.5,
            ha="right")
    tm = arc(ax, (xc0, 0), 0.3, 0, ph0)
    ax.text(xc0 + 0.4 * np.cos(tm), 0.4 * np.sin(tm), r"$\phi$", fontsize=13.5, va="center")
    rr = np.array([np.cos(ph0), np.sin(ph0)])
    tt = np.array([-np.sin(ph0), np.cos(ph0)])
    p0 = pr + 0.05 * rr                      # start just outside the ring
    arrow(ax, p0, p0 + 0.5 * rr, color=PURPLE, lw=2.2)
    ax.text(*(p0 + 0.56 * rr), r"$V_r$", color=PURPLE, fontsize=13.5)
    arrow(ax, p0, p0 + 0.55 * tt, color=ORANGE, lw=2.2)
    ax.text(*(p0 + 0.6 * tt), r"$V_\phi$  ", color=ORANGE, fontsize=13.5, ha="right")
    ax.text(-0.08, 1.4, "far side", ha="right", va="bottom",
            fontsize=10.1, color=RED)
    ax.text(xc0, -1.08, "near side (toward the observer)", ha="center",
            va="top", fontsize=10.1, color=RED)
    a_ = np.radians(np.linspace(170, 212, 30))
    if dsense < 0:
        a_ = a_[::-1]
    rD = 1.3
    ax.plot(rD * np.cos(a_), rD * np.sin(a_), color="0.4", lw=1.3)
    arrow(ax, (rD * np.cos(a_[-3]), rD * np.sin(a_[-3])), (rD * np.cos(a_[-1]), rD * np.sin(a_[-1])),
          color="0.4", lw=1.3)
    ax.text(rD * np.cos(np.radians(212)) - 0.12, rD * np.sin(np.radians(212)) - 0.1,
            "disk\nrotation", fontsize=10.1, color="0.4", ha="center", va="top")
    ax.set_xlim(-1.5, 1.88)
    ax.set_ylim(-1.4, 1.65)
    ax.set_aspect("equal")
    ax.axis("off")
    titles[2] = "(c) one ring, seen from the tip of the lobe"

    # ---------------- (d) PV curve of the cut ----------------
    ax = axs[3]
    phi = np.linspace(0, 2 * np.pi, 721)
    z = cut_z(phi)
    R = Rz(z)
    if VM is not None:
        Vzz, Vrr = VM.vel(np.where(np.isfinite(z), z, 1.0))
    else:
        f = np.clip(z / ZREF, 1e-3, None) ** P["qv"]
        Vzz, Vrr = P["Vz"] * f, P["Vr"] * f
    X = (P["xoff_au"] + P["xslope"] * z + R * np.cos(phi)) / dist
    v = vsys - s * Vzz * ci + Vrr * si * np.sin(phi) + P["j"] / R * si * np.cos(phi)
    ok = np.isfinite(z)
    far = ok & (np.sin(phi) >= 0)
    near = ok & (np.sin(phi) < 0)
    ax.plot(np.where(far, X, np.nan), np.where(far, v, np.nan), color=RED, lw=2,
            label="far half")
    ax.plot(np.where(near, X, np.nan), np.where(near, v, np.nan), color=RED, lw=1.2,
            ls=(0, (3, 2)), label="near half")
    zbar = Ycut / si
    vcen = vsys - s * (VM.vel(np.array([zbar]))[0][0] if VM is not None
                       else P["Vz"] * (zbar / ZREF) ** P["qv"]) * ci
    xc = (P["xoff_au"] + P["xslope"] * zbar) / dist
    xr = np.nanmax(np.abs(X - xc))
    ax.set_xlim(xc - 1.35 * xr, xc + 1.35 * xr)
    vlo, vhi = np.nanmin(v), np.nanmax(v)
    ax.set_ylim(min(vlo, vsys) - 0.40 * (vhi - vlo), vhi + 0.30 * (vhi - vlo))
    ax.axhline(vsys, color="0.5", lw=0.8, ls=":")
    ax.text(xc - 1.3 * xr, vsys, r"$V_{\rm sys}$", fontsize=12.2, va="bottom", color="0.4")
    ax.axhline(vcen, color="k", lw=0.7, ls="--")
    ax.axvline(xc, color="k", lw=0.7, ls="--")
    # V_cent below its line, so it stays clear of the 180-degree label above it
    ax.annotate(r"$V_{\rm cent}$", (xc - 1.3 * xr, vcen), xytext=(0, -4), textcoords="offset points",
                fontsize=12.2, va="top", ha="left")
    ax.text(0.03, 0.03, "ring-center velocity:\n" r"$V_{\rm cent}=V_{\rm sys}-s\,V_z\cos i$",
            transform=ax.transAxes, fontsize=10.8, ha="left", va="bottom")
    ax.annotate(r"$x_c$", (xc, ax.get_ylim()[0]), xytext=(4, 5), textcoords="offset points",
                fontsize=12.2, va="bottom", ha="left")
    for deg in (0, 90, 180, 270):
        k = int(round(deg / 360 * (len(phi) - 1)))
        if ok[k]:
            ax.plot(X[k], v[k], "o", color=RED, ms=5, zorder=4)
            off = {0: (22, 14), 90: (48, 22), 180: (-28, 4), 270: (24, -12)}[deg]
            lab_ = rf"$\phi$ = {deg}°" if deg == 90 else f"{deg}°"
            ax.annotate(lab_, (X[k], v[k]), xytext=off,
                        textcoords="offset points", ha="center", va="center", fontsize=11.5)
    iv = np.where(far)[0]
    if len(iv) > 10:
        # synthetic trace points along the detected half, to show what is fitted
        rng_ = np.random.default_rng(4)
        kk = iv[np.linspace(3, len(iv) - 4, 16).astype(int)]
        sv = 0.035 * (vhi - vlo)
        ax.errorbar(X[kk] + 0.01 * xr * rng_.standard_normal(len(kk)),
                    v[kk] + sv * rng_.standard_normal(len(kk)), yerr=sv, fmt="o", color="k",
                    ms=3, lw=0.7, capsize=1.2, zorder=3, label=r"trace points $(x_k, v_k)$")
        k0 = iv[int(0.28 * len(iv))]
        px, pv = X[k0] - 0.1 * xr, v[k0] + 0.12 * (vhi - vlo)
        d = np.hypot((X[far] - px) / (0.05 * xr), (v[far] - pv) / (0.05 * (vhi - vlo)))
        m = np.argmin(d)
        ax.errorbar([px], [pv], yerr=0.04 * (vhi - vlo), fmt="o", color=BLUE, ms=4, capsize=1.5,
                    lw=0.8)
        ax.plot([px, X[far][m]], [pv, v[far][m]], color=BLUE, lw=1.3)
        ax.text(px, pv, r"$d_k$  ", color=BLUE, fontsize=13.5, ha="right", va="bottom")
    ax.set_xlabel(r"position along the cut, $X$ (arcsec)")
    ax.set_ylabel(r"line-of-sight velocity, $v$ (km s$^{-1}$)")
    ax.legend(fontsize=9.5, loc="lower right", frameon=False, borderaxespad=0.2,
              handlelength=1.6, labelspacing=0.25)
    titles[3] = rf"(d) PV diagram of the cut at $y_0$ = {y0:g}$''$"
    fig.canvas.draw()
    for row in ((0, 1), (2, 3)):
        ytop = max(axs[k].get_position().y1 for k in row) + 0.012
        for k in row:
            fig.text(axs[k].get_position().x0, ytop, titles[k], fontsize=13.5, ha="left",
                     va="bottom")

    out = wk.figpath(cfg, "model_definitions")
    fig.savefig(out, dpi=cfg.FIG_DPI)
    print(f"[defs] cone from {src} values; wrote {out}")


if __name__ == "__main__":
    main()
