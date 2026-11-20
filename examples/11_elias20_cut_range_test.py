"""
11 -- Robustness of the Elias 2-20 fit to the highest, faintest cuts.

    python 11_elias20_cut_range_test.py      # needs Elias20_CO.fits next to configs/

Fits the same configuration (inner-arch window, cuts from y0 = 0.45") with the cuts
restricted to y0 <= YMAX, for each YMAX in YMAXES (default: 2.0" and 2.4", the
highest cut), and compares the results.

Every fit starts from scratch, as in example 02: an unguided single power law
fitted to the same cuts, then the broken power law guided by it.

Every file carries its label:

  broken-power-law_ymax<YMAX>.json/.txt                one per YMAX
  power-law-first_ymax<YMAX>.json/.txt                 its starting fit
  broken-power-law_ymax<YMAX>_cuts, _profiles, _angular_momentum, _bootstrap,
                    _rotation_per_cut (.txt and figure), _cut_y<y0>
  compare_ymax.txt                                     parameters side by side, derived
                                                       quantities and the sign of j
  compare_ymax_profiles, compare_ymax_angular_momentum  profiles and j(z) of all fits
"""
import matplotlib.pyplot as plt
import numpy as np

from diskwind import WindAnalysis
from diskwind.plotting import plot_angular_momentum, plot_profiles

YMAXES = (2.0, 2.4)                          # highest cut of each fit (arcsec)
CONFIG = "configs/elias20.toml"
GM_SUN = 886.9                               # G M_sun in au km^2 s^-2

wa = WindAnalysis(CONFIG)
ys_all = [float(y) for y in wa.config.y_offsets]

# --- one fit per YMAX, each from scratch ---------------------------------------
fits, rots = {}, {}
for ymax in YMAXES:
    tag = f"ymax{ymax:.2f}"
    ys = [y for y in ys_all if y <= ymax + 1e-6]
    print(f"\n=== {tag}: cuts {ys[0]:.2f}-{ys[-1]:.2f} arcsec ===")
    try:
        start = wa.fit("power-law", n_boot=0, y_offsets=ys, label=f"power-law-first_{tag}")
        r = wa.fit("broken-power-law", guide=start, start_from=start, iterations=2,
                   y_offsets=ys, label=f"broken-power-law_{tag}")
    except Exception as e:                   # e.g. too few usable cuts
        print(f"{tag}: fit failed ({type(e).__name__}: {e})")
        continue
    print(r.summary())
    wa.plot_cuts(r)
    for y0 in r.traces:
        wa.plot_cut(y0, r)
        plt.close("all")
    wa.plot_profiles(r)
    wa.plot_angular_momentum(r)
    try:
        rots[tag] = wa.rotation_per_cut(r)
        wa.plot_rotation(rots[tag], result=r)
    except ValueError as e:
        print(f"{tag}: rotation per cut failed ({e})")
    if r.boot is not None and len(r.boot):
        wa.plot_bootstrap(r)
    fits[tag] = r


# --- comparison ------------------------------------------------------------
def derived(r):
    """A few derived quantities of one fit."""
    g, p, e = r.geometry, r.params, r.errors("max")
    th, eth = r.theta_deg
    z = np.array([g.zc_au(min(r.traces)), g.zc_au(max(r.traces))])
    P = r.profiles(z)
    gm = GM_SUN * g.mstar_msun
    jz = r.model.j_of_z(z)                   # j at the lowest and highest cut (j(z) if qj free)
    Vphi = P["Vphi"]
    out = {"cuts used": f"{len(r.traces)} ({min(r.traces):.2f}-{max(r.traces):.2f}\")",
           "trace points": f"{r.n}",
           "chi2_red": f"{r.chi2r:.2f}",
           "theta (deg)": f"{th:.2f} +/- {eth:.2f}",
           "Vp at lowest / highest cut (km/s)": f"{P['Vp'][0]:.2f} / {P['Vp'][1]:.2f}",
           "j at lowest / highest cut": f"{jz[0]:.1f} / {jz[1]:.1f}",
           "Vphi at lowest / highest cut (km/s)": f"{Vphi[0]:.2f} / {Vphi[1]:.2f}",
           "r_L = j^2/GM* at lowest cut (au)": f"{jz[0] ** 2 / gm:.1f}",
           "worst cut (mean d^2)": "{:.2f}\" ({:.1f})".format(
               *max(r.per_cut_d2.items(), key=lambda kv: kv[1]))}
    rot = rots.get(r.label.split("_")[-1])
    if rot is not None:
        med = rot.median()
        out["median j, well-fitted cuts"] = (f"{med['value']:.1f} "
                                             f"[{med['lo']:.0f}, {med['hi']:.0f}]")
        det = rot.detected
        out["rotation detected (cuts)"] = (", ".join(f"{y:.2f}" for y in rot.y0[det])
                                           or "none")
        d, u = rot.joint.get("detected"), rot.joint.get("undetected")
        if d:
            out["joint j, detected cuts"] = f"{d['j']:.1f} +/- {d['err']:.1f}"
        if u:
            out["|j| limit elsewhere (2 sig)"] = f"< {abs(u['j']) + 2 * u['err']:.0f}"
    if r.boot is not None and len(r.boot):
        jb = r.boot[:, r.names.index("j")]
        out["j < 0 in samples (%)"] = f"{100 * np.mean(jb < 0):.0f}"   # same sign at all z
        lo, hi = np.percentile(jb, [2.5, 97.5])
        out["j(z_ref) 95% range (au km/s)"] = f"{lo:.0f} to {hi:.0f}"
    return out


if fits:
    tags = list(fits)
    w = 26
    lines = ["# Elias 2-20: the adopted analysis with the cuts restricted to y0 <= YMAX",
             "# errors: larger of covariance and bootstrap",
             f"{'parameter':30s}" + "".join(f"{t:>{w}s}" for t in tags)]
    names = fits[tags[0]].names
    for k in names:
        row = f"{k:30s}"
        for t in tags:
            r = fits[t]
            row += f"{r.params[k]:>{w - 13}.4g} +/- {r.errors('max')[k]:<8.3g}"
        lines.append(row)
    lines.append("")
    ds = {t: derived(fits[t]) for t in tags}
    for k in ds[tags[0]]:
        lines.append(f"{k:30s}" + "".join(f"{ds[t].get(k, '-'):>{w}s}" for t in tags))
    lines += ["",
              "NB: chi2 values are not comparable between fits (different points).",
              "    q2 and z_b rest on the highest cuts; with fewer cuts they are poorly "
              "constrained, and z_b sits at its upper bound (1.1 x the highest cut height)."]
    text = "\n".join(lines)
    print("\n" + text)
    with open(wa.config.outpath("compare_ymax", "txt"), "w") as f:
        f.write(text + "\n")
    rs = [fits[t] for t in tags]
    plot_profiles(rs[-1], others=rs[:-1], save=wa.config.outpath("compare_ymax_profiles"))
    plot_angular_momentum(rs[-1], others=rs[:-1],
                          save=wa.config.outpath("compare_ymax_angular_momentum"))
    print(f"\nwrote {wa.config.outpath('compare_ymax', 'txt')}, compare_ymax_profiles and "
          "compare_ymax_angular_momentum")
