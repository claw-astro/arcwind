"""
07 -- An HH 30-like bipolar wind: nearly edge-on, both lobes fitted.

    python 07_hh30_like_both_lobes.py        (~5 minutes; no data needed)

Makes a synthetic disk and bipolar conical wind with the parameters of HH 30
(Louvet et al. 2018), then

  1. fits one mirror-symmetric wind to the cuts of BOTH lobes;
  2. fits each lobe on its own and compares them;
  3. runs the independent per-cut ring fits of Louvet et al. (full rings);
  4. computes the mass of each lobe and makes the figures of both lobes.
"""
import numpy as np

from diskwind import WindAnalysis
from diskwind.mock import MockWind, make_mock_cube

# --- 1. the synthetic cube --------------------------------------------------
truth = MockWind.hh30_like()             # i = 89 deg (91 deg), theta = 35 deg,
make_mock_cube("hh30_mock.fits", truth)  # Vp = 9.3 km/s, j = 38 au km/s, both lobes
cfg = truth.config("hh30_mock.fits")     # sets source.lobes = "both"
cfg.output.dir = "out_07_hh30"
cfg.fit.n_boot = 0
cfg.save("hh30_mock.toml")
t = truth.truth()

wa = WindAnalysis(cfg)
print(wa.cube)                           # main lobe (NE, PA 31.2 deg)
print(wa.counter_cube)                   # counter lobe: frame turned by 180 deg
wa.plot_slit_layout()                    # cuts on both sides of the disk

# --- 2. one wind, both lobes ------------------------------------------------
both = wa.fit("constant", lobe="both")   # Louvet+18: constant poloidal speed
print(both.summary())
Vp = both.params["Vz"] / np.cos(np.arctan(both.params["tan_theta"]))
print(f"theta {both.theta_deg[0]:.1f} (true {t['theta_deg']:.0f}) deg, "
      f"Vp {Vp:.2f} (true {truth.vp}) km/s, j {both.params['j']:.1f} (true {t['j']:.0f})")

# does the speed change with height? (the same comparison as for Elias 2-20)
pl = wa.fit("power-law", lobe="both", start_from=both)
wa.compare(both, pl)

# --- 3. each lobe on its own --------------------------------------------------
north = wa.fit("constant", lobe="main")
south = wa.fit("constant", lobe="counter")
south_in_main_frame = south.view("main")     # same sign conventions as `north`
for name, r in (("main", north), ("counter", south_in_main_frame)):
    print(f"{name:8s}: theta {r.theta_deg[0]:.1f} deg, Vz {r.params['Vz']:.2f} km/s, "
          f"j {r.params['j']:+.1f} au km/s")

# --- 4. per-cut ring fits, as in Louvet et al. (2018) -----------------------
rings = wa.ring_fits(lobe="main")
print("cone from R(z):", {k: round(float(v), 2) for k, v in rings.cone.items()})

# --- 5. mass and figures of both lobes --------------------------------------
for lobe in ("main", "counter"):
    wa.mass(both, lobe=lobe)                      # prints the total mass of the lobe
    wa.plot_cuts(both, lobe=lobe)                 # out_07_hh30/constant-both_cuts*.pdf
    wa.plot_gallery(both, lobe=lobe)
wa.plot_profiles(both, others=[pl])
