"""
05 -- Independent ring fits to each cut (Louvet et al. 2018, Sect. 4.4).

For winds where both halves of every ring are detected (e.g. edge-on HH 30),
each cut can be fitted on its own; the cone is then fitted to R(z).

Data: the synthetic cube of example 01 (mock_cube.fits + mock.toml; run
01_quickstart_mock.py first): i = 49 deg, full rings, theta = 30 deg,
Vp = 8 km/s, j = -40 au km/s.
"""
from arcwind import WindAnalysis
from arcwind.plotting import plot_ring_fits

wa = WindAnalysis("mock.toml", **{"output.dir": "out_05_rings"})
rings = wa.ring_fits()                     # two passes: free theta, then theta fixed
print(rings.cone)                          # R0, theta from R(z)
for f in rings.good():
    print(f"y0 = {f.y0:.2f}: z = {f.z_au:.0f} au, Vz = {f.Vz:.2f} +/- {f.errors['Vz']:.2f}, "
          f"Vphi = {f.Vphi:.2f}, j = {f.j:.0f}")
plot_ring_fits(rings, save="out_05_rings/ring_fits.pdf")

# constrain the per-cut fits by a global fit (priors on Vcent and Vr sin i)
g = wa.fit(n_boot=0)
rings2 = wa.ring_fits(prior_from=g)
