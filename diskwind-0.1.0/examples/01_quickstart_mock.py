"""
01 -- Quick start on a synthetic cube with known parameters (about 2 minutes).

    python 01_quickstart_mock.py

Makes a mock disk + conical wind, fits it, and compares with the truth.
"""
import numpy as np

from diskwind import WindAnalysis
from diskwind.mock import MockWind, make_mock_cube

# 1. a synthetic cube (full rings, i = 49 deg) and a matching configuration
truth = MockWind(theta_deg=30.0, vp=8.0, j=-40.0)
make_mock_cube("mock_cube.fits", truth)
cfg = truth.config("mock_cube.fits")
cfg.output.dir = "out_01_mock"
cfg.fit.n_boot = 0                       # covariance errors only: fast
cfg.save("mock.toml")                    # look at it: every option is listed

# 2. the analysis
wa = WindAnalysis(cfg)
print(wa.cube)                           # Vsys, lobe orientation, disk sense estimated
res = wa.fit("constant", radial="cone")  # one cone, constant speed, flow along the cone
print(res.summary())

# 3. compare with the truth
t = truth.truth()
print(f"theta : {res.theta_deg[0]:.1f}  (true {t['theta_deg']:.1f}) deg")
print(f"Vz    : {res.params['Vz']:.2f}  (true {t['Vz']:.2f}) km/s")
print(f"j     : {res.params['j']:.1f}  (true {t['j']:.1f}) au km/s")

# 4. figures (written to out_01_mock/)
wa.plot_cuts(res)
wa.plot_profiles(res)
