"""
06 -- The building blocks, without WindAnalysis.

Useful to embed the method in your own pipeline or to try new ideas.
"""
import numpy as np

from arcwind import (ConicalWindModel, GlobalFitter, TraceExtractor, WindConfig,
                      WindCube)
from arcwind.constants import FWHM2SIG

cfg = WindConfig.load("mock.toml")
cube = WindCube(cfg)                                  # resampled cube, noise, Vsys ...
geom = cube.geometry                                  # inclination, s, distance, Vsys

pv = cube.pv(1.0)                                     # PVDiagram: data, rms, mask
print(pv.data.shape, pv.rms, pv.excluded.mean())

ext = TraceExtractor(cube)
traces = ext.extract_all(np.arange(0.3, 1.6, 0.15))   # {y0: Trace}

model = ConicalWindModel("power-law", radial="cone")
model.fix(xslope=0.0).set_bounds(Vz=(0, 20))
fitter = GlobalFitter(model, geom, traces, sig_x=0.5 * cube.beam_x * FWHM2SIG,
                      loss="soft_l1", equal_cut_weight=True)
res = fitter.fit()
fitter.bootstrap(res, n_boot=10)
print(res.summary())

curve = res.curve(1.0)                                # the model in the cut y0 = 1"
print(curve.X[curve.valid][:5], curve.v[curve.valid][:5])
print(res.profiles(np.array([50.0, 200.0]))["Vp"])
res.save("out_06_result.json")
