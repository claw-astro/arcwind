"""
Helpers for the schematic examples (08-10): read the configuration and the
saved fit in the form these drawing scripts use.

* ``load_config(path)`` reads a diskwind configuration and exposes the old
  upper-case names (``cfg.INCLINATION_DEG`` ...).
* ``VelocityModel(cfg, model)`` wraps the saved :class:`diskwind.FitResult`
  ``<output>/<label>.json`` (label: the velocity law, e.g. broken-power-law).
"""
import os
from types import SimpleNamespace

import numpy as np

from diskwind import FitResult, WindConfig
from diskwind.config import LEGACY_KEYS

GLOBAL_ZREF_AU = 100.0
VELOCITY_MODELS = ("single", "bpl-cone")
_LABEL = {"single": "power-law", "bpl-cone": "broken-power-law"}


def load_config(path):
    cfg = WindConfig.load(path)
    d = {}
    for old, new in LEGACY_KEYS.items():
        if new:
            d[old] = cfg.get(new)
    d["OUTPUT_DIR"] = cfg.output_dir
    d["WIND_APPROACHING"] = cfg.source.wind_approaching
    d["FIT_Y_OFFSETS_ARCSEC"] = cfg.y_offsets
    d["VELOCITY_MODEL"] = "bpl-cone" if cfg.fit.velocity_law == "broken-power-law" else "single"
    d["_cfg"] = cfg
    return SimpleNamespace(**d)


def figpath(cfg, name):
    return cfg._cfg.outpath(name)


def setup_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    from diskwind.plotting._style import plt
    return plt()


class VelocityModel:
    """The subset of the old VelocityModel API used by the schematic scripts."""

    def __init__(self, cfg, model="bpl-cone", n_draw=300):
        label = _LABEL.get(model, model)
        path = os.path.join(cfg.OUTPUT_DIR, f"{label}.json")
        if not os.path.exists(path):
            raise FileNotFoundError(f"{path} not found: run the fit first")
        self.res = FitResult.load(path)
        self.model = model
        self.name = self.res.label
        self.names = self.res.names
        self.p = self.res.model.vector()
        self.draws = self.res.draws(n_draw)

    def par(self, key, p=None):
        p = self.p if p is None else np.asarray(p)
        return p[..., self.names.index(key)]

    def vel(self, z, p=None):
        return self.res.model.velocities(z, None if p is None else p)

    def jz(self, z, p=None):
        return self.res.model.j_of_z(z, None if p is None else p)


def result_path(cfg, model=None):
    """Path of the saved fit used by the schematics."""
    label = _LABEL.get(model or cfg.VELOCITY_MODEL, model or cfg.VELOCITY_MODEL)
    return os.path.join(cfg.OUTPUT_DIR, f"{label}.json")


def geometry_meta(cfg, model=None):
    """Geometry stored with the saved fit."""
    r = FitResult.load(result_path(cfg, model))
    g = r.geometry
    return dict(inc_deg=g.inclination_deg, dist_pc=g.distance_pc, s=g.s,
                disk_sense=g.disk_sense, vsys=g.vsys, y_offsets=sorted(r.traces))
