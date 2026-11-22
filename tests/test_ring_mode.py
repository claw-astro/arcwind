"""trace.mode = "ring": walls traced channel by channel."""
import numpy as np
import pytest

from arcwind.trace import Trace


def test_trace_xerr_defaults_and_roundtrip(tmp_path):
    t = Trace(1.0, np.array([0.1, 0.2]), np.array([1.0, 2.0]), np.array([0.1, 0.1]),
              np.array([-1, 2]), np.array([5.0, 6.0]), np.array([0.5, 0.2]))
    assert np.isnan(t.xerr).all()                          # old 6-field construction
    t.xerr = np.array([np.nan, 0.03])
    assert np.allclose(t.x_errors(0.05), [0.05, 0.03])
    Trace.save_all({1.0: t}, tmp_path / "t.npz")
    back = Trace.load_all(tmp_path / "t.npz")[1.0]
    assert np.allclose(back.x_errors(0.05), [0.05, 0.03])
    # thinning keeps wall points (sampled in velocity)
    assert set(t.thinned(1.0, 0.0, 0.1).side) == {2}


@pytest.fixture(scope="module")
def tall(tmp_path_factory):
    from arcwind.mock import MockWind, make_mock_cube
    d = tmp_path_factory.mktemp("tall")
    m = MockWind.tall_ring_like()
    fits = make_mock_cube(str(d / "tall.fits"), m, n_wind=500_000, n_disk=300_000)
    cfg = m.config(fits)
    cfg.output.dir = str(d / "out")
    cfg.output.verbose = False
    cfg.fit.n_boot = 0
    cfg.cuts.y_offsets_arcsec = [0.6, 0.9, 1.2, 1.5]
    return m, cfg


def _analysis(cfg, mode):
    from copy import deepcopy
    from arcwind import WindAnalysis
    c = deepcopy(cfg)
    c.trace.mode = mode
    return WindAnalysis(c)


def test_ring_mode_traces_the_walls(tall):
    from arcwind import TraceExtractor
    m, cfg = tall
    n = {}
    for mode in ("wings", "ring"):
        tr = TraceExtractor(_analysis(cfg, mode).cube).extract_all(min_points=0)
        n[mode] = sum(t.n for t in tr.values())
        if mode == "ring":
            sides = np.concatenate([t.side for t in tr.values()])
            assert {-2, 2} <= set(sides) and {-1, 1} & set(sides)
    assert n["ring"] > 3 * n["wings"]


def test_ring_mode_recovers_tall_ring_wind(tall):
    m, cfg = tall
    t = m.truth()
    r = _analysis(cfg, "ring").fit("constant")
    p = r.params
    assert abs(p["tan_theta"] - t["tan_theta"]) < 0.03
    assert abs(p["Vz"] - t["Vz"]) < 1.0
    assert np.sign(p["j"]) == np.sign(t["j"]) and abs(p["j"] - t["j"]) < 20
    assert (np.abs(np.concatenate([tr.side for tr in r.traces.values()])) == 2).any()
    # the result (with the wall points' position errors) survives a save/load
    from arcwind.fitting import FitResult
    path = r.save(str(tmp_path_of(cfg) / "ring.json"))
    back = FitResult.load(path)
    y = next(iter(r.traces))
    assert np.allclose(back.traces[y].x_errors(1.0), r.traces[y].x_errors(1.0))


def tmp_path_of(cfg):
    from pathlib import Path
    p = Path(cfg.output.dir)
    p.mkdir(parents=True, exist_ok=True)
    return p
