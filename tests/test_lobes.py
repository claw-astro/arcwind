"""Two-lobe (bipolar) geometry and fits."""
import numpy as np
import pytest

from diskwind import ConicalWindModel, ViewGeometry


def test_mirror_is_an_involution():
    m = ConicalWindModel("power-law", R0_au=10, tan_theta=0.5, xoff_au=5, xslope=0.1,
                         Vz=7, qv=0.2, j=40).fix(xoff_au=5).set_bounds(j=(-100, 300))
    mm = m.mirrored()
    assert mm.values["xoff_au"] == -5 and mm.values["j"] == -40
    assert mm.values["xslope"] == 0.1 and mm.fixed["xoff_au"] == -5
    assert mm.bounds["j"] == (-300, 100)
    assert mm.mirrored().values == m.values
    g = ViewGeometry(80.0, -1, 140.0, 6.9, disk_sense=1)
    assert g.mirrored().s == 1 and g.mirrored().disk_sense == -1
    assert g.mirrored().mirrored() == g


def test_counter_lobe_sees_the_same_physical_wind():
    # one point of the wind, seen from both frames: same line-of-sight velocity
    g = ViewGeometry(70.0, 1, 140.0, 0.0)
    m = ConicalWindModel("constant", radial="cone", R0_au=10, tan_theta=0.4, Vz=8, j=60)
    c = m.curve(1.0, g)
    cm = m.mirrored().curve(1.0, g.mirrored())
    # the counter lobe recedes where the main lobe approaches: ring-center
    # velocities are opposite, and the curves are mirror images in X
    ok, okm = c.valid, cm.valid
    assert np.isclose(np.nanmean(c.v[ok]), -np.nanmean(cm.v[okm]), atol=0.3)


@pytest.fixture(scope="module")
def hh30(tmp_path_factory):
    from diskwind import WindAnalysis
    from diskwind.mock import MockWind, make_mock_cube
    d = tmp_path_factory.mktemp("hh30")
    m = MockWind.hh30_like()
    fits = make_mock_cube(str(d / "hh30.fits"), m, n_wind=500_000, n_disk=300_000)
    cfg = m.config(fits)
    cfg.output.dir = str(d / "out")
    cfg.output.verbose = False
    cfg.fit.n_boot = 0
    cfg.cuts.y_offsets_arcsec = [0.6, 1.2, 1.8, 2.4]
    return m, WindAnalysis(cfg)


def test_joint_fit_recovers_hh30_like_wind(hh30):
    m, wa = hh30
    t = m.truth()
    r = wa.fit("constant", lobe="both")
    assert r.lobes == "both" and len(r.other_traces) >= 3
    assert abs(r.theta_deg[0] - t["theta_deg"]) < 2.0
    assert abs(r.params["Vz"] - t["Vz"]) / t["Vz"] < 0.1
    assert abs(r.params["j"] - t["j"]) < 10
    c = r.view("counter")
    assert c.frame == "counter" and np.isclose(c.params["j"], -r.params["j"])
    assert c.traces is r.other_traces


def test_single_counter_lobe_fit_matches(hh30):
    m, wa = hh30
    south = wa.fit("constant", lobe="counter")
    assert south.frame == "counter" and south.params["j"] < 0      # own frame
    assert south.view("main").params["j"] > 0
    assert abs(south.theta_deg[0] - m.theta_deg) < 3.0


def test_two_lobe_result_roundtrip(hh30, tmp_path):
    from diskwind import FitResult
    m, wa = hh30
    r = wa.results["constant-both"]
    r2 = FitResult.load(r.save(str(tmp_path / "both.json")))
    assert r2.lobes == "both" and len(r2.other_traces) == len(r.other_traces)
    assert np.isclose(r2.chi2_on(r2.traces, other_traces=r2.other_traces),
                      r.chi2_on(r.traces, other_traces=r.other_traces))
