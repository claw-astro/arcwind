"""Rotation measured cut by cut, and fixed parameters outside the default bounds."""
import matplotlib
matplotlib.use("Agg")
import numpy as np


def _const(analysis):
    if "t-const" not in analysis.results:
        analysis.fit("constant", radial="cone", n_boot=0, label="t-const")
    return analysis.results["t-const"]


def test_rotation_per_cut_recovers_constant_j(analysis):
    r = _const(analysis)
    from arcwind import rotation_per_cut
    rot = rotation_per_cut(r)
    assert len(rot.y0) >= 3
    good = rot.detected
    assert good.any()
    # each cut alone agrees with the global constant j within 3 sigma
    jg = r.params["j"]
    assert np.all(np.abs(rot.j[good] - jg) < 3 * rot.j_err[good] + 1e-6)
    assert "detected" in rot.joint
    txt = rot.summary()
    assert "rotation per cut" in txt
    import arcwind.plotting as P
    P.plot_rotation(rot, r)


def test_fixed_parameter_outside_default_bounds(analysis):
    """A fixed value is kept even if the default bounds of these cuts exclude it."""
    from arcwind import ConicalWindModel, GlobalFitter
    r = _const(analysis)
    y0 = sorted(r.traces)[0]
    m = ConicalWindModel(r.model.law, radial=r.model.radial, rotation="constant")
    fix = {k: v for k, v in r.params.items() if k != "j"}
    m.fix(**fix)
    f = GlobalFitter(m, r.geometry, {y0: r.traces[y0]}, r.sig_x, verbose=False, n_jobs=1)
    for i, k in enumerate(f.names):
        if k in fix:
            assert f.lo[i] <= fix[k] <= f.hi[i]
    res = f.fit(init={"j": r.params["j"]})
    for k, v in fix.items():
        assert abs(res.params[k] - v) < 1e-9


def test_contiguous_flag():
    from arcwind.rotation import RotationProfile
    n = 5
    a = np.arange(n, dtype=float)
    j = np.array([0, 50, 0, 50, 0.0])
    p = RotationProfile("x", a, a * 50 + 50, a * 10 + 30, j, np.full(n, 10.0),
                        np.full(n, 20), np.ones(n), np.zeros(n), 0.75, 0.75)
    assert not p.contiguous
    p.j = np.array([0, 50, 50, 0, 0.0])
    assert p.contiguous


def test_median():
    from arcwind.rotation import RotationProfile
    n = 5
    a = np.arange(n, dtype=float)
    p = RotationProfile("x", a, a * 50 + 50, a * 10 + 30, np.array([1, 2, 3, 4, 100.0]),
                        np.full(n, 10.0), np.full(n, 20), np.array([1, 1, 1, 1, 5.0]),
                        np.zeros(n), 0.75, 0.75)
    m = p.median()
    assert m["n"] == 4 and m["value"] == 2.5          # the poor cut (d2 = 5) is left out
    assert m["lo"] <= m["value"] <= m["hi"]
