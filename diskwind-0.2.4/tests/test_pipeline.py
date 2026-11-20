"""End-to-end tests on a small synthetic cube with known parameters."""
import numpy as np

from diskwind import FitResult, compare_fits
from diskwind.mass import column_density_co, column_density_up


def test_cube_estimates(analysis):
    cube = analysis.cube
    assert abs(cube.vsys - 3.0) < 0.2
    assert cube.s == +1 and cube.dsense == -1
    pv = cube.pv(1.0)
    assert pv.data.shape == (len(cube.v), len(cube.x)) and pv.rms > 0
    assert pv.excluded.shape == pv.data.shape


def test_noise_and_rotation_checks(analysis):
    rep = analysis.check_noise(plot=False)
    assert rep.ok
    chk = analysis.check_disk_rotation(plot=False)
    assert chk.best_sense == -1


def test_traces(analysis):
    tr = analysis.extract_traces()
    assert len(tr) == 4 and all(t.n > 5 for t in tr.values())
    t = next(iter(tr.values()))
    assert set(np.unique(t.side)) <= {-1, 1}


def test_global_fit_recovers_truth(analysis, mock_files, tmp_path):
    truth = mock_files[0].truth()
    res = analysis.fit("constant", radial="cone", n_boot=0, label="t-const")
    p = res.params
    assert abs(np.degrees(np.arctan(p["tan_theta"])) - truth["theta_deg"]) < 2.0
    assert abs(p["Vz"] - truth["Vz"]) / truth["Vz"] < 0.1
    assert np.sign(p["j"]) == np.sign(truth["j"]) and abs(p["j"] - truth["j"]) < 20
    path = res.save(str(tmp_path / "r.json"))
    r2 = FitResult.load(path)
    assert np.allclose(r2.model.vector(), res.model.vector())
    assert len(r2.traces) == len(res.traces)
    assert "theta" in res.summary()


def test_bootstrap_and_comparison(analysis):
    a = analysis.results["t-const"]
    b = analysis.fit("power-law", guide=a, start_from=a, n_boot=3, fixed={"xslope": 0.0},
                     label="t-pl")
    assert b.boot is not None and b.boot.shape == (3, len(b.names))
    rows = compare_fits([a, b])
    assert {r["label"] for r in rows} == {"t-const", "t-pl"}
    prof = b.profiles(np.linspace(50, 300, 10))
    assert prof["Vp_lo"].shape == (10,)


def test_resampling_parallel_serial_and_jackknife(analysis):
    a = analysis.results["t-const"]
    ser = analysis.fit("constant", start_from=a, n_boot=4, n_jobs=1, label="t-ser")
    par = analysis.fit("constant", start_from=a, n_boot=4, n_jobs=2, label="t-par")
    # same resampled cuts (same seed): same samples in parallel and serially
    assert np.allclose(ser.boot, par.boot, rtol=1e-6, atol=1e-8)
    jk = analysis.fit("constant", start_from=a, n_boot=1, error_method="jackknife",
                      label="t-jk")
    assert jk.boot.shape[0] == len(jk.traces)
    assert jk.settings["error_method"] == "jackknife"
    assert np.allclose(jk.boot.mean(axis=0), jk.model.vector(), atol=0.5)


def test_mass_and_derived(analysis):
    r = analysis.results["t-const"]
    m = analysis.mass(r, plot=False)
    assert m.total_mass > 0 and np.isfinite(m.median_mdot)
    props = analysis.derived(r, m)
    assert props.launch_radius()[0] > 0
    assert "HH 30" in props.compare_with("HH 30")


def test_co_column_densities_match_louvet2018():
    assert np.isclose(column_density_up(38.7), 5.8e15, rtol=0.02)
    assert np.isclose(column_density_co(38.7, 30.0) / 1e-4, 2.25e20, rtol=0.01)


def test_ring_fit_full_rings(analysis):
    from diskwind import RingFitter
    f = RingFitter(analysis.cube).fit_cut(0.75, nmc=50)
    assert f.ok and abs(f.Vr - 4.0) < 0.6


def test_figures(analysis, tmp_path):
    import diskwind.plotting as P
    r = analysis.results["t-const"]
    for fig in (P.plot_cuts(analysis.cube, r), P.plot_profiles(r),
                P.plot_gallery(analysis.cube, result=r), P.plot_slit_layout(analysis.cube)):
        fig.savefig(tmp_path / "f.png")
