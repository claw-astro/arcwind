
import matplotlib
import pytest

matplotlib.use("Agg")


@pytest.fixture(scope="session")
def mock_files(tmp_path_factory):
    """A small synthetic cube with known parameters and a matching configuration."""
    from arcwind.mock import MockWind, make_mock_cube
    d = tmp_path_factory.mktemp("mock")
    m = MockWind()
    fits = make_mock_cube(str(d / "mock.fits"), m, n_wind=400_000, n_disk=400_000)
    cfg = m.config(fits)
    cfg.output.dir = str(d / "out")
    cfg.output.verbose = False
    cfg.cuts.y_offsets_arcsec = [0.45, 0.6, 0.75, 1.05]
    cfg.fit.n_boot = 0
    return m, cfg


@pytest.fixture(scope="session")
def analysis(mock_files):
    from arcwind import WindAnalysis
    return WindAnalysis(mock_files[1])
