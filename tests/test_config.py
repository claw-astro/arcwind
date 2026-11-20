import pytest

from arcwind import WindConfig
from arcwind.config import ConfigError


def test_defaults_and_set_get():
    c = WindConfig()
    assert c.trace.mode == "wings" and c.fit.velocity_law == "power-law"
    c.set(snr_min=3, **{"fit.n_boot": 0})
    assert c.get("trace.snr_min") == 3 and c.fit.n_boot == 0


def test_unknown_option_suggests():
    c = WindConfig()
    with pytest.raises(ConfigError, match="did you mean"):
        c.set(snr_mni=3)
    with pytest.raises(ConfigError, match="unknown section"):
        WindConfig.from_dict({"trac": {"snr_min": 3}})


def test_per_cut_override():
    c = WindConfig()
    c.mask.x_cone = (0.1, 0.3, 0.5)
    c.cuts.overrides = {2.0: {"snr_min": 3.5, "x_range_arcsec": (-1.8, 1.8)}}
    c._normalize()
    cc = c.for_cut(2.0)
    assert cc.trace.snr_min == 3.5 and cc.mask.x_cone is None
    assert c.for_cut(1.0).trace.snr_min == 5.0 and c.trace.snr_min == 5.0


def test_toml_roundtrip(tmp_path):
    c = WindConfig.minimal("cube.fits", 140, 0.5, 45, 150, 240, snr_min=4,
                           **{"cuts.y_offsets_arcsec": {"start": 0.3, "stop": 0.9, "step": 0.3}})
    c.cuts.overrides = {"0.60": {"trace.snr_min": 3.0}}
    p = c.save(str(tmp_path / "c.toml"))
    c2 = WindConfig.load(p)
    assert c2.to_dict() == c.to_dict()
    assert list(c2.y_offsets) == [0.3, 0.6, 0.9]
    assert c2.for_cut(0.6).trace.snr_min == 3.0


def test_legacy_conversion(tmp_path):
    p = tmp_path / "old_config.py"
    p.write_text("FITS_FILE='a.fits'\nDISTANCE_PC=140\nMSTAR_MSUN=0.5\nINCLINATION_DEG=45\n"
                 "DISK_PA_DEG=150\nWIND_PA_DEG=240\nTRACE_MODE='single'\nVELOCITY_MODEL='bpl-cone'\n"
                 "GLOBAL_FIXED={'qj': 0.0}\nTRACE_SLIT_OVERRIDES={2.0: dict(TRACE_SNR_MIN=3.5)}\n")
    c = WindConfig.load(str(p))
    assert c.trace.mode == "single" and c.fit.velocity_law == "broken-power-law"
    assert c.fit.rotation == "constant" and c.for_cut(2.0).trace.snr_min == 3.5


def test_describe_lists_all_sections():
    txt = WindConfig.describe()
    for sec in ("[source]", "[trace]", "[fit]", "[mass]"):
        assert sec in txt
    assert "list-table" in WindConfig.describe("fit", fmt="rst")


def test_validate_requires_source():
    with pytest.raises(ConfigError, match="missing required"):
        WindConfig().validate()
