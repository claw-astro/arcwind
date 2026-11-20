from arcwind.cli import main
from arcwind import WindConfig


def test_init_and_options(tmp_path, capsys):
    p = tmp_path / "x.toml"
    main(["init", str(p)])
    c = WindConfig.load(str(p))
    assert c.source.name == "My source"
    main(["options", "trace"])
    assert "snr_min" in capsys.readouterr().out
