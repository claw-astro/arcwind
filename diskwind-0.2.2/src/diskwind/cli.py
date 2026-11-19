"""
Command-line interface: ``diskwind <command> ...``.

Commands
--------
``init``         write a template configuration (``diskwind init mysource.toml``)
``options``      list every option with its default (``diskwind options [section]``)
``convert``      convert an old ``*_config.py`` to TOML
``center``       fit the disk center to a continuum image
``check``        noise and disk-rotation checks
``fit``          one global fit (``--law``, ``--radial``, ``--guide``, ``--nboot`` ...)
``run``          the full standard analysis
``mock``         write a synthetic cube and a matching configuration
                 (``--preset elias20`` or ``--preset hh30``)
"""
from __future__ import annotations

import argparse
import sys

__all__ = ["main"]


def _cmd_init(a):
    import os
    with open(os.path.join(os.path.dirname(__file__), "data", "template.toml")) as fh:
        txt = fh.read()
    with open(a.path, "w") as fh:
        fh.write(txt)
    print(f"wrote {a.path}: fill in the [source] block, then `diskwind check {a.path}`")


def _cmd_options(a):
    from .config import WindConfig
    print(WindConfig.describe(a.section))


def _cmd_convert(a):
    from .config import WindConfig
    cfg = WindConfig.from_legacy(a.legacy)
    out = a.out or a.legacy.rsplit(".", 1)[0] + ".toml"
    cfg.save(out)
    print(f"wrote {out}")


def _cmd_center(a):
    from .continuum import fit_disk_center
    c = fit_disk_center(a.continuum, box_arcsec=a.box)
    print(c.summary())
    print("\n# for the [source] block of the configuration:")
    print(c.toml())


def _cmd_check(a):
    from .analysis import WindAnalysis
    wa = WindAnalysis(a.config)
    print(wa.check_noise().summary())
    print(wa.check_disk_rotation().summary())


def _cmd_fit(a):
    from .analysis import WindAnalysis
    wa = WindAnalysis(a.config)
    guide = wa.load_result(a.guide) if a.guide else None
    start = wa.load_result(a.start_from) if a.start_from else None
    fixed = dict(kv.split("=") for kv in (a.fix or []))
    res = wa.fit(a.law, radial=a.radial, rotation=a.rotation, guide=guide, start_from=start,
                 iterations=a.iterations, n_boot=a.nboot, label=a.label,
                 fixed={k: float(v) for k, v in fixed.items()})
    print(res.summary())
    if not a.no_plots:
        wa.plot_cuts(res)
        wa.plot_profiles(res)


def _cmd_run(a):
    from .analysis import WindAnalysis
    WindAnalysis(a.config).run(n_boot=a.nboot, figures=not a.no_plots)


def _cmd_mock(a):
    from .mock import MockWind, make_mock_cube
    m = MockWind.hh30_like() if a.preset == "hh30" else MockWind.elias20_like()
    make_mock_cube(a.fits, m, n_wind=a.n_wind, n_disk=a.n_disk)
    cfg = m.config(a.fits)
    cfg.output.dir = "mock_output"
    cfg.save(a.config)
    print(f"wrote {a.fits} and {a.config}; truth: {m.truth()}")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="diskwind", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("init", help="write a template configuration")
    p.add_argument("path", nargs="?", default="diskwind.toml")
    p.set_defaults(func=_cmd_init)
    p = sub.add_parser("options", help="list all options")
    p.add_argument("section", nargs="?", default=None)
    p.set_defaults(func=_cmd_options)
    p = sub.add_parser("convert", help="convert an old *_config.py")
    p.add_argument("legacy")
    p.add_argument("-o", "--out", default=None)
    p.set_defaults(func=_cmd_convert)
    p = sub.add_parser("center", help="fit the disk center to a continuum image")
    p.add_argument("continuum")
    p.add_argument("--box", type=float, default=None, help="half-size of the fitted region (arcsec)")
    p.set_defaults(func=_cmd_center)
    p = sub.add_parser("check", help="noise and disk-rotation checks")
    p.add_argument("config")
    p.set_defaults(func=_cmd_check)
    p = sub.add_parser("fit", help="one global fit")
    p.add_argument("config")
    p.add_argument("--law", default=None, help="velocity law (default: [fit] velocity_law)")
    p.add_argument("--radial", default=None, choices=("cone", "proportional", "independent"))
    p.add_argument("--rotation", default=None, choices=("constant", "power-law"))
    p.add_argument("--guide", default=None, help="label or path of a result to guide the traces")
    p.add_argument("--start-from", default=None, help="label or path of a result to start from")
    p.add_argument("--iterations", type=int, default=None)
    p.add_argument("--nboot", type=int, default=None)
    p.add_argument("--fix", nargs="*", help="name=value pairs to hold fixed")
    p.add_argument("--label", default=None)
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=_cmd_fit)
    p = sub.add_parser("run", help="the full standard analysis")
    p.add_argument("config")
    p.add_argument("--nboot", type=int, default=None)
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=_cmd_run)
    p = sub.add_parser("mock", help="write a synthetic cube + configuration")
    p.add_argument("--preset", default="elias20", choices=("elias20", "hh30"),
                   help="elias20: inclined, one lobe; hh30: nearly edge-on, both lobes")
    p.add_argument("--fits", default="mock_wind_cube.fits")
    p.add_argument("--config", default="mock.toml")
    p.add_argument("--n-wind", type=int, default=1_500_000)
    p.add_argument("--n-disk", type=int, default=2_000_000)
    p.set_defaults(func=_cmd_mock)
    a = ap.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    sys.exit(main())
