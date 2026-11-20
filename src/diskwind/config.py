"""
Configuration of a wind analysis.

All settings live in one :class:`WindConfig`, organized in sections
(``source``, ``grid``, ``noise``, ``cuts``, ``mask``, ``trace``, ``fit``,
``ring``, ``mass``, ``output``, ``plot``). A configuration can be built in
Python, read from a TOML / YAML / JSON file, or converted from an old-style
``*_config.py`` file of the original scripts.

Examples
--------
Build a configuration in Python and change options with attribute access:

>>> from diskwind import WindConfig
>>> cfg = WindConfig.minimal(fits_file="Elias20_CO.fits", distance_pc=137.5,
...                          mstar_msun=0.75, inclination_deg=49.0,
...                          disk_pa_deg=153.2, wind_pa_deg=243.2)
>>> cfg.trace.snr_min = 4.0
>>> cfg.fit.velocity_law = "broken-power-law"
>>> cfg.save("mysource.toml")                         # doctest: +SKIP

Load a file and see what a cut actually uses after per-cut overrides:

>>> cfg = WindConfig.load("mysource.toml")            # doctest: +SKIP
>>> cfg.for_cut(2.0).trace.snr_min                    # doctest: +SKIP
3.5

Print every option with its default and meaning:

>>> print(WindConfig.describe("trace"))               # doctest: +SKIP

Values that are "unset" (``None`` in Python) are written as ``"auto"`` in
TOML/YAML files, because TOML has no null value.
"""
from __future__ import annotations

import copy
import dataclasses
import difflib
import importlib.util
import json
import os
import warnings
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

__all__ = ["WindConfig", "SourceConfig", "GridConfig", "NoiseConfig", "CutsConfig",
           "MaskConfig", "TraceConfig", "FitConfig", "RingConfig", "MassConfig",
           "OutputConfig", "PlotConfig", "ConfigError"]


class ConfigError(ValueError):
    """Raised for invalid or inconsistent configuration values."""


def opt(default=None, doc="", unit="", choices=None, factory=None):
    """A documented configuration field (used to build the reference tables)."""
    md = dict(doc=doc, unit=unit, choices=choices)
    if factory is not None:
        return field(default_factory=factory, metadata=md)
    return field(default=default, metadata=md)


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------
@dataclass
class SourceConfig:
    """Data file and the fixed geometry of the source (required block)."""
    fits_file: Optional[str] = opt(None, "Path of the line cube (FITS). Relative paths are "
                                   "taken relative to the configuration file.")
    name: str = opt("source", "Name used in titles and output files.")
    distance_pc: Optional[float] = opt(None, "Distance.", "pc")
    mstar_msun: Optional[float] = opt(None, "Stellar mass (Keplerian disk mask, launch "
                                      "radius, escape speed).", "Msun")
    inclination_deg: Optional[float] = opt(None, "Disk inclination, 0 (face-on) to 90 "
                                           "(edge-on). The wind axis is assumed to share it.",
                                           "deg")
    disk_pa_deg: Optional[float] = opt(None, "Position angle (E of N) of the disk major axis.",
                                       "deg")
    wind_pa_deg: Optional[float] = opt(None, "Position angle of the projected wind lobe; "
                                       "must be disk_pa_deg +/- 90.", "deg")
    wind_approaching: Optional[bool] = opt(None, "True if the axis of the main lobe is tilted "
                                           "toward the observer (s = +1), False if away "
                                           "(s = -1); "
                                           "auto: from the mean velocity of the lobe.")
    lobes: str = opt("main", "Lobes to analyze: 'main' (along wind_pa_deg), 'counter' (the "
                     "opposite lobe, along wind_pa_deg + 180) or 'both' (one model fitted to "
                     "both lobes).", choices=("main", "counter", "both"))
    disk_rotation_sense: Optional[int] = opt(None, "+1 if the disk rotates counter-clockwise "
                                             "about the wind axis (+X side red-shifted), -1 "
                                             "otherwise; auto: estimated from the cube.",
                                             choices=(1, -1, None))
    vsys_kms: Optional[float] = opt(None, "Systemic velocity (LSRK); auto: from the line wings "
                                    "of the disk half without wind.", "km/s")
    vsys_err_kms: float = opt(0.1, "Uncertainty of vsys, propagated into Vz.", "km/s")
    vsys_aperture_arcsec: float = opt(1.5, "Aperture radius for the automatic vsys.", "arcsec")
    center_ra: Optional[Any] = opt(None, "RA of the disk center (deg or 'hh:mm:ss'); auto: "
                                   "the reference pixel (CRVAL) of the cube, or the "
                                   "continuum fit if center_from_continuum.")
    center_dec: Optional[Any] = opt(None, "Dec of the disk center (deg or 'dd:mm:ss').")
    center_offset_arcsec: Tuple[float, float] = opt((0.0, 0.0), "Offset (dRA cos dec, dDec) "
                                                    "of the disk center from center_ra/dec (or "
                                                    "from the reference pixel).", "arcsec")
    continuum_fits: Optional[str] = opt(None, "Continuum image (disk center fit, slit-layout "
                                        "figure). Relative to the configuration file.")
    center_from_continuum: bool = opt(False, "Fit the disk center to continuum_fits "
                                      "(elliptical Gaussian) instead of center_ra/dec.")
    restfreq_hz: float = opt(230.538e9, "Rest frequency if not in the header (12CO 2-1).",
                             "Hz")


@dataclass
class GridConfig:
    """Resampling of the cube onto the frame aligned with the wind."""
    x_half_width_arcsec: float = opt(3.0, "Half-length of the transverse cuts (|X| range).",
                                     "arcsec")
    y_range_arcsec: Tuple[float, float] = opt((-1.5, 4.0), "Range along the wind axis to "
                                              "resample (negative: side without wind).",
                                              "arcsec")
    step_arcsec: Optional[float] = opt(None, "Grid step; auto: the native pixel size.",
                                       "arcsec")
    v_range_kms: Optional[Tuple[float, float]] = opt(None, "LSRK velocity range kept in memory; "
                                                     "auto: the whole cube.", "km/s")
    units: str = opt("mJy/beam", "Intensity units of PV diagrams.",
                     choices=("mJy/beam", "Jy/beam", "K"))


@dataclass
class NoiseConfig:
    """How the rms noise of every PV diagram is measured.

    The noise channels are resampled, slit-averaged and smoothed exactly like the
    data, and the rms is 1.4826 x the median absolute deviation."""
    n_edge_channels: int = opt(10, "Line-free channels at EACH end of the cube.")
    vel_ranges_kms: Optional[List[Tuple[float, float]]] = opt(
        None, "LSRK velocity ranges of the noise channels (overrides n_edge_channels), "
        "e.g. [[-12, -8]].", "km/s")
    rms_jy_beam: Optional[float] = opt(None, "Force the channel rms instead of measuring it.",
                                       "Jy/beam")


@dataclass
class CutsConfig:
    """Positions of the transverse cuts and per-cut overrides."""
    y_offsets_arcsec: Any = opt(None, "Offsets of the cuts along the wind axis: a list, or "
                                "{start, stop, step}.", "arcsec")
    slit_width_arcsec: Optional[float] = opt(None, "Width of a cut along Y; auto: half the "
                                             "beam FWHM projected on Y.", "arcsec")
    overrides: Dict[str, Dict[str, Any]] = opt(
        doc="Per-cut overrides: {y0: {option: value}}, e.g. {'2.00': {'snr_min': 3.5, "
        "'slit_width_arcsec': 0.1}}. Option names may be bare (unique) or 'section.name'.",
        factory=dict)


@dataclass
class MaskConfig:
    """Which PV pixels may be used by the trace (contamination masks)."""
    vrel_range_kms: Optional[Tuple[float, float]] = opt(None, "Window of V - Vsys kept for "
                                                        "the trace; auto: all.", "km/s")
    v_exclude_kms: List[Tuple[float, float]] = opt(doc="Ranges of V - Vsys excluded (e.g. "
                                                   "cloud absorption).", unit="km/s",
                                                   factory=list)
    x_range_arcsec: Optional[Tuple[float, float]] = opt(None, "Fixed window in X kept for the "
                                                        "trace.", "arcsec")
    x_cone: Optional[Tuple[float, ...]] = opt(None, "Window that follows the wind: (xc, w0, "
                                              "slope[, xc_slope]) -> |X - (xc + xc_slope y)| "
                                              "<= w0 + slope y. Takes priority over "
                                              "x_range_arcsec.", "arcsec")
    x_max_arcsec: float = opt(2.5, "|X| limit when no window is given.", "arcsec")
    kepler: bool = opt(True, "Mask channels where Keplerian disk gas can emit.")
    kepler_dv_kms: float = opt(0.8, "Margin around the Keplerian velocities.", "km/s")
    disk_rout_au: float = opt(250.0, "Outer radius of the line-emitting disk.", "au")
    disk_z_over_r: float = opt(0.35, "Height z/r of the emitting surfaces (gas between "
                                     "them is masked).")
    disk_height_func: Any = opt(None, "Python only: function r_au -> H_au overriding "
                                "disk_z_over_r.")


@dataclass
class TraceConfig:
    """Extraction of the trace (locus of the emission) in each PV diagram."""
    mode: str = opt("wings", "'wings': blue and red wing fitted separately (full rings, "
                    "Louvet+2018); 'single': one peak per position (partial rings).",
                    choices=("wings", "single"))
    snr_min: float = opt(5.0, "Minimum peak S/N of a trace point.")
    fit_frac: float = opt(0.3, "Gaussian fitted to channels above this fraction of the peak.")
    err_inflate: float = opt(2.0, "Velocity error = this x the formal centroid error.")
    dx_arcsec: Optional[float] = opt(None, "Spacing of trace points along the cut; auto: "
                                     "half the beam FWHM along X (about two points per "
                                     "beam: faster fits, honest chi2); 0: every pixel.",
                                     "arcsec")
    min_points: int = opt(8, "Cuts with fewer trace points (counted every pixel, before "
                          "thinning to dx_arcsec) are skipped.")
    sigx_beam: float = opt(0.5, "Position error of a point, in units of the beam dispersion "
                                "along X (fit statistic).")
    drop_clipped: bool = opt(True, "Drop points whose line profile is cut by a mask edge.")
    guide_dv_kms: Any = opt(1.0, "Guided extraction: channels within +-dv of the model, or "
                            "(below, above).", "km/s")
    peak_choice: str = opt("brightest", "Guided extraction: trace the brightest allowed peak "
                           "or the significant peak nearest to the model.",
                           choices=("brightest", "nearest"))
    xsmooth_arcsec: Optional[float] = opt(None, "Boxcar smoothing of the PV along X (faint "
                                          "cuts).", "arcsec")
    vsmooth_chan: Optional[int] = opt(None, "Boxcar smoothing of the PV in velocity "
                                      "(channels).")


@dataclass
class FitConfig:
    """Global fit of one conical wind to all cuts."""
    velocity_law: str = opt("power-law", "Axial velocity law Vz(z); see "
                            "diskwind.laws.available_laws().")
    radial: str = opt("cone", "Radial velocity: 'cone' (Vr = Vz tan theta, flow along the "
                      "cone wall), 'proportional' (Vr = free amplitude x same law), "
                      "'independent' (own law).",
                      choices=("cone", "proportional", "independent"))
    rotation: str = opt("constant", "Specific angular momentum: 'constant' (j) or "
                        "'power-law' (j (z/zref)^qj).", choices=("constant", "power-law"))
    fixed: Dict[str, float] = opt(doc="Parameters held fixed, {name: value}.", factory=dict)
    bounds: Dict[str, Tuple[float, float]] = opt(doc="Bounds, {name: (lo, hi)}.",
                                                 factory=dict)
    equal_cut_weight: bool = opt(False, "Give every cut the same total weight "
                                 "(sqrt(<n>/n_k)).")
    loss: str = opt("linear", "scipy least_squares loss.",
                    choices=("linear", "soft_l1", "huber", "cauchy", "arctan"))
    n_boot: int = opt(100, "Bootstrap resamplings of the cuts (0: none). With "
                      "error_method = 'jackknife', any value > 0 switches it on.")
    error_method: str = opt("bootstrap", "Resampling of the cuts for the errors: "
                            "'bootstrap' (n_boot draws with replacement) or 'jackknife' "
                            "(one refit per left-out cut: far fewer refits).",
                            choices=("bootstrap", "jackknife"))
    n_jobs: int = opt(-1, "Processes for the resampling refits (-1: all cores; 1: serial). "
                      "Falls back to serial if the parallel run fails.")
    resample_tol: float = opt(1e-6, "Tolerance of the resampling refits (they start at "
                              "the best fit); the main fit uses 1e-8.")
    n_refine: int = opt(4, "Multistart: every starting point is solved coarsely, and "
                        "only this many of the best are solved to full precision (both "
                        "steps in parallel on n_jobs processes).")
    check_vr_sign: Any = opt(True, "If Vr is free and comes out < 0, refit with Vr >= 0 and "
                             "keep the better (True) or always Vr >= 0 ('prefer_positive').")
    scale_errors_by_chi2: bool = opt(True, "Scale covariance errors by chi2_red when > 1.")
    zref_au: float = opt(100.0, "Reference height of power laws.", "au")
    bpl_delta: float = opt(0.1, "Smoothness of the broken power law.")
    iterations: int = opt(1, "Guided iterations: re-extract traces around the new model "
                          "and refit.")
    j_multistart: bool = opt(True, "Restart the fit with j = +-j_s, +-3 j_s.")
    max_nfev: int = opt(3000, "Maximum function evaluations per start.")
    nphi: int = opt(181, "Azimuths per model curve (every 2 deg).")
    seed: int = opt(42, "Random seed (bootstrap, Monte Carlo).")


@dataclass
class RingConfig:
    """Independent ring fits to each cut (Louvet+2018 Sect. 4.4); best for full rings."""
    rotation_law: str = opt("j", "Rotation across one cut: 'j' (V_phi ~ 1/R) or 'const'.",
                            choices=("j", "const"))
    cone_half_angle_deg: Optional[float] = opt(None, "Fix the cone angle; auto: two-pass "
                                               "estimate.", "deg")
    n_passes: int = opt(2, "Passes (fit, fit R(z), refit with theta fixed).")
    n_mc: int = opt(1000, "Monte Carlo draws for the derived quantities.")
    min_phi_coverage: float = opt(0.6, "Rings with less azimuthal coverage are flagged.")
    cone_fit_zrange_au: Tuple[float, float] = opt((0.0, float("inf")), "Heights used in the "
                                                  "R(z) fit.", "au")
    prior_sigma_kms: float = opt(0.5, "Width of priors from a global fit.", "km/s")


@dataclass
class MassConfig:
    """Brightness, mass and mass-loss rate in slices (Louvet+2018 Sect. 4.2)."""
    tex_k: float = opt(30.0, "Excitation temperature.", "K")
    tex_range_k: Tuple[float, float] = opt((20.0, 50.0), "T_ex range for the systematic band.",
                                           "K")
    xco: float = opt(1e-4, "CO abundance relative to H2.")
    mu: float = opt(1.4, "Mean mass per H2 in units of m_H2.")
    dz_au: float = opt(30.0, "Slice thickness along the axis.", "au")
    nsigma: float = opt(5.0, "Pixels with W > nsigma sigma_W are kept.")
    flux_cal: float = opt(0.10, "Fractional flux-calibration uncertainty.")
    nmc: int = opt(1000, "Monte Carlo draws for the totals.")
    y_range_arcsec: Optional[Tuple[float, float]] = opt(None, "Range of the slices; auto: the "
                                                        "cuts.", "arcsec")
    nu_hz: float = opt(230.538e9, "Line frequency (12CO 2-1).", "Hz")
    a_ul: float = opt(6.910e-7, "Einstein A of the line.", "s^-1")
    g_up: float = opt(5.0, "Upper-level degeneracy.")
    e_up_k: float = opt(16.596, "Upper-level energy E_up/k.", "K")
    b_rot_hz: float = opt(57.635968e9, "Rotational constant (partition function).", "Hz")


@dataclass
class OutputConfig:
    """Where results go."""
    dir: str = opt("diskwind_output", "Output folder (created if needed).")
    fig_format: str = opt("pdf", "Figure format.", choices=("pdf", "png", "svg", "jpg"))
    dpi: int = opt(200, "Figure resolution.")
    verbose: bool = opt(True, "Print progress messages.")


@dataclass
class PlotConfig:
    """Defaults for figures."""
    xlim: Optional[Tuple[float, float]] = opt(None, "X range (left, right).", "arcsec")
    vlim: Optional[Tuple[float, float]] = opt(None, "Velocity range.", "km/s")
    contour_start_sigma: float = opt(5.0, "First contour.")
    contour_step_sigma: float = opt(5.0, "Contour step.")
    cmap: str = opt("magma", "Color map of PV diagrams.")
    mom0_fits: Optional[str] = opt(None, "Moment-0 image for the slit-layout figure; auto: "
                                   "integrated from the cube.")
    map_cmap: str = opt("cmr.freeze", "Color map of the moment-0 image (cmasher maps need "
                        "the cmasher package; falls back to 'bone').")
    line_label: Optional[str] = opt(None, "Name of the line shown on maps, e.g. "
                                    "'$^{12}$CO $J$=2-1'.")
    map_min: Optional[float] = opt(None, "Floor of the moment-0 image (black below); auto: "
                                   "its median.", "mJy/beam km/s")
    continuum_min: Optional[float] = opt(None, "Continuum shown above this level; auto: its "
                                         "90th percentile.", "mJy/beam")
    pv_grid_offsets_arcsec: Any = opt(None, "Cuts of the PV grid figure; auto: the fit cuts.",
                                      "arcsec")
    ncols: int = opt(4, "Columns of multi-panel figures.")


SECTIONS = {"source": SourceConfig, "grid": GridConfig, "noise": NoiseConfig,
            "cuts": CutsConfig, "mask": MaskConfig, "trace": TraceConfig, "fit": FitConfig,
            "ring": RingConfig, "mass": MassConfig, "output": OutputConfig,
            "plot": PlotConfig}
REQUIRED = ["source.fits_file", "source.distance_pc", "source.mstar_msun",
            "source.inclination_deg", "source.disk_pa_deg", "source.wind_pa_deg"]


def _field_index():
    idx = {}
    for sec, cls in SECTIONS.items():
        for f in dataclasses.fields(cls):
            idx.setdefault(f.name, []).append(sec)
    return idx


_FIELDS = _field_index()


def _resolve(key):
    """'snr_min' or 'trace.snr_min' -> ('trace', 'snr_min')."""
    if "." in key:
        sec, name = key.split(".", 1)
        if sec not in SECTIONS or name not in {f.name for f in dataclasses.fields(SECTIONS[sec])}:
            raise ConfigError(f"unknown option {key!r}{_suggest(key, _all_keys())}")
        return sec, name
    if key not in _FIELDS:
        raise ConfigError(f"unknown option {key!r}{_suggest(key, list(_FIELDS))}")
    if len(_FIELDS[key]) > 1:
        raise ConfigError(f"option {key!r} is ambiguous: use one of "
                          f"{[s + '.' + key for s in _FIELDS[key]]}")
    return _FIELDS[key][0], key


def _all_keys():
    return [f"{s}.{f.name}" for s, c in SECTIONS.items() for f in dataclasses.fields(c)]


def _suggest(word, choices):
    m = difflib.get_close_matches(word, choices, n=2)
    return f" (did you mean {' or '.join(repr(x) for x in m)}?)" if m else ""


def _none(v):
    """'auto' / 'none' / '' in a file mean None."""
    if isinstance(v, str) and v.strip().lower() in ("auto", "none", "null", ""):
        return None
    return v


def _clean(v):
    """Recursively convert file values (lists, 'auto') to Python values."""
    v = _none(v)
    if isinstance(v, dict):
        return {k: _clean(x) for k, x in v.items()}
    if isinstance(v, list):
        return [_clean(x) for x in v]
    return v


def _plain(v):
    """Python value -> something JSON/TOML can hold (None -> 'auto')."""
    if v is None:
        return "auto"
    if isinstance(v, (np.floating, np.integer)):
        return v.item()
    if isinstance(v, np.ndarray):
        return [_plain(x) for x in v.tolist()]
    if isinstance(v, (list, tuple)):
        return [_plain(x) for x in v]
    if isinstance(v, dict):
        return {str(k): _plain(x) for k, x in v.items()}
    if isinstance(v, float) and not np.isfinite(v):
        return "inf" if v > 0 else "-inf"
    if callable(v):
        return "auto"
    return v


# ---------------------------------------------------------------------------
# The configuration
# ---------------------------------------------------------------------------
@dataclass
class WindConfig:
    """All settings of a wind analysis, in sections.

    Attributes
    ----------
    source, grid, noise, cuts, mask, trace, fit, ring, mass, output, plot
        The sections (see :meth:`describe` for every option).
    path : str or None
        File the configuration was read from (relative paths are resolved
        against its folder).

    Examples
    --------
    >>> cfg = WindConfig()                 # all defaults
    >>> cfg.trace.mode = "single"          # partial rings
    >>> _ = cfg.set(snr_min=4, **{"fit.n_boot": 0})   # returns the config
    >>> cfg.get("snr_min")
    4
    """
    source: SourceConfig = field(default_factory=SourceConfig)
    grid: GridConfig = field(default_factory=GridConfig)
    noise: NoiseConfig = field(default_factory=NoiseConfig)
    cuts: CutsConfig = field(default_factory=CutsConfig)
    mask: MaskConfig = field(default_factory=MaskConfig)
    trace: TraceConfig = field(default_factory=TraceConfig)
    fit: FitConfig = field(default_factory=FitConfig)
    ring: RingConfig = field(default_factory=RingConfig)
    mass: MassConfig = field(default_factory=MassConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    plot: PlotConfig = field(default_factory=PlotConfig)
    path: Optional[str] = None

    # ----------------------------------------------------------- construction
    @classmethod
    def minimal(cls, fits_file, distance_pc, mstar_msun, inclination_deg, disk_pa_deg,
                wind_pa_deg, **options):
        """Configuration from the required source values plus any options.

        Parameters
        ----------
        fits_file, distance_pc, mstar_msun, inclination_deg, disk_pa_deg, wind_pa_deg
            The required values (see :class:`SourceConfig`).
        **options
            Any other option, by bare or dotted name (``snr_min=4``,
            ``**{"fit.n_boot": 0}``).

        Returns
        -------
        WindConfig
        """
        cfg = cls()
        cfg.source = SourceConfig(fits_file=fits_file, distance_pc=distance_pc,
                                  mstar_msun=mstar_msun, inclination_deg=inclination_deg,
                                  disk_pa_deg=disk_pa_deg, wind_pa_deg=wind_pa_deg)
        cfg.set(**options)
        return cfg

    @classmethod
    def from_dict(cls, d, path=None):
        """Build from a nested dict ``{section: {option: value}}``.

        Unknown sections or options raise :class:`ConfigError` with a
        suggestion; ``"auto"`` means unset (None).
        """
        cfg = cls(path=path)
        d = _clean(dict(d))
        for sec, vals in d.items():
            if sec in ("path",):
                continue
            if sec not in SECTIONS:
                raise ConfigError(f"unknown section [{sec}]{_suggest(sec, list(SECTIONS))}")
            if not isinstance(vals, dict):
                raise ConfigError(f"section [{sec}] must be a table of options")
            obj = getattr(cfg, sec)
            names = {f.name for f in dataclasses.fields(obj)}
            for k, v in vals.items():
                if k not in names:
                    raise ConfigError(f"unknown option {sec}.{k}{_suggest(k, sorted(names))}")
                setattr(obj, k, v)
        cfg._normalize()
        return cfg

    @classmethod
    def load(cls, path):
        """Read a configuration file.

        Parameters
        ----------
        path : str
            ``.toml`` (recommended), ``.yaml``/``.yml`` (needs PyYAML), ``.json``,
            or ``.py``. A ``.py`` file either defines ``CONFIG = {...}`` (nested
            dict) or is an old-style configuration with upper-case keys
            (``FITS_FILE = ...``), which is converted with :meth:`from_legacy`.

        Returns
        -------
        WindConfig
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"configuration file not found: {path}")
        ext = os.path.splitext(path)[1].lower()
        if ext == ".toml":
            try:
                import tomllib
            except ImportError:            # Python < 3.11
                import tomli as tomllib
            with open(path, "rb") as fh:
                d = tomllib.load(fh)
        elif ext in (".yaml", ".yml"):
            try:
                import yaml
            except ImportError as exc:
                raise ImportError("reading YAML needs PyYAML (pip install pyyaml); "
                                  "or use a .toml file") from exc
            with open(path) as fh:
                d = yaml.safe_load(fh) or {}
        elif ext == ".json":
            with open(path) as fh:
                d = json.load(fh)
        elif ext == ".py":
            mod = _exec_py(path)
            if hasattr(mod, "CONFIG"):
                d = mod.CONFIG
            else:
                return cls.from_legacy(path)
        else:
            raise ConfigError(f"unsupported configuration format {ext!r}")
        return cls.from_dict(d, path=os.path.abspath(path))

    @classmethod
    def from_legacy(cls, path):
        """Convert an old-style ``*_config.py`` (``FITS_FILE = ...`` etc.).

        Every upper-case key of the original scripts is mapped onto the new
        options (see :data:`diskwind.config.LEGACY_KEYS`); unknown keys are
        reported. Save the result with :meth:`save` to switch to the new format.
        """
        mod = _exec_py(path)
        old = {k: getattr(mod, k) for k in dir(mod) if k.isupper()}
        cfg = cls(path=os.path.abspath(path))
        unknown = []
        for k, v in old.items():
            if k in _LEGACY_SPECIAL:
                continue
            tgt = LEGACY_KEYS.get(k)
            if tgt is None:
                unknown.append(k)
                continue
            if tgt == "":
                continue
            sec, name = tgt.split(".")
            setattr(getattr(cfg, sec), name, v)
        # special conversions
        ov = {}
        for y, d in (old.get("TRACE_SLIT_OVERRIDES") or {}).items():
            new = {}
            for k, v in d.items():
                if LEGACY_KEYS.get(k):
                    new[LEGACY_KEYS[k]] = v
                else:
                    unknown.append(f"TRACE_SLIT_OVERRIDES[{y}].{k}")
            ov[f"{float(y):.2f}"] = new
        cfg.cuts.overrides = ov
        vm = old.get("VELOCITY_MODEL", "single")
        flow = bool(old.get("GLOBAL_FLOW_ALONG_CONE", False))
        fixed = dict(old.get("GLOBAL_FIXED", {"xslope": 0.0, "qv": 0.0}) or {})
        if vm == "single":
            cfg.fit.velocity_law = "power-law"
            cfg.fit.radial = "cone" if flow else "proportional"
        else:
            cfg.fit.velocity_law = "broken-power-law"
            cfg.fit.radial = "independent" if vm == "bpl-separate" else "cone"
        if "qj" in fixed and float(fixed["qj"]) == 0.0:
            cfg.fit.rotation = "constant"
            fixed.pop("qj")
        else:
            cfg.fit.rotation = "power-law"
        if cfg.fit.radial == "cone":
            fixed.pop("Vr", None)
        cfg.fit.fixed = fixed
        if "WIND_APPROACHING" in old:
            cfg.source.wind_approaching = old["WIND_APPROACHING"]
        if unknown:
            warnings.warn(f"legacy keys not converted: {sorted(unknown)}")
        cfg._normalize()
        return cfg

    # ----------------------------------------------------------------- output
    def to_dict(self, plain=True):
        """Nested dict of all options (``plain``: file-safe values, None -> 'auto')."""
        out = {}
        for sec in SECTIONS:
            obj = getattr(self, sec)
            d = {f.name: getattr(obj, f.name) for f in dataclasses.fields(obj)}
            if sec == "mask":
                d.pop("disk_height_func")
            out[sec] = {k: _plain(v) for k, v in d.items()} if plain else d
        return out

    def save(self, path):
        """Write the configuration (``.toml``, ``.json`` or ``.yaml``)."""
        ext = os.path.splitext(path)[1].lower()
        if self.mask.disk_height_func is not None:
            warnings.warn("mask.disk_height_func is a Python function and is not saved")
        d = self.to_dict()
        if ext == ".toml":
            txt = _to_toml(d)
        elif ext == ".json":
            txt = json.dumps(d, indent=2)
        elif ext in (".yaml", ".yml"):
            import yaml
            txt = yaml.safe_dump(d, sort_keys=False)
        else:
            raise ConfigError(f"unsupported format {ext!r}")
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(txt)
        return path

    @staticmethod
    def describe(section=None, fmt="text"):
        """Reference of all options: name, default, unit and meaning.

        Parameters
        ----------
        section : str, optional
            Only this section.
        fmt : {'text', 'rst'}
            Plain text, or a reStructuredText list-table (used by the docs).
        """
        secs = [section] if section else list(SECTIONS)
        lines = []
        for sec in secs:
            cls = SECTIONS[sec]
            doc = (cls.__doc__ or "").strip().splitlines()[0]
            if fmt == "rst":
                lines += [f".. _config-{sec}:", "", f"``[{sec}]``", "~" * (len(sec) + 6), "",
                          (cls.__doc__ or "").strip(), "",
                          ".. list-table::", "   :header-rows: 1", "   :widths: 22 14 8 56", "",
                          "   * - Option", "     - Default", "     - Unit", "     - Meaning"]
            else:
                lines += [f"[{sec}]  {doc}"]
            for f in dataclasses.fields(cls):
                d = f.default if f.default is not dataclasses.MISSING else f.default_factory()
                dv = "auto" if d is None else repr(d)
                md = f.metadata
                txt = md.get("doc", "")
                if md.get("choices"):
                    txt += " Choices: " + ", ".join("auto" if c is None else repr(c)
                                                    for c in md["choices"]) + "."
                if fmt == "rst":
                    lines += [f"   * - ``{f.name}``", f"     - ``{dv}``",
                              f"     - {md.get('unit') or '-'}", f"     - {txt}"]
                else:
                    u = f" [{md['unit']}]" if md.get("unit") else ""
                    lines.append(f"  {f.name:24s} = {dv:14s}{u}  {txt}")
            lines.append("")
        return "\n".join(lines)

    # --------------------------------------------------------------- access
    def set(self, **options):
        """Set options by bare (``snr_min``) or dotted (``trace.snr_min``) name."""
        for k, v in options.items():
            sec, name = _resolve(k)
            setattr(getattr(self, sec), name, v)
        self._normalize()
        return self

    def get(self, key):
        """Value of an option by bare or dotted name."""
        sec, name = _resolve(key)
        return getattr(getattr(self, sec), name)

    def copy(self):
        """Deep copy."""
        return copy.deepcopy(self)

    def for_cut(self, y0):
        """The configuration seen by the cut at ``y0`` (overrides applied).

        A per-cut ``x_range_arcsec`` replaces the ``x_cone`` window for that cut.
        """
        ov = self.override_for(y0)
        if not ov:
            return self
        c = copy.copy(self)
        for sec in SECTIONS:
            setattr(c, sec, copy.copy(getattr(self, sec)))
        for k, v in ov.items():
            sec, name = _resolve(k)
            setattr(getattr(c, sec), name, v)
            if name == "x_range_arcsec":
                c.mask.x_cone = None
        return c

    def override_for(self, y0):
        """The override dict of the cut at ``y0`` (or {})."""
        for k, d in (self.cuts.overrides or {}).items():
            if abs(float(k) - float(y0)) < 1e-3:
                return d
        return {}

    @property
    def y_offsets(self):
        """Cut offsets as a sorted numpy array (arcsec)."""
        return _offsets(self.cuts.y_offsets_arcsec)

    @property
    def output_dir(self):
        """Absolute output folder (created on first use)."""
        d = self.output.dir
        if not os.path.isabs(d) and self.path:
            d = os.path.join(os.path.dirname(self.path), d)
        os.makedirs(d, exist_ok=True)
        return d

    def outpath(self, name, ext=None):
        """Path ``<output_dir>/<name>.<ext>`` (default ext: the figure format)."""
        return os.path.join(self.output_dir, f"{name}.{ext or self.output.fig_format}")

    @property
    def fits_path(self):
        """The cube path, resolved against the configuration's folder."""
        f = self.source.fits_file
        if f is None:
            raise ConfigError("source.fits_file is not set")
        if not os.path.isabs(f) and not os.path.exists(f) and self.path:
            alt = os.path.join(os.path.dirname(self.path), f)
            if os.path.exists(alt):
                return alt
            raise ConfigError(f"cube not found: {f!r} (looked in the current folder and in "
                              f"{os.path.dirname(os.path.abspath(self.path))}); put it there "
                              "(e.g. examples/data/<source>/) or set source.fits_file")
        return f

    def resolve(self, f):
        """A file path resolved like ``fits_file`` (relative to the configuration's
        folder if not found from the current folder); None stays None."""
        if f is None or os.path.isabs(f) or os.path.exists(f) or not self.path:
            return f
        alt = os.path.join(os.path.dirname(self.path), f)
        return alt if os.path.exists(alt) else f

    def validate(self, need_data=True):
        """Check required values and consistency; returns a list of warnings."""
        miss = [k for k in REQUIRED if self.get(k) is None and (need_data or "fits" not in k)]
        if miss:
            raise ConfigError(f"missing required options: {miss}")
        msgs = []
        i = self.source.inclination_deg
        if not 0.0 <= i <= 90.0:
            raise ConfigError("inclination_deg must be in [0, 90]: for i > 90 use 180 - i "
                              "and set wind_approaching accordingly")
        dpa = (self.source.wind_pa_deg - self.source.disk_pa_deg) % 180.0
        if abs(dpa - 90.0) > 5.0:
            msgs.append(f"wind_pa - disk_pa = {dpa:.1f} deg (mod 180), not ~90: the wind axis "
                        "is not along the disk minor axis")
        if self.source.lobes == "both" and (self.source.vsys_kms is None
                                            or self.source.disk_rotation_sense is None):
            msgs.append("lobes = 'both': the automatic Vsys and disk-rotation estimates use "
                        "the side of the disk without wind; set source.vsys_kms and "
                        "source.disk_rotation_sense")
        from .laws import available_laws
        if self.fit.velocity_law not in available_laws():
            raise ConfigError(f"unknown velocity_law {self.fit.velocity_law!r}; available: "
                              f"{available_laws()}")
        for sec in SECTIONS:
            obj = getattr(self, sec)
            for f in dataclasses.fields(obj):
                ch = f.metadata.get("choices")
                v = getattr(obj, f.name)
                if ch and v not in ch and not (f.name == "check_vr_sign"):
                    raise ConfigError(f"{sec}.{f.name} = {v!r}: must be one of {ch}")
        for m in msgs:
            warnings.warn(m)
        return msgs

    def _normalize(self):
        s = self.source
        if s.disk_rotation_sense is not None:
            s.disk_rotation_sense = int(s.disk_rotation_sense)
        if isinstance(self.cuts.overrides, dict):
            self.cuts.overrides = {f"{float(k):.2f}": dict(v)
                                   for k, v in self.cuts.overrides.items()}
            for d in self.cuts.overrides.values():
                for k in d:
                    _resolve(k)
        r = self.ring.cone_fit_zrange_au
        if r is not None:
            self.ring.cone_fit_zrange_au = tuple(float(x) for x in r)

    def __repr__(self):
        s = self.source
        return (f"WindConfig(source={s.name!r}, fits={s.fits_file!r}, i={s.inclination_deg}, "
                f"d={s.distance_pc} pc, law={self.fit.velocity_law!r}, "
                f"trace.mode={self.trace.mode!r}, output={self.output.dir!r})")


def _offsets(v):
    if v is None:
        return None
    if isinstance(v, dict):
        a = np.arange(float(v["start"]), float(v["stop"]) + 0.5 * float(v["step"]),
                      float(v["step"]))
        return np.round(a, 4)
    return np.sort(np.atleast_1d(np.asarray(v, float)))


def _exec_py(path):
    spec = importlib.util.spec_from_file_location("diskwind_user_config", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _toml_val(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, str):
        return json.dumps(v)
    if isinstance(v, list):
        return "[" + ", ".join(_toml_val(x) for x in v) + "]"
    if isinstance(v, dict):
        return "{" + ", ".join(f"{_toml_key(k)} = {_toml_val(x)}" for k, x in v.items()) + "}"
    raise TypeError(f"cannot write {type(v)} to TOML")


def _toml_key(k):
    k = str(k)
    return k if k.replace("_", "").replace("-", "").isalnum() else json.dumps(k)


def _to_toml(d):
    lines = ["# diskwind configuration (values 'auto' = unset / automatic)", ""]
    for sec, vals in d.items():
        lines.append(f"[{sec}]")
        sub = []
        for k, v in vals.items():
            if sec == "cuts" and k == "overrides" and isinstance(v, dict):
                sub = [(f"cuts.overrides.{_toml_key(y)}", o) for y, o in v.items()]
                continue
            lines.append(f"{k} = {_toml_val(v)}")
        lines.append("")
        for name, o in sub:
            lines.append(f"[{name}]")
            lines += [f"{_toml_key(k)} = {_toml_val(x)}" for k, x in o.items()]
            lines.append("")
    return "\n".join(lines)


#: Mapping from the upper-case keys of the original scripts to the new options
#: (``""``: no longer used).
LEGACY_KEYS = {
    "FITS_FILE": "source.fits_file", "SOURCE_NAME": "source.name",
    "OUTPUT_DIR": "output.dir", "RESTFREQ_HZ": "source.restfreq_hz", "UNITS": "grid.units",
    "FIG_FORMAT": "output.fig_format", "FIG_DPI": "output.dpi",
    "DISTANCE_PC": "source.distance_pc", "MSTAR_MSUN": "source.mstar_msun",
    "INCLINATION_DEG": "source.inclination_deg", "DISK_PA_DEG": "source.disk_pa_deg",
    "WIND_PA_DEG": "source.wind_pa_deg", "DISK_ROTATION_SENSE": "source.disk_rotation_sense",
    "VSYS_KMS": "source.vsys_kms", "VSYS_ERR_KMS": "source.vsys_err_kms",
    "VSYS_APERTURE_ARCSEC": "source.vsys_aperture_arcsec", "CENTER_RA_DEG": "source.center_ra",
    "CENTER_DEC_DEG": "source.center_dec", "CENTER_OFFSET_ARCSEC": "source.center_offset_arcsec",
    "X_HALF_WIDTH_ARCSEC": "grid.x_half_width_arcsec", "Y_RANGE_ARCSEC": "grid.y_range_arcsec",
    "GRID_STEP_ARCSEC": "grid.step_arcsec", "V_RANGE_KMS": "grid.v_range_kms",
    "NOISE_CHANNELS": "noise.n_edge_channels", "NOISE_VEL_RANGES_KMS": "noise.vel_ranges_kms",
    "RMS_JY_BEAM": "noise.rms_jy_beam", "SLIT_WIDTH_ARCSEC": "cuts.slit_width_arcsec",
    "V_EXCLUDE_KMS": "mask.v_exclude_kms", "USE_KEPLER_MASK": "mask.kepler",
    "KEPLER_MASK_DV_KMS": "mask.kepler_dv_kms", "DISK_ROUT_AU": "mask.disk_rout_au",
    "DISK_MAX_Z_OVER_R": "mask.disk_z_over_r", "DISK_HEIGHT_FUNC": "mask.disk_height_func",
    "TRACE_SNR_MIN": "trace.snr_min", "TRACE_FIT_FRAC": "trace.fit_frac",
    "TRACE_ERR_INFLATE": "trace.err_inflate", "TRACE_XMAX_ARCSEC": "mask.x_max_arcsec",
    "TRACE_DX_ARCSEC": "trace.dx_arcsec", "TRACE_MIN_POINTS": "trace.min_points",
    "TRACE_SIGX_BEAM": "trace.sigx_beam", "MIN_PHI_COVERAGE": "ring.min_phi_coverage",
    "TRACE_MODE": "trace.mode", "TRACE_X_RANGE_ARCSEC": "mask.x_range_arcsec",
    "TRACE_VREL_RANGE_KMS": "mask.vrel_range_kms", "TRACE_DROP_CLIPPED": "trace.drop_clipped",
    "TRACE_X_CONE": "mask.x_cone", "TRACE_GUIDE_DV_KMS": "trace.guide_dv_kms",
    "TRACE_XSMOOTH_ARCSEC": "trace.xsmooth_arcsec", "TRACE_VSMOOTH_CHAN": "trace.vsmooth_chan",
    "TRACE_PEAK_CHOICE": "trace.peak_choice", "GLOBAL_EQUAL_SLIT_WEIGHT": "fit.equal_cut_weight",
    "GLOBAL_N_BOOT": "fit.n_boot", "GLOBAL_BOUNDS": "fit.bounds",
    "GLOBAL_CHECK_VR_SIGN": "fit.check_vr_sign", "FIT_Y_OFFSETS_ARCSEC": "cuts.y_offsets_arcsec",
    "CONE_HALF_ANGLE_DEG": "ring.cone_half_angle_deg", "FIT_N_PASSES": "ring.n_passes",
    "CONE_FIT_ZRANGE_AU": "ring.cone_fit_zrange_au", "N_MC": "ring.n_mc",
    "SCALE_ERRORS_BY_CHI2": "fit.scale_errors_by_chi2", "FIT_LOSS": "fit.loss",
    "RANDOM_SEED": "fit.seed", "ROTATION_LAW": "ring.rotation_law",
    "PRIOR_SIGMA_KMS": "ring.prior_sigma_kms", "FIG7_XLIM": "plot.xlim",
    "FIG7_VLIM": "plot.vlim", "CONTOUR_START_SIGMA": "plot.contour_start_sigma",
    "CONTOUR_STEP_SIGMA": "plot.contour_step_sigma", "FIG7_CMAP": "plot.cmap",
    "FIG7_Y_OFFSETS_ARCSEC": "plot.pv_grid_offsets_arcsec", "FIG7_NCOLS": "plot.ncols",
    # no longer used / handled elsewhere
    "GLOBAL_Y_OFFSETS_ARCSEC": "", "RING_FIT_MODE": "", "RING_GUIDED": "", "FLIP_X_AXIS": "",
    "FIG8_Y_ARCSEC": "", "FIG8_USE_FIT": "", "FIG8_SHELLS": "", "FIG8_GAMMA": "",
    "FIG8_CONTOUR_SIGMA": "", "FIG11_Y_ARCSEC": "", "FIG11_N_DRAW": "",
}
_LEGACY_SPECIAL = {"TRACE_SLIT_OVERRIDES", "VELOCITY_MODEL", "GLOBAL_FLOW_ALONG_CONE",
                   "GLOBAL_FIXED", "WIND_APPROACHING"}
