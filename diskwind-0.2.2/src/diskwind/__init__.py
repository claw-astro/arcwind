"""
diskwind -- geometry and kinematics of disk winds from transverse
position-velocity diagrams.

The method of Louvet et al. (2018, A&A 618, A120) for the edge-on HH 30
outflow, generalized to any inclination, to partially detected rings and to
flexible velocity laws: a conical wind is fitted to the traces of many cuts at
once.

Quick start
-----------
>>> from diskwind import WindAnalysis
>>> wa = WindAnalysis("mysource.toml")                       # doctest: +SKIP
>>> res = wa.fit()                                           # doctest: +SKIP
>>> print(res.summary())                                     # doctest: +SKIP

Main classes
------------
:class:`WindConfig`        all settings (TOML/YAML/JSON/Python)
:class:`WindAnalysis`      the analysis, step by step or ``run()``
:class:`WindCube`          the resampled cube and its PV diagrams
:class:`TraceExtractor`    traces of the emission in the PV diagrams
:class:`ConicalWindModel`  the model (geometry + velocity laws)
:class:`GlobalFitter`      fit of the model to all cuts; returns :class:`FitResult`
:class:`RingFitter`        independent per-cut ring fits (full rings)
:class:`MassAnalysis`      mass and mass-loss rate
:func:`fit_disk_center`    disk center from a continuum image
"""
from ._version import __version__
from .analysis import WindAnalysis
from .config import WindConfig
from .continuum import fit_disk_center
from .cube import WindCube
from .derived import WindProperties
from .fitting import FitResult, GlobalFitter, compare_fits
from .geometry import ViewGeometry
from .laws import VelocityLaw, available_laws, register_law
from .mass import MassAnalysis, MassResult
from .model import ConicalWindModel, ModelCurve
from .pv import PVDiagram
from .ringfit import RingFitter
from .trace import Trace, TraceExtractor

__all__ = ["__version__", "WindAnalysis", "WindConfig", "WindCube", "PVDiagram", "Trace",
           "TraceExtractor", "ConicalWindModel", "ModelCurve", "ViewGeometry",
           "GlobalFitter", "FitResult", "compare_fits", "RingFitter", "MassAnalysis",
           "MassResult", "WindProperties", "VelocityLaw", "register_law", "available_laws",
           "fit_disk_center"]
