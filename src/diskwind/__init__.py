"""Deprecated alias: ``diskwind`` was renamed ``arcwind`` in version 0.3.0.

``import diskwind`` and ``import diskwind.<module>`` still work and give the
``arcwind`` modules themselves (same objects, so results and classes are
interchangeable), with a ``FutureWarning``. The alias will be removed in a later
release; replace ``diskwind`` by ``arcwind`` in your scripts.
"""
import importlib
import importlib.abc
import importlib.util
import sys
import warnings

warnings.warn("'diskwind' has been renamed 'arcwind'; import arcwind instead "
              "(the diskwind alias will be removed in a later release).",
              FutureWarning, stacklevel=2)


class _AliasFinder(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    """Serve ``diskwind.X`` as the already-imported ``arcwind.X``."""

    def find_spec(self, name, path=None, target=None):
        if name.startswith("diskwind."):
            return importlib.util.spec_from_loader(name, self)
        return None

    def create_module(self, spec):
        return importlib.import_module("arcwind" + spec.name[len("diskwind"):])

    def exec_module(self, module):
        pass


if not any(isinstance(f, _AliasFinder) for f in sys.meta_path):
    sys.meta_path.insert(0, _AliasFinder())

import arcwind as _arcwind  # noqa: E402

sys.modules[__name__] = _arcwind
