Installation
============

diskwind needs Python 3.8 or later (3.10+ recommended), with numpy, scipy and matplotlib.

.. code-block:: bash

   pip install diskwind                 # from PyPI (once released)
   pip install "diskwind[all]"          # + astropy (FITS), PyYAML, cmasher

From the source folder:

.. code-block:: bash

   git clone <repository-url> diskwind
   cd diskwind
   pip install -e ".[all,test,docs]"
   pytest                               # ~20 s, on a synthetic cube

Optional packages
-----------------

``astropy``
   Reading FITS files. Without it a built-in reader handles standard
   primary-HDU cubes.
``pyyaml``
   YAML configuration files (TOML works out of the box).
``cmasher``
   Color maps of the slit-layout figure (default ``cmr.freeze``).

Check the installation:

.. code-block:: bash

   diskwind --help
   diskwind options fit        # every option of the [fit] section
