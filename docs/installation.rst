Installation
============

arcwind needs Python 3.10 or later, with numpy, scipy and matplotlib.

.. code-block:: bash

   pip install arcwind                 # from PyPI
   pip install "arcwind[all]"          # + astropy (FITS), PyYAML, cmasher
   pip install --upgrade arcwind       # update to the latest release

The latest development version, straight from GitHub:

.. code-block:: bash

   pip install "arcwind[all] @ git+https://github.com/claw-astro/arcwind.git"

If pip answers ``Could not find a version that satisfies the requirement arcwind
(from versions: none)``, the Python in use is older than 3.10 (check with
``python --version``): pip then hides arcwind. Install into an environment with
Python 3.10 or later, e.g. ``conda create -n arcwind python=3.12`` and
``conda activate arcwind``.

From the source folder (for development):

.. code-block:: bash

   git clone https://github.com/claw-astro/arcwind.git
   cd arcwind
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

   arcwind --help
   arcwind options fit        # every option of the [fit] section
