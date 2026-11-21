Installation
============

arcwind needs Python 3.10 or later, with numpy, scipy and matplotlib.

.. code-block:: bash

   pip install arcwind                 # from PyPI (once released)
   pip install "arcwind[all]"          # + astropy (FITS), PyYAML, cmasher

From the source folder:

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
