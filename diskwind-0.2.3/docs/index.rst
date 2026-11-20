diskwind
========

**diskwind** measures the geometry and kinematics of protoplanetary disk winds
from line cubes (e.g. ALMA CO), by fitting a rotating, expanding conical wind to
the position--velocity (PV) diagrams of many cuts perpendicular to the wind axis.

It generalizes the method of `Louvet et al. (2018)
<https://ui.adsabs.harvard.edu/abs/2018A%26A...618A.120L>`_ for the edge-on HH 30
outflow to

* **any inclination**: a cut no longer samples a single height, so the exact
  intersection of the cone with every cut is computed;
* **partially detected rings** (e.g. only the red-shifted half): one model is
  fitted to all cuts at once;
* **bipolar winds**: one mirror-symmetric wind fitted to both lobes, or each lobe
  on its own;
* **flexible velocity laws** (constant, power law, broken power law, or your own),
  radial-velocity and rotation laws, all switchable from a configuration file.

.. code-block:: python

   from diskwind import WindAnalysis

   wa = WindAnalysis("mysource.toml")
   wa.check_noise(); wa.check_disk_rotation()
   first = wa.fit("power-law", n_boot=0)
   final = wa.fit("broken-power-law", guide=first, start_from=first)
   print(final.summary())
   wa.plot_cuts(final)

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   installation
   quickstart
   concepts

.. toctree::
   :maxdepth: 2
   :caption: User guide

   workflow
   fitting_options
   configuration
   outputs
   faq

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   tutorials/mock
   tutorials/elias20
   tutorials/hh30_like
   tutorials/custom_law
   tutorials/ring_fits

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api/index
   cli
   citing
   changelog
