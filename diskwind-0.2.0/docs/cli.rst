Command line
============

.. code-block:: text

   diskwind init [file.toml]          annotated configuration template
   diskwind options [section]         every option with its default and meaning
   diskwind convert old_config.py     old-style *_config.py -> TOML
   diskwind center continuum.fits     fit the disk center to a continuum image
   diskwind check config.toml         noise and disk-rotation checks (+ figures)
   diskwind fit config.toml [options] one global fit
   diskwind run config.toml           the standard analysis (checks, fits, mass, figures)
   diskwind mock [--preset hh30]      synthetic cube + configuration
                                      (elias20: inclined, one lobe; hh30: edge-on, both lobes)

``diskwind fit`` options:

.. code-block:: text

   --law broken-power-law             velocity law
   --radial cone|proportional|independent
   --rotation constant|power-law
   --guide power-law                  label (or path) of a saved result guiding the traces
   --start-from power-law             start from a saved result
   --iterations 2                     guided re-extraction cycles
   --nboot 100                        bootstrap resamplings
   --fix xslope=0 xoff_au=0           fixed parameters
   --label mylabel                    name of the outputs
   --no-plots

Example:

.. code-block:: bash

   diskwind fit elias20.toml --law power-law --nboot 0
   diskwind fit elias20.toml --law broken-power-law --guide power-law \
       --start-from power-law --iterations 2 --nboot 100
