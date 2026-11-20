Quick start
===========

1. Write a configuration
------------------------

.. code-block:: bash

   diskwind init mysource.toml

and fill in the ``[source]`` block -- the cube, distance, stellar mass,
inclination, and the position angles of the disk and of the wind lobe:

.. code-block:: toml

   [source]
   fits_file = "mysource_CO.fits"
   distance_pc = 140.0
   mstar_msun = 0.5
   inclination_deg = 45.0
   disk_pa_deg = 150.0
   wind_pa_deg = 240.0       # disk_pa +/- 90, toward the lobe

   [cuts]
   y_offsets_arcsec = {start = 0.3, stop = 2.0, step = 0.15}

Everything else has a default (see :doc:`configuration`).

2. Check the data
-----------------

.. code-block:: python

   from diskwind import WindAnalysis

   wa = WindAnalysis("mysource.toml")
   print(wa.cube)                           # Vsys, lobe orientation, disk sense
   print(wa.check_noise().summary())        # noise channels line-free?
   print(wa.check_disk_rotation().summary())
   wa.plot_pv_grid(show_mask=True)          # PV diagrams with the masks

3. Fit
------

.. code-block:: python

   res = wa.fit()                    # model as set in [fit]
   print(res.summary())
   res.params["tan_theta"], res.errors()["j"], res.theta_deg

   wa.plot_cuts(res)                 # every cut: trace + model
   wa.plot_profiles(res)             # Vz, Vr, Vp, flow angle, R, Vphi vs height
   wa.plot_angular_momentum(res)     # j(z) and Vphi(z)

4. Mass and derived quantities
------------------------------

.. code-block:: python

   mass = wa.mass(res)
   props = wa.derived(res, mass)
   print(props.summary())            # Vp, r_L, travel time, Mdot, ...

Or all of it from the command line:

.. code-block:: bash

   diskwind check mysource.toml
   diskwind run mysource.toml

Try it without data
-------------------

.. code-block:: bash

   diskwind mock                     # Elias 2-20-like: inclined, one lobe
   diskwind mock --preset hh30       # HH 30-like: nearly edge-on, both lobes
   diskwind run mock.toml --nboot 0
