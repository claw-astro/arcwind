"""
02 -- The Elias 2-20 analysis, step by step (the adopted result).

    python 02_elias20_analysis.py            # needs Elias20_CO.fits next to configs/

Only the red-shifted (far) half of each ring is detected, so:
  * trace.mode = "single" (one peak per position),
  * the flow is tied to the cone (fit.radial = "cone"),
  * the traces are extracted guided by a first fit (peak_choice = "nearest"),
  * the final velocity law is a smoothly broken power law.
All of this is set in configs/elias20.toml.
"""
from diskwind import WindAnalysis

wa = WindAnalysis("configs/elias20.toml")

# --- 1. checks before fitting -------------------------------------------
print(wa.check_noise().summary())          # are the noise channels line-free?
print(wa.check_disk_rotation().summary())  # disk rotation sense (+1 expected)
wa.plot_slit_layout()                      # where the cuts are
wa.plot_pv_grid(show_mask=True)            # PV diagrams with the masks (cyan)

# --- 2. first fit: single power law, unguided, no bootstrap (fast) --------
pl = wa.fit("power-law", n_boot=0, label="power-law-first")

# --- 3. adopted fit: broken power law, traces guided by the first fit -----
bpl = wa.fit("broken-power-law", guide=pl, start_from=pl, iterations=2, n_boot=100)
print(bpl.summary())

# --- 4. is the break justified? same traces for both models ---------------
wa.compare(pl, bpl)

# --- 5. mass, mass-loss rate, derived quantities --------------------------
mass = wa.mass(bpl)
props = wa.derived(bpl, mass)
print(props.summary())
print(props.compare_with("HH 30"))

# --- 6. figures ----------------------------------------------------------
wa.plot_cuts(bpl, others=[pl])
wa.plot_cut(2.0, bpl, others=[pl])
wa.plot_profiles(bpl, others=[pl])
wa.plot_gallery(bpl, n_boot_curves=30)
wa.plot_bootstrap(bpl)
