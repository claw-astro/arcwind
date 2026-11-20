"""
03 -- Turning fitting options on and off.

Every option can be set (a) in the configuration file, (b) on the config
object, or (c) as an argument of WindAnalysis.fit() for one fit only.
Run 01_quickstart_mock.py first (it writes mock.toml and the cube).
"""
from arcwind import WindAnalysis

wa = WindAnalysis("mock.toml", **{"output.dir": "out_03_options", "fit.n_boot": 0})

# Velocity law ------------------------------------------------------------
r_const = wa.fit("constant")                 # Vz constant
r_pl = wa.fit("power-law")                   # Vz ~ (z/100 au)^qv
r_bpl = wa.fit("broken-power-law", start_from=r_pl)   # start from the simpler fit

# Radial velocity -----------------------------------------------------------
r_cone = wa.fit(radial="cone")               # Vr = Vz tan(theta)  (default)
r_free = wa.fit(radial="proportional")       # Vr free (same height dependence)
r_ind = wa.fit("broken-power-law", radial="independent", start_from=r_bpl)  # own law

# Rotation --------------------------------------------------------------------
r_jz = wa.fit(rotation="power-law")          # j(z) = j (z/100 au)^qj

# Fix parameters, set bounds ------------------------------------------------
r_fix = wa.fit(fixed={"xslope": 0.0, "xoff_au": 0.0}, label="no-axis-offset")
r_b = wa.fit(bounds={"tan_theta": (0.3, 0.8)}, label="bounded-theta")

# Statistic -----------------------------------------------------------------
r_l1 = wa.fit(loss="soft_l1", equal_cut_weight=True, label="robust-equal-weights")

# Traces --------------------------------------------------------------------
r_g = wa.fit(guide=r_cone, iterations=2, label="guided")      # traces around a model
tr = wa.extract_traces(snr_min=8)            # stricter traces, for this extraction
r_strict = wa.fit(traces=tr, label="snr8")

# Errors --------------------------------------------------------------------
r_boot = wa.fit(n_boot=30, label="bootstrap30")
print(r_boot.errors("cov"), r_boot.errors("boot"), r_boot.errors("max"), sep="\n")

# Compare everything on the same traces ---------------------------------------
wa.compare(r_const, r_pl, r_bpl, r_free, r_jz, traces=r_cone.traces)

# The same through the configuration object -----------------------------------
cfg = wa.config.copy()
cfg.fit.velocity_law = "broken-power-law"
cfg.fit.fixed = {"xslope": 0.0}
cfg.trace.peak_choice = "nearest"
cfg.save("mock_bpl.toml")
