"""
02 -- The Elias 2-20 analysis, step by step (the adopted result).

    python 02_elias20_analysis.py            # needs the cube in data/elias20/ (data/README.md)

Only the red-shifted (far) half of each ring is detected, so:
  * trace.mode = "single" (one peak per position),
  * the flow is tied to the cone (fit.radial = "cone"),
  * the traces follow the INNER arch (mask.x_cone, half-width 0.10 + 0.50 y); a wider
    window also takes in a fainter outer arc and mixes the two shells,
  * the traces are extracted guided by a fit (peak_choice = "nearest"),
  * the final velocity law is a smoothly broken power law.
  * j(z) = j (z/100 au)^qj in the adopted global fit (fit.rotation = "power-law"),
    and a constant j fitted to the same traces for comparison (CONSTANT_J); the
    rotation itself is measured cut by cut (wa.rotation_per_cut).
All of this is set in configs/elias20.toml (cuts y0 = 0.45-2.40", every 0.15").

Every fit starts from scratch: an unguided single power law, then the broken power
law guided by it. YMAX restricts the cuts (example 11 compares several YMAX).
"""
import matplotlib.pyplot as plt

from arcwind import WindAnalysis

YMAX = None
CONSTANT_J = True     # also fit a constant j to the same traces and compare           # e.g. 2.0: only cuts with y0 <= 2.0" (None: all cuts)

wa = WindAnalysis("configs/elias20.toml")
ys = [float(y) for y in wa.config.y_offsets if YMAX is None or y <= YMAX + 1e-6]

# --- 1. checks before fitting -------------------------------------------
print(wa.check_noise().summary())          # are the noise channels line-free?
print(wa.check_disk_rotation().summary())  # disk rotation sense (+1 expected)
wa.plot_slit_layout()                      # where the cuts are
wa.plot_pv_grid(show_mask=True)            # PV diagrams with the masks (cyan)

# --- 2. starting point: a single power law, unguided ---------------------
start = wa.fit("power-law", n_boot=0, y_offsets=ys, label="power-law-first")

# --- 3. adopted fit: broken power law, traces guided by the starting fit ---
bpl = wa.fit("broken-power-law", guide=start, start_from=start, iterations=2,
             n_boot=100, y_offsets=ys, label="broken-power-law")
print(bpl.summary())

# --- 3b. rotation: j(z) or a constant j, and measured cut by cut ----------
if CONSTANT_J:                              # same traces, rotation = constant
    bpl_c = wa.fit("broken-power-law", rotation="constant", traces=bpl.traces,
                   start_from=bpl, n_boot=100, label="broken-power-law-constj")
    print(bpl_c.summary())
    wa.compare(bpl_c, bpl)
# The global law only describes the rotation; it is measured in each cut with the
# geometry and velocity law of bpl held fixed.
rot = wa.rotation_per_cut(bpl)              # -> <label>_rotation_per_cut.txt
print(rot.summary())

# --- 4. is the break justified? a single power law on the same traces ----
pl = wa.fit("power-law", traces=bpl.traces, start_from=bpl, n_boot=0,
            label="power-law-same-traces")
wa.compare(pl, bpl)

# --- 5. mass, mass-loss rate, derived quantities --------------------------
mass = wa.mass(bpl)
props = wa.derived(bpl, mass)
print(props.summary())
print(props.compare_with("HH 30"))

# --- 6. figures ----------------------------------------------------------
wa.plot_cuts(bpl, others=[pl])
for y0 in bpl.traces:                      # one figure per cut: broken-power-law_cut_y<y0>
    wa.plot_cut(y0, bpl, others=[pl])
    plt.close("all")
wa.plot_profiles(bpl, others=[pl])
wa.plot_angular_momentum(bpl, others=[pl])  # j(z) and V_phi(z) of the global fits
wa.plot_rotation(rot, result=bpl)           # j and V_phi measured cut by cut
if CONSTANT_J:
    wa.plot_cuts(bpl_c, others=[bpl])       # constant j (red) vs j(z) (dashed)
    wa.plot_rotation(wa.rotation_per_cut(bpl_c), result=bpl_c)
wa.plot_gallery(bpl, n_boot_curves=30)
wa.plot_bootstrap(bpl)

# --- 7. robustness to the highest cuts: 11_elias20_cut_range_test.py fits the
#        cuts up to y0 = 2.0" and up to 2.4" from the same start and compares them.
