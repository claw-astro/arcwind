# diskwind

Geometry and kinematics of protoplanetary disk winds from line cubes, by fitting
a rotating, expanding **conical wind** to the position–velocity (PV) diagrams of
many cuts perpendicular to the wind axis.

diskwind generalises the method of [Louvet et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...618A.120L)
for the edge-on HH 30 outflow to:

- **any inclination**: the exact intersection of the cone with every cut is computed;
- **partially detected rings** (e.g. only the red-shifted half): one model fitted to all cuts at once;
- **flexible models**: constant / power-law / broken power-law velocity laws (or your own),
  flow along the cone or free radial velocity, constant or height-dependent angular momentum.

## Install

```bash
pip install -e ".[all]"        # from this folder; extras: astropy, pyyaml, cmasher
pytest                         # ~20 s on a synthetic cube
```

## Use

```bash
diskwind init mysource.toml    # fill in [source]
diskwind check mysource.toml   # noise channels, disk rotation sense
diskwind run mysource.toml     # fits, mass, derived quantities, figures
```

```python
from diskwind import WindAnalysis

wa = WindAnalysis("mysource.toml")
first = wa.fit("power-law", n_boot=0)
final = wa.fit("broken-power-law", guide=first, start_from=first, iterations=2)
wa.compare(first, final)
print(final.summary())
mass = wa.mass(final)
print(wa.derived(final, mass).summary())
wa.plot_cuts(final, others=[first])
```

Every option (`diskwind options` lists them all) can be set in the configuration
file, on the configuration object, or for one call:

```python
wa.fit(radial="proportional", fixed={"xslope": 0.0}, loss="soft_l1", n_boot=0)
```

## Contents

| | |
|---|---|
| `src/diskwind/` | the package |
| `examples/` | runnable examples (`01_quickstart_mock.py` needs no data) and configurations |
| `examples/paper_elias20/` | schematic figures of the Elias 2-20 paper |
| `docs/` | Sphinx documentation (Read the Docs: `.readthedocs.yaml`) |
| `tests/` | pytest suite on a synthetic cube |

Build the documentation locally with `pip install -e ".[docs]"` and
`sphinx-build docs docs/_build/html`.

## Citing

Please cite the Elias 2-20 paper (in preparation) and Louvet et al. (2018, A&A 618, A120).
