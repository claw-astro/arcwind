# Examples

| script | what | data |
|---|---|---|
| `01_quickstart_mock.py` | synthetic cube, fit, compare with the truth (writes `mock.toml`) | none |
| `02_elias20_analysis.py` | the adopted Elias 2-20 analysis, step by step | `Elias20_CO.fits` |
| `03_fitting_options.py` | every model / statistic / trace / error option, on and off | run 01 first |
| `04_custom_velocity_law.py` | register and fit your own velocity law | run 01 first |
| `05_ring_fits_full_rings.py` | per-cut ring fits (full rings, e.g. edge-on) | run 01 first |
| `06_low_level_api.py` | the building blocks without `WindAnalysis` | run 01 first |

`configs/elias20.toml` and `configs/mock.toml` are complete configurations;
`paper_elias20/` has the schematic figures of the paper (run after the fit:
`python make_model_definitions.py --config ../configs/elias20.toml`).
