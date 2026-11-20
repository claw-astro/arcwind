# Examples

Run the examples from this folder. 01 and 07 need no data; 03-06 reuse the
synthetic cube written by 01; 02, 08-11 need the Elias 2-20 data.

**Where the data go.** Each source has an input folder `data/<source>/` and an output
folder `output/<source>/`; the configuration in `configs/` points to both (paths are
relative to the configuration file). For Elias 2-20: put `Elias20_CO.fits`,
`Elias20_continuum.fits` and (optionally) `Elias20_CO_M0.fits` in `data/elias20/`;
results go to `output/elias20/`. FITS files and `output/` are ignored by git. See
[`data/README.md`](data/README.md), which also explains how to set up another source.

| script | what it shows |
|---|---|
| `01_quickstart_mock.py` | a synthetic wind: fit it and compare with the true values (writes `mock.toml`) |
| `02_elias20_analysis.py` | the adopted Elias 2-20 analysis, step by step |
| `03_fitting_options.py` | every model, statistic, trace and error option, switched on and off |
| `04_custom_velocity_law.py` | define and fit your own velocity law |
| `05_ring_fits_full_rings.py` | independent ring fits to each cut (full rings) |
| `06_low_level_api.py` | the building blocks, without `WindAnalysis` |
| `07_hh30_like_both_lobes.py` | an HH 30-like bipolar wind, nearly edge-on: both lobes fitted |
| `08_model_definitions.py` | schematic of the model definitions (four panels) |
| `09_model_definitions_meridional.py` | the same with the meridional plane (five panels) |
| `10_ring_schematic.py` | 3-D schematic of the fitted rings and velocities |
| `11_elias20_cut_range_test.py` | Elias 2-20 fit with all cuts vs cuts up to y = 2.0" only, compared |

`configs/` holds complete configurations: `elias20.toml` (the adopted Elias 2-20
settings) and `mock.toml`. `_schematic_helpers.py` is used by 08-10.
