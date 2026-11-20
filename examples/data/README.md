# Data for the source-specific examples

Each source has its own folder here, and its results go to `../output/<source>/`:

```
examples/
  configs/elias20.toml        paths below are relative to this file
  data/elias20/               inputs (not in git)
      Elias20_CO.fits             12CO 2-1 line cube        (source.fits_file)
      Elias20_continuum.fits      continuum image           (source.continuum_fits)
      Elias20_CO_M0.fits          moment-0 map, optional    (plot.mom0_fits)
  output/elias20/             results of 02 and 11 (created on the first run; not in git)
```

FITS files are ignored by git (`.gitignore`), so the cubes stay on your disk only.

## Another source

1. Make `data/<source>/` and put its cube (and continuum, moment-0 map) there.
2. Copy a configuration, e.g. `cp configs/elias20.toml configs/<source>.toml`, and set
   `fits_file = "../data/<source>/<cube>.fits"` (and `continuum_fits`, `plot.mom0_fits`),
   `output.dir = "../output/<source>"` and the source parameters.
3. `arcwind check configs/<source>.toml`, then run as in `02_elias20_analysis.py`.

Absolute paths also work, e.g. to keep large cubes elsewhere.
