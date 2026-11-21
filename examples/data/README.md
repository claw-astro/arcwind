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

## Getting the Elias 2-20 data

The Elias 2-20 files are the fiducial images of the DSHARP ALMA Large Program
(Andrews et al. 2018, ApJL 869, L41), public in the
[DSHARP Data Release](https://almascience.org/alma-data/lp/DSHARP)
(mirror: <https://bulk.cv.nrao.edu/almadata/lp/DSHARP/>). From the row for
Elias 20, download the 12CO J=2-1 image cube and the 1.25 mm continuum image, and
save them in `data/elias20/` as `Elias20_CO.fits` and `Elias20_continuum.fits`.
The moment-0 map (`Elias20_CO_M0.fits`) is optional; it is only used as the
background of the slit-layout figure and can be made from the cube (e.g. with
`bettermoments`). The calibrated measurement sets are also in the data release
and in the ALMA archive (project 2016.1.00484.L).

## Another source

1. Make `data/<source>/` and put its cube (and continuum, moment-0 map) there.
2. Copy a configuration, e.g. `cp configs/elias20.toml configs/<source>.toml`, and set
   `fits_file = "../data/<source>/<cube>.fits"` (and `continuum_fits`, `plot.mom0_fits`),
   `output.dir = "../output/<source>"` and the source parameters.
3. `arcwind check configs/<source>.toml`, then run as in `02_elias20_analysis.py`.

Absolute paths also work, e.g. to keep large cubes elsewhere.
