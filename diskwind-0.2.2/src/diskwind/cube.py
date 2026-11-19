"""
The line cube, resampled onto the frame aligned with the wind.

:class:`WindCube` reads the FITS cube, resamples every channel onto a regular
``(X, Y)`` grid (``Y`` along the projected wind axis), measures the noise, and
provides transverse PV diagrams (:meth:`WindCube.pv`), the contamination masks
(:meth:`WindCube.excluded_mask`) and the automatic estimates of Vsys, the
disk-rotation sense and the lobe orientation.
"""
from __future__ import annotations

import warnings

import numpy as np
from scipy.ndimage import map_coordinates, uniform_filter

from .config import WindConfig
from .constants import C_KMS
from .geometry import ViewGeometry, keplerian_speed, wind_to_sky
from .io import CelestialWCS, parse_angle, read_fits
from .pv import PVDiagram

__all__ = ["WindCube", "robust_std"]


def robust_std(a):
    """1.4826 x the median absolute deviation of the finite values of ``a``."""
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return np.nan
    med = np.median(a)
    return 1.4826 * np.median(np.abs(a - med))


class WindCube:
    """A line cube resampled onto the wind frame.

    Parameters
    ----------
    config : WindConfig or str
        Configuration (or path of a configuration file).
    verbose : bool, optional
        Print what is measured (default: ``config.output.verbose``).
    lobe : {'main', 'counter'}
        Which lobe the wind frame follows. ``'counter'`` turns the frame by 180 deg
        (``+Y`` along ``wind_pa_deg + 180``, ``X`` reversed): the opposite lobe is
        then analyzed exactly like the main one, with ``s`` and the disk rotation
        sense reversed (:meth:`ViewGeometry.mirrored`).

    Attributes
    ----------
    x, y : ndarray
        Grid positions across and along the wind axis (arcsec).
    v : ndarray
        Channel velocities kept (LSRK, km/s, increasing).
    cube : ndarray, shape (nv, ny, nx)
        Resampled intensities (``config.grid.units``).
    rms : float
        Channel rms of the resampled noise channels.
    beam_x, beam_y : float
        Beam FWHM projected on X and Y (arcsec).
    geometry : ViewGeometry
        Inclination, lobe orientation ``s``, distance, Vsys, disk sense.

    Examples
    --------
    >>> from diskwind import WindCube
    >>> cube = WindCube("elias20.toml")                       # doctest: +SKIP
    [diskwind] cube: (151, 640, 640) ...; channel rms = 1.2 mJy/beam ...
    >>> pv = cube.pv(1.0)                # transverse PV diagram at y0 = 1"   # doctest: +SKIP
    >>> mask = cube.excluded_mask(1.0)   # pixels the trace may not use       # doctest: +SKIP
    >>> cube.geometry.s, cube.geometry.disk_sense                             # doctest: +SKIP
    (-1, 1)
    """

    def __init__(self, config, verbose=None, lobe="main"):
        if not isinstance(config, WindConfig):
            config = WindConfig.load(config)
        config.validate()
        self.config = cfg = config
        if lobe not in ("main", "counter"):
            raise ValueError("lobe must be 'main' or 'counter'")
        self.lobe = lobe
        self.wind_pa = config.source.wind_pa_deg + (180.0 if lobe == "counter" else 0.0)
        self.verbose = cfg.output.verbose if verbose is None else verbose
        hdr, data, beams = read_fits(cfg.fits_path)
        self.header = hdr
        self._load_axes(hdr, data, beams)

        src = cfg.source
        self.inc = np.radians(src.inclination_deg)
        self.sini, self.cosi = np.sin(self.inc), np.cos(self.inc)
        self.dist = float(src.distance_pc)
        self.mstar = float(src.mstar_msun)
        self.beam_x = self._beam_fwhm_along(self.wind_pa - 90.0)
        self.beam_y = self._beam_fwhm_along(self.wind_pa)

        # disk center
        ra0, de0 = parse_angle(src.center_ra, True), parse_angle(src.center_dec, False)
        if src.center_from_continuum:
            from .continuum import fit_disk_center
            if not src.continuum_fits:
                raise ValueError("center_from_continuum needs source.continuum_fits")
            cc = fit_disk_center(cfg.resolve(src.continuum_fits))
            ra0, de0 = cc.ra_deg, cc.dec_deg
            self.continuum_center = cc
            self._log("disk center from the continuum fit: " + cc.summary().splitlines()[0])
            if any(abs(float(o)) > 0 for o in src.center_offset_arcsec):
                warnings.warn("center_from_continuum is on and center_offset_arcsec is not "
                              "zero: the offset is added to the fitted center")
        elif ra0 is None or de0 is None:
            ra0, de0 = self.wcs.crval
            self._log(f"center = reference pixel (CRVAL) + center_offset_arcsec: "
                      f"RA={ra0:.7f} Dec={de0:.7f}")
        px, py = self.wcs.world_to_pixel0(ra0, de0)
        off = self.wcs.cdinv @ (np.array(src.center_offset_arcsec, float) / 3600.0)
        self.cpix = np.array([px + off[0], py + off[1]])
        self.center_radec = (ra0, de0)

        # velocity crop
        order_v = np.argsort(self.full_vel)
        if cfg.grid.v_range_kms is not None:
            vlo, vhi = cfg.grid.v_range_kms
            order_v = order_v[(self.full_vel[order_v] >= vlo) & (self.full_vel[order_v] <= vhi)]
        self.chan_idx = order_v
        self.v = self.full_vel[order_v]

        # grid
        g = cfg.grid
        self.step = g.step_arcsec or self.pixscale
        xh = g.x_half_width_arcsec
        self.x = np.arange(-xh, xh + 0.5 * self.step, self.step)
        y0, y1 = g.y_range_arcsec
        self.y = np.arange(y0, y1 + 0.5 * self.step, self.step)
        XX, YY = np.meshgrid(self.x, self.y)
        self._grid_pix = self.xy_to_pix(XX, YY)
        self._log(f"cube: {self._data.shape} (nv, ny, nx); pixel {self.pixscale:.4f}\"; beam "
                  f"{self.bmaj:.3f}\"x{self.bmin:.3f}\" PA {self.bpa:.1f}; dv {self.dv:.3f} km/s")
        self.cube = self.regrid(self.chan_idx)

        # noise
        self.noise_idx = self.noise_channel_indices()
        self.v_noise = self.full_vel[self.noise_idx]
        self.noise_cube = self.regrid(self.noise_idx)
        if cfg.noise.rms_jy_beam is not None:
            self.rms = cfg.noise.rms_jy_beam * self.ufac
        else:
            self.rms = robust_std(self.noise_cube)
        self._log(f"channel rms = {self.rms:.4g} {self.units} from {len(self.noise_idx)} "
                  f"channels at v = {np.array2string(np.sort(self.v_noise), precision=2)} km/s")
        self.slit_width = cfg.cuts.slit_width_arcsec or 0.5 * self.beam_y

        # systemic velocity, disk sense, lobe orientation
        self.vsys = src.vsys_kms if src.vsys_kms is not None else self.estimate_vsys()
        self.vsys_err = src.vsys_err_kms
        # the configuration describes the main lobe; the counter lobe's frame
        # reverses both signs (estimates are made in this cube's own frame)
        flip = -1 if lobe == "counter" else 1
        self.dsense = (flip * int(src.disk_rotation_sense) if src.disk_rotation_sense is not None
                       else self.estimate_disk_sense())
        if src.wind_approaching is None:
            self.s = self.estimate_lobe_orientation()
        else:
            self.s = flip * (1 if src.wind_approaching else -1)

    # ------------------------------------------------------------------ setup
    def _load_axes(self, hdr, data, beams):
        cfg = self.config
        naxis = int(hdr["NAXIS"])
        ctypes = {k: str(hdr.get(f"CTYPE{k}", "")).upper() for k in range(1, naxis + 1)}
        ax_spec = [k for k, c in ctypes.items()
                   if c.startswith(("FREQ", "VRAD", "VELO", "VOPT", "FELO", "VLSR"))]
        ax_lon = [k for k, c in ctypes.items() if c.startswith("RA")]
        ax_lat = [k for k, c in ctypes.items() if c.startswith("DEC")]
        if not (ax_spec and ax_lon and ax_lat):
            raise ValueError(f"cannot identify RA/DEC/spectral axes in {ctypes}")
        ax_spec, ax_lon, ax_lat = ax_spec[0], ax_lon[0], ax_lat[0]
        sl = []
        for npax in range(naxis):
            fax = naxis - npax
            sl.append(slice(None) if fax in (ax_spec, ax_lon, ax_lat) else 0)
        d = data[tuple(sl)]
        order = sorted([ax_spec, ax_lat, ax_lon], reverse=True)
        perm = [order.index(a) for a in (ax_spec, ax_lat, ax_lon)]
        self._data = np.transpose(d, perm)               # (nv, ny, nx)
        self.wcs = CelestialWCS(hdr, ax_lon, ax_lat)
        self.pixscale = self.wcs.pixscale_arcsec
        nchan = self._data.shape[0]
        cdl = (float(hdr[f"CDELT{ax_spec}"]) if f"CDELT{ax_spec}" in hdr
               else float(hdr[f"CD{ax_spec}_{ax_spec}"]))
        world = float(hdr[f"CRVAL{ax_spec}"]) + (np.arange(nchan) + 1
                                                 - float(hdr[f"CRPIX{ax_spec}"])) * cdl
        ct = ctypes[ax_spec]
        unit = str(hdr.get(f"CUNIT{ax_spec}", "")).strip()
        self.restfreq = float(hdr.get("RESTFRQ", hdr.get("RESTFREQ", cfg.source.restfreq_hz)))
        if ct.startswith("FREQ"):
            vel = C_KMS * (1.0 - world / self.restfreq)
        else:
            if ct.startswith("VOPT"):
                warnings.warn("optical velocity axis; treated as-is")
            vel = world * (1.0 if unit.lower() in ("km/s", "km s-1") else 1e-3)
        self.full_vel = vel
        self.dv = abs(np.median(np.diff(vel)))
        if "BMAJ" in hdr:
            self.bmaj, self.bmin = hdr["BMAJ"] * 3600, hdr["BMIN"] * 3600
            self.bpa = float(hdr.get("BPA", 0.0))
        elif beams is not None:
            self.bmaj = float(np.median(beams["BMAJ"]))
            self.bmin = float(np.median(beams["BMIN"]))
            self.bpa = float(np.median(beams["BPA"]))
        else:
            raise ValueError("no beam information (BMAJ/BMIN) in the header")
        bunit = str(hdr.get("BUNIT", "Jy/beam")).lower().replace(" ", "")
        if "jy" not in bunit:
            warnings.warn(f"BUNIT={bunit}: assuming Jy/beam")
        self.units = cfg.grid.units
        if self.units == "K":
            self.ufac = self.k_per_jy
        elif self.units == "mJy/beam":
            self.ufac = 1e3
        else:
            self.ufac = 1.0

    @property
    def k_per_jy(self):
        """Rayleigh-Jeans brightness temperature (K) of 1 Jy/beam."""
        return 1.222e6 / ((self.restfreq / 1e9) ** 2 * self.bmaj * self.bmin)

    @property
    def geometry(self):
        """The :class:`~diskwind.geometry.ViewGeometry` of this cube."""
        return ViewGeometry(inclination_deg=float(np.degrees(self.inc)), s=int(self.s),
                            distance_pc=self.dist, vsys=float(self.vsys),
                            vsys_err=float(self.vsys_err), disk_sense=int(self.dsense),
                            mstar_msun=self.mstar)

    def _log(self, msg):
        if self.verbose:
            print(f"[diskwind] {msg}")

    def _beam_fwhm_along(self, pa_deg):
        a = np.radians(pa_deg - self.bpa)
        return float(np.sqrt((self.bmaj * np.cos(a)) ** 2 + (self.bmin * np.sin(a)) ** 2))

    # --------------------------------------------------------- coordinates
    def xy_to_sky(self, X, Y):
        """Wind-frame (X, Y) arcsec -> sky (dRA cos dec, dDec) arcsec from the center."""
        return wind_to_sky(X, Y, self.wind_pa)

    def sky_to_pix(self, dE, dN):
        q = np.tensordot(self.wcs.cdinv, np.array([dE, dN]) / 3600.0, axes=1)
        return self.cpix[0] + q[0], self.cpix[1] + q[1]

    def xy_to_pix(self, X, Y):
        return self.sky_to_pix(*self.xy_to_sky(X, Y))

    def zc_au(self, y0):
        """Height (au) of the ring centered on the cut at ``y0``."""
        return y0 / self.sini * self.dist

    def regrid(self, chans, ppx=None, ppy=None):
        """Resample channels ``chans`` (indices into the FITS cube) onto the grid
        (bilinear interpolation). Returns an array (n, ny, nx) in ``units``."""
        if ppx is None:
            ppx, ppy = self._grid_pix
        ny, nx = self._data.shape[1:]
        x0 = int(max(np.floor(np.nanmin(ppx)) - 2, 0))
        x1 = int(min(np.ceil(np.nanmax(ppx)) + 3, nx))
        y0 = int(max(np.floor(np.nanmin(ppy)) - 2, 0))
        y1 = int(min(np.ceil(np.nanmax(ppy)) + 3, ny))
        if x1 <= x0 or y1 <= y0:
            raise ValueError("resampling window falls outside the image: check the center")
        out = np.empty((len(chans),) + ppx.shape, dtype=np.float32)
        coords = np.array([ppy - y0, ppx - x0])
        for k, c in enumerate(chans):
            plane = np.asarray(self._data[c, y0:y1, x0:x1], dtype=np.float64)
            out[k] = map_coordinates(plane, coords, order=1, cval=np.nan, mode="constant")
        return out * self.ufac

    @property
    def blank_channels(self):
        """Indices (into the FITS cube) of channels without data: all zero or all
        NaN (e.g. edge channels blanked by the imaging). They are never used for
        the noise."""
        if not hasattr(self, "_blank"):
            blank = []
            for c in range(self._data.shape[0]):
                plane = np.asarray(self._data[c])
                if not np.any(np.isfinite(plane) & (plane != 0)):
                    blank.append(c)
            self._blank = np.array(blank, int)
            if len(blank):
                self._log(f"blank channels (all zero/NaN), excluded from the noise: v = "
                          f"{np.array2string(self.full_vel[self._blank], precision=2)} km/s")
        return self._blank

    def noise_channel_indices(self):
        """Indices (into the FITS cube) of the noise channels.

        All channels in ``noise.vel_ranges_kms`` if set, otherwise
        ``noise.n_edge_channels`` at each end of the cube. Blank channels (all
        zero or NaN) are skipped: with ``n_edge_channels`` the next channels
        inwards are taken instead.
        """
        vel = np.asarray(self.full_vel)
        blank = set(self.blank_channels.tolist())
        rngs = self.config.noise.vel_ranges_kms
        if rngs:
            ok = np.zeros(vel.size, bool)
            for lo, hi in rngs:
                ok |= (vel >= min(lo, hi)) & (vel <= max(lo, hi))
            idx = np.array([i for i in np.where(ok)[0] if i not in blank], int)
            if idx.size < 2:
                raise ValueError(f"noise.vel_ranges_kms = {rngs} selects {idx.size} usable "
                                 f"channels; the cube covers {vel.min():.2f} to "
                                 f"{vel.max():.2f} km/s")
            return idx
        nn = int(self.config.noise.n_edge_channels)
        o = [i for i in np.argsort(vel) if i not in blank]
        return np.r_[o[:nn], o[-nn:]].astype(int)

    def sky_map(self, half_size, func="max"):
        """Peak (``func='max'``) or integrated (``'sum'``) map on a regular sky grid,
        east to the left. Returns (dE, dN, map)."""
        st = self.pixscale
        e = np.arange(half_size, -half_size - 0.5 * st, -st)
        n = np.arange(-half_size, half_size + 0.5 * st, st)
        EE, NN = np.meshgrid(e, n)
        ppx, ppy = self.sky_to_pix(EE, NN)
        c = self.regrid(self.chan_idx, ppx, ppy)
        m = np.nanmax(c, axis=0) if func == "max" else np.nansum(c, axis=0) * self.dv
        return e, n, m

    # --------------------------------------------------------------- cuts
    def rows(self, y0, width=None):
        """Boolean selection of grid rows inside the cut at ``y0``."""
        w = self.slit_width if width is None else width
        sel = np.abs(self.y - y0) <= 0.5 * w + 1e-9
        if not np.any(sel):
            sel = np.abs(self.y - y0) == np.min(np.abs(self.y - y0))
        return sel

    def cut_config(self, y0):
        """The configuration seen by the cut at ``y0`` (per-cut overrides applied)."""
        return self.config.for_cut(y0)

    def xsmooth_pix(self, y0):
        xs = self.cut_config(y0).trace.xsmooth_arcsec
        return max(1, int(round(xs / self.step))) if xs else 1

    def pv(self, y0, width=None, with_mask=True, smooth=True):
        """Transverse PV diagram at offset ``y0`` along the wind axis.

        Parameters
        ----------
        y0 : float
            Offset of the cut (arcsec).
        width : float, optional
            Width along Y; default: per-cut override, ``cuts.slit_width_arcsec`` or
            half the projected beam.
        with_mask : bool
            Also compute the contamination mask (``PVDiagram.excluded``).
        smooth : bool
            Apply the per-cut smoothing ``trace.xsmooth_arcsec`` / ``vsmooth_chan``
            (the noise channels are smoothed identically).

        Returns
        -------
        PVDiagram
        """
        c = self.cut_config(y0)
        if width is None:
            width = c.cuts.slit_width_arcsec or self.slit_width
        sel = self.rows(y0, width)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            pv = np.nanmean(self.cube[:, sel, :], axis=1)
            npv = np.nanmean(self.noise_cube[:, sel, :], axis=1)
        nx = self.xsmooth_pix(y0) if smooth else 1
        nv = (c.trace.vsmooth_chan or 1) if smooth else 1
        if nx > 1 or nv > 1:
            pv, npv = _boxcar(pv, nv, nx), _boxcar(npv, nv, nx)
        out = PVDiagram(y0=float(y0), x=self.x, v=self.v, data=pv, rms=robust_std(npv),
                        width=float(width), units=self.units, zc_au=float(self.zc_au(y0)),
                        smoothing=(int(nv), int(nx)))
        if with_mask:
            out.excluded = self.excluded_mask(y0, width=width)
        return out

    def excluded_mask(self, y0, width=None):
        """Boolean (nv, nx): True where a PV pixel must not be used by the trace.
        See :class:`diskwind.masks.ContaminationMask`."""
        from .masks import ContaminationMask
        return ContaminationMask(self).excluded(y0, width=width)

    def x_window(self, y0):
        """Allowed X range (arcsec) of the trace at ``y0``."""
        from .masks import x_window
        return x_window(self.cut_config(y0), y0)

    # ------------------------------------------------- automatic estimates
    def estimate_vsys(self):
        """Vsys from the 20%-of-peak line wings of the spectrum of the disk half
        without wind (Y <= 0), which is symmetric about Vsys."""
        r_ap = self.config.source.vsys_aperture_arcsec
        XX, YY = np.meshgrid(self.x, self.y)
        ap = (np.hypot(XX, YY) <= r_ap) & (YY <= 0)
        if ap.sum() < 4:
            ap = np.hypot(XX, YY) <= r_ap
        spec = np.nansum(self.cube[:, ap], axis=1)
        pk = np.nanmax(spec)
        good = np.where(spec > 0.2 * pk)[0]
        v_lo = _cross(self.v, spec, 0.2 * pk, good[0], -1)
        v_hi = _cross(self.v, spec, 0.2 * pk, good[-1], +1)
        vs = 0.5 * (v_lo + v_hi)
        self._log(f"Vsys estimated from the 20%-peak line wings: {vs:.2f} km/s "
                  "<-- check it and set source.vsys_kms")
        return float(vs)

    def estimate_disk_sense(self):
        """+1 if the +X side of the disk is red-shifted (counter-clockwise rotation
        about the wind axis), from the intensity-weighted velocity on each side."""
        sel_y = (self.y <= 0.0) & (self.y >= -0.5)
        thr = 5 * self.rms
        vm = []
        for sgn in (+1, -1):
            sel_x = (sgn * self.x > max(2 * self.beam_x, 0.15)) & (sgn * self.x < 2.0)
            sub = self.cube[:, sel_y][:, :, sel_x]
            w = np.where(sub > thr, sub, 0.0)
            vm.append(np.nansum(w * self.v[:, None, None]) / max(np.nansum(w), 1e-30))
        d = 1 if vm[0] > vm[1] else -1
        self._log(f"disk rotation: <v>(+X)={vm[0]:.2f}, <v>(-X)={vm[1]:.2f} km/s -> "
                  f"disk_rotation_sense = {d:+d}")
        return d

    def estimate_lobe_orientation(self):
        """+1 if the lobe is blue-shifted on average (axis tilted toward us)."""
        cfg = self.config
        ys = cfg.y_offsets
        ylo = np.min(ys) if ys is not None else 0.3
        yhi = np.max(ys) if ys is not None else self.y.max()
        sel_y = (self.y >= ylo) & (self.y <= yhi)
        sel_x = np.abs(self.x) <= cfg.mask.x_max_arcsec
        sub = self.cube[:, sel_y][:, :, sel_x]
        vrel = (self.v - self.vsys)[:, None, None]
        ok = sub > 5 * self.rms
        for lo, hi in cfg.mask.v_exclude_kms:
            ok &= ~((vrel >= lo) & (vrel <= hi))
        if cfg.source.inclination_deg < 85:
            XX, YY = np.meshgrid(self.x[sel_x], self.y[sel_y])
            vk = self._thin_disk_v(XX, YY)
            ok &= ~(np.abs(vrel - vk[None]) < 1.5)
        w = np.where(ok, sub, 0.0)
        if np.nansum(w) <= 0:
            self._log("no wind emission found to set wind_approaching; using +1")
            return 1
        vmean = np.nansum(w * vrel) / np.nansum(w)
        s = 1 if vmean < 0 else -1
        msg = f"lobe mean velocity (v - Vsys) = {vmean:+.2f} km/s -> wind_approaching = {s > 0}"
        if abs(vmean) < 0.3:
            msg += "  [AMBIGUOUS: set source.wind_approaching explicitly]"
        self._log(msg)
        return s

    def _thin_disk_v(self, X, Y):
        x = X * self.dist
        y = Y * self.dist / max(self.cosi, 1e-3)
        r = np.hypot(x, y)
        v = self.dsense * keplerian_speed(r, self.mstar) * self.sini * x / np.maximum(r, 1e-3)
        return np.where(r <= self.config.mask.disk_rout_au, v, np.nan)

    def disk_height(self, r_au):
        """Height (au) of the line-emitting surfaces at radius ``r_au``."""
        f = self.config.mask.disk_height_func
        return f(r_au) if f is not None else self.config.mask.disk_z_over_r * r_au

    # ------------------------------------------------------- diagnostics
    def check_noise(self, n_edge=None):
        """Are the noise channels line-free? See :func:`diskwind.diagnostics.check_noise`."""
        from .diagnostics import check_noise
        return check_noise(self, n_edge=n_edge)

    def check_disk_rotation(self, y_offsets=(0.0, -0.15, -0.3)):
        """Which disk-rotation sense fits the data? See
        :func:`diskwind.diagnostics.check_disk_rotation`."""
        from .diagnostics import check_disk_rotation
        return check_disk_rotation(self, y_offsets)

    def __repr__(self):
        return (f"WindCube({self.config.source.name!r}, {self.lobe} lobe: {len(self.v)} channels "
                f"{self.v[0]:.2f}..{self.v[-1]:.2f} km/s, grid {len(self.y)}x{len(self.x)}, "
                f"rms {self.rms:.3g} {self.units}, Vsys {self.vsys:.2f}, s {self.s:+d}, "
                f"disk sense {self.dsense:+d})")


def _boxcar(a, nv, nx):
    bad = ~np.isfinite(a)
    out = uniform_filter(np.where(bad, 0.0, a), size=(int(nv), int(nx)), mode="nearest")
    wgt = uniform_filter((~bad).astype(float), size=(int(nv), int(nx)), mode="nearest")
    out = np.where(wgt > 0, out / np.maximum(wgt, 1e-9), np.nan)
    out[bad] = np.nan
    return out


def _cross(v, spec, lev, i, direction):
    j = i + direction
    if j < 0 or j >= len(v):
        return v[i]
    f = (spec[i] - lev) / (spec[i] - spec[j]) if spec[i] != spec[j] else 0.0
    return v[i] + f * (v[j] - v[i])
