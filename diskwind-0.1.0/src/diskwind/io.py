"""FITS input and a minimal celestial WCS (identical behaviour with or without astropy)."""
from __future__ import annotations

import warnings

import numpy as np

__all__ = ["read_fits", "write_fits", "CelestialWCS", "parse_angle"]


def _parse_card_value(s):
    s = s.strip()
    if s.startswith("'"):
        end = s.find("'", 1)
        while end != -1 and end + 1 < len(s) and s[end + 1] == "'":
            end = s.find("'", end + 2)
        return s[1:end].replace("''", "'").rstrip()
    s = s.split("/")[0].strip()
    if s == "T":
        return True
    if s == "F":
        return False
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s.replace("D", "E"))
    except ValueError:
        return s


def _read_fits_minimal(path):
    hdr = {}
    with open(path, "rb") as fh:
        nblock = 0
        done = False
        while not done:
            block = fh.read(2880)
            if len(block) < 2880:
                raise IOError("truncated FITS header")
            nblock += 1
            for k in range(36):
                card = block[80 * k:80 * (k + 1)].decode("ascii", "replace")
                key = card[:8].strip()
                if key == "END":
                    done = True
                    break
                if card[8:10] == "= ":
                    hdr[key] = _parse_card_value(card[10:])
    bitpix = hdr["BITPIX"]
    naxis = hdr["NAXIS"]
    shape = tuple(hdr[f"NAXIS{k}"] for k in range(naxis, 0, -1))
    dt = {8: ">u1", 16: ">i2", 32: ">i4", 64: ">i8", -32: ">f4", -64: ">f8"}[bitpix]
    data = np.memmap(path, dtype=dt, mode="r", offset=2880 * nblock, shape=shape)
    bs, bz = hdr.get("BSCALE", 1.0), hdr.get("BZERO", 0.0)
    if bs != 1.0 or bz != 0.0:
        data = data * bs + bz
    return hdr, data, None


def read_fits(path):
    """Read a FITS cube.

    Uses astropy if installed, otherwise a minimal primary-HDU reader.

    Returns
    -------
    header : mapping
    data : ndarray (memory-mapped)
    beams : table or None
        The BEAMS extension of per-channel beams, if any.
    """
    try:
        from astropy.io import fits
    except ImportError:
        return _read_fits_minimal(path)
    hdul = fits.open(path, memmap=True)
    hdu = next(h for h in hdul if h.data is not None and h.data.ndim >= 2)
    beams = None
    for h in hdul:
        if h.name == "BEAMS":
            beams = h.data
    return hdu.header, hdu.data, beams


# --------------------------------------------------------------------------
# Celestial projection (SIN/TAN), done locally so behaviour is identical
# with or without astropy.
# --------------------------------------------------------------------------
def parse_angle(val, is_ra):
    """Angle in degrees from a number or a sexagesimal string ('16:26:18.87')."""
    if val is None or isinstance(val, (int, float, np.floating)):
        return val
    parts = str(val).replace("h", ":").replace("d", ":").replace("m", ":") \
        .replace("s", "").replace(" ", ":").split(":")
    parts = [p for p in parts if p != ""]
    sign = -1.0 if parts[0].strip().startswith("-") else 1.0
    a = abs(float(parts[0])) + float(parts[1]) / 60 + float(parts[2]) / 3600
    return sign * a * (15.0 if is_ra else 1.0)


class CelestialWCS:
    """Minimal celestial WCS (SIN/TAN) of the two sky axes of a cube header.

    Parameters
    ----------
    hdr : mapping
        FITS header.
    ax_lon, ax_lat : int
        FITS axis numbers of RA and Dec.
    """
    def __init__(self, hdr, ax_lon, ax_lat):
        self.ax = (ax_lon, ax_lat)
        g = lambda k, d=None: hdr.get(k, d)
        self.crval = np.array([g(f"CRVAL{ax_lon}"), g(f"CRVAL{ax_lat}")], float)
        self.crpix = np.array([g(f"CRPIX{ax_lon}"), g(f"CRPIX{ax_lat}")], float)
        proj = str(g(f"CTYPE{ax_lon}"))[-3:]
        if proj not in ("SIN", "TAN"):
            warnings.warn(f"projection {proj} treated as SIN (fine for small fields)")
            proj = "SIN"
        self.proj = proj
        a, b = ax_lon, ax_lat
        if f"CD{a}_{a}" in hdr:
            cd = np.array([[g(f"CD{a}_{a}", 0.0), g(f"CD{a}_{b}", 0.0)],
                           [g(f"CD{b}_{a}", 0.0), g(f"CD{b}_{b}", 0.0)]], float)
        else:
            cdelt = np.array([g(f"CDELT{a}"), g(f"CDELT{b}")], float)
            pc = np.array([[g(f"PC{a}_{a}", g(f"PC0{a}_0{a}", 1.0)),
                            g(f"PC{a}_{b}", g(f"PC0{a}_0{b}", 0.0))],
                           [g(f"PC{b}_{a}", g(f"PC0{b}_0{a}", 0.0)),
                            g(f"PC{b}_{b}", g(f"PC0{b}_0{b}", 1.0))]], float)
            crota = g(f"CROTA{b}", 0.0) or 0.0
            if crota and not any(k.startswith("PC") for k in hdr.keys()):
                c, s = np.cos(np.radians(crota)), np.sin(np.radians(crota))
                pc = np.array([[c, -s], [s, c]])
            cd = np.diag(cdelt) @ pc
        self.cd = cd
        self.cdinv = np.linalg.inv(cd)
        self.pixscale_arcsec = np.sqrt(abs(np.linalg.det(cd))) * 3600.0

    def world_to_pixel0(self, ra, dec):
        """(ra, dec) deg -> 0-based pixel (x, y)."""
        a0, d0 = np.radians(self.crval)
        a, d = np.radians(ra), np.radians(dec)
        l_ = np.cos(d) * np.sin(a - a0)
        m_ = np.sin(d) * np.cos(d0) - np.cos(d) * np.sin(d0) * np.cos(a - a0)
        n_ = np.sin(d) * np.sin(d0) + np.cos(d) * np.cos(d0) * np.cos(a - a0)
        if self.proj == "TAN":
            l_, m_ = l_ / n_, m_ / n_
        xy = np.degrees(np.array([l_, m_]))
        p = self.cdinv @ xy + self.crpix - 1.0
        return p[0], p[1]


def write_fits(path, data, hdr):
    """Write ``data`` with the header dict ``hdr`` (astropy if available)."""
    try:
        from astropy.io import fits
        h = fits.Header()
        for k, v in hdr.items():
            if k not in ("SIMPLE", "BITPIX", "NAXIS") and not k.startswith("NAXIS"):
                h[k] = v
        fits.PrimaryHDU(data, header=h).writeto(path, overwrite=True)
        return
    except ImportError:
        pass
    cards = []
    for k, v in hdr.items():
        if isinstance(v, bool):
            vs = "T" if v else "F"
            cards.append(f"{k:<8}= {vs:>20}")
        elif isinstance(v, (int, np.integer)):
            cards.append(f"{k:<8}= {v:>20d}")
        elif isinstance(v, (float, np.floating)):
            cards.append(f"{k:<8}= {v:>20.12E}")
        else:
            cards.append(f"{k:<8}= '{v:<8}'")
    cards.append("END")
    txt = "".join(c.ljust(80) for c in cards)
    txt += " " * ((-len(txt)) % 2880)
    raw = data.astype(">f4").tobytes()
    raw += b"\0" * ((-len(raw)) % 2880)
    with open(path, "wb") as fh:
        fh.write(txt.encode("ascii"))
        fh.write(raw)
