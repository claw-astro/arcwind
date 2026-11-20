"""Generate the arcwind logo set as self-contained SVGs (text outlined).

Usage: python make_logos.py OUTDIR   (needs fontTools and the Inter fonts)
"""
import math
import os
import sys

from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
from fontTools.ttLib import TTFont

OUT = sys.argv[1] if len(sys.argv) > 1 else "out"
os.makedirs(OUT, exist_ok=True)

NAVY, CREAM, ORANGE, BLUE = "#0F1B2D", "#F7F5F0", "#D9622B", "#2F6DB5"
GRAY_L, GRAY_D = "#5B6472", "#AEB7C4"
FONTDIR = "/usr/share/fonts/opentype/inter/"
BOLD = TTFont(FONTDIR + "InterDisplay-Bold.otf")
MED = TTFont(FONTDIR + "Inter-Medium.otf")


def text_path(font, s, size, x0=0.0, y0=0.0, track=0.0):
    """Outline `s`; baseline at y0, left edge at x0. Returns (d, width)."""
    gs = font.getGlyphSet()
    cmap = font.getBestCmap()
    upm = font["head"].unitsPerEm
    k = size / upm
    pen = SVGPathPen(gs)
    x = 0.0
    for ch in s:
        g = cmap[ord(ch)]
        tp = TransformPen(pen, (k, 0, 0, -k, x0 + x, y0))
        gs[g].draw(tp)
        x += gs[g].width * k + track * size
    x -= track * size
    return pen.getCommands(), x


def cap_height(font, size):
    return font["OS/2"].sCapHeight * size / font["head"].unitsPerEm


def xheight(font, size):
    return font["OS/2"].sxHeight * size / font["head"].unitsPerEm


# ---------------------------------------------------------------- marks
def arch_mark(cx, cy, r, ring=14, dash=3, ndots=9, dot=5.5, base=True,
              ink=NAVY, base_op=0.25, simple=False):
    """Red-shifted half ring of a PV cut, with trace points and blue lower half."""
    s = ring / 14
    out = [f'<path d="M{cx-r:.2f} {cy:.2f} A{r} {r} 0 0 1 {cx+r:.2f} {cy:.2f}" '
           f'stroke="{ORANGE}" stroke-width="{ring}" stroke-linecap="round" fill="none"/>']
    if not simple:
        out.append(f'<path d="M{cx-r:.2f} {cy:.2f} A{r} {r} 0 0 0 {cx+r:.2f} {cy:.2f}" '
                   f'stroke="{BLUE}" stroke-width="{dash}" stroke-dasharray="{3*s:.2f} {8*s:.2f}" '
                   f'stroke-linecap="round" fill="none" opacity="0.75"/>')
    for i in range(ndots):
        a = math.radians(10 + 160 * i / (ndots - 1))
        x = cx + r * math.cos(a)
        y = cy - r * math.sin(a)
        out.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{dot}" fill="{ink}"/>')
    if base:
        out.append(f'<line x1="{cx-r-18*s:.2f}" y1="{cy}" x2="{cx+r+18*s:.2f}" y2="{cy}" '
                   f'stroke="{ink}" stroke-width="{2*s:.2f}" opacity="{base_op}"/>')
    return "\n".join(out)


def helix_mark(ox, oy, k=1.0, disk=CREAM):
    """Two streamlines spiraling up a widening cone (180x170 design box)."""
    red = ("M90 158 C70 150 66 138 90 132 C118 125 122 112 90 104 C52 95 46 80 90 70 "
           "C140 60 146 42 90 30 C40 20 30 12 40 6")
    blu = ("M90 158 C110 150 114 138 90 132 C62 125 58 112 90 104 C128 95 134 80 90 70 "
           "C40 60 34 42 90 30 C140 20 150 12 140 6")
    return (f'<g transform="translate({ox} {oy}) scale({k})">'
            f'<path d="{red}" stroke="{ORANGE}" stroke-width="7" stroke-linecap="round" fill="none"/>'
            f'<path d="{blu}" stroke="{BLUE}" stroke-width="4" stroke-linecap="round" fill="none" opacity="0.85"/>'
            f'<ellipse cx="90" cy="160" rx="60" ry="9" fill="{disk}" opacity="0.9"/></g>')


def svg(w, h, body, bg=None, title="arcwind"):
    rect = f'<rect width="{w}" height="{h}" fill="{bg}"/>' if bg else ""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w:.0f} {h:.0f}" '
            f'width="{w:.0f}" height="{h:.0f}" role="img" aria-label="{title}">'
            f'<title>{title}</title>{rect}{body}</svg>\n')


def write(name, content):
    with open(os.path.join(OUT, name), "w") as f:
        f.write(content)


# ---------------------------------------------------------------- wordmarks
def wordmark(s, size, x0, y0, colors, track=-0.02, font=BOLD):
    """Draw s with per-segment colors: colors = [(text, fill), ...]."""
    parts, x = [], x0
    for seg, fill in colors:
        d, w = text_path(font, seg, size, x, y0, track)
        parts.append(f'<path d="{d}" fill="{fill}"/>')
        x += w + track * size
    return "".join(parts), x - x0 - track * size


def measure(s, size, track=-0.02, font=BOLD):
    return text_path(font, s, size, 0, 0, track)[1]


TAG = "conical disk-wind fitting"
EXPAND = [("A", 1), ("nalysis of ", 0), ("R", 1), ("esolved ", 0), ("C", 1), ("onical ", 0),
          ("W", 1), ("inds ", 0), ("I", 1), ("n ", 0), ("N", 1), ("ascent ", 0), ("D", 1), ("isks", 0)]


def tagline(parts, size, cx, y, base, hi):
    """Centered tagline; parts = [(text, highlight?)]."""
    total = sum(measure(t, size, 0.0, MED) for t, _ in parts)
    x = cx - total / 2
    out = []
    for t, h in parts:
        d, w = text_path(BOLD if h else MED, t, size, x, y, 0.0)
        out.append(f'<path d="{d}" fill="{hi if h else base}"/>')
        x += w
    return "".join(out)


def stacked(name, mark, word_colors, bg, ink_sub, tag_parts, hi=ORANGE,
            W=520, H=400, word=64, title="arcwind logo"):
    wtxt = "".join(t for t, _ in word_colors)
    ww = measure(wtxt, word)
    body = mark
    base_y = 318
    wd, _ = wordmark(wtxt, word, W / 2 - ww / 2, base_y, word_colors)
    body += wd
    if tag_parts:
        body += tagline(tag_parts, 15, W / 2, base_y + 38, ink_sub, hi)
    write(name, svg(W, H, body, bg, title))



def helix_thick(ox, oy, k, disk, red_w=7, blue_w=4):
    red = ("M90 158 C70 150 66 138 90 132 C118 125 122 112 90 104 C52 95 46 80 90 70 "
           "C140 60 146 42 90 30 C40 20 30 12 40 6")
    blu = ("M90 158 C110 150 114 138 90 132 C62 125 58 112 90 104 C128 95 134 80 90 70 "
           "C40 60 34 42 90 30 C140 20 150 12 140 6")
    return (f'<g transform="translate({ox:.2f} {oy:.2f}) scale({k:.4f})">'
            f'<path d="{red}" stroke="{ORANGE}" stroke-width="{red_w}" stroke-linecap="round" fill="none"/>'
            f'<path d="{blu}" stroke="{BLUE}" stroke-width="{blue_w}" stroke-linecap="round" fill="none" opacity="0.9"/>'
            f'<ellipse cx="90" cy="160" rx="60" ry="9" fill="{disk}" opacity="0.9"/></g>')


THEMES = {"light": (CREAM, NAVY), "dark": (NAVY, CREAM)}
for theme, (bg, ink) in THEMES.items():
    # stacked logo: helix over two-tone wordmark, no tagline
    W, H, word = 520, 400, 72
    ww = measure("arcwind", word)
    body = helix_thick(260 - 90, 52, 1.0, ink)
    wd, _ = wordmark("arcwind", word, W / 2 - ww / 2, 312, [("arc", ORANGE), ("wind", ink)])
    for b, nm in ((bg, f"arcwind-logo-{theme}.svg"), (None, f"arcwind-logo-{theme}-transparent.svg")):
        write(nm, svg(W, H, body + wd, b, "arcwind"))

    # horizontal banner: helix left, wordmark right
    word = 124
    ww = measure("arcwind", word)
    k = 1.25
    W, H = 50 + 180 * k + 22 + ww + 60, 260
    body = helix_thick(50, (H - 170 * k) / 2, k, ink)
    xh = cap_height(BOLD, word)
    wd, _ = wordmark("arcwind", word, 50 + 180 * k + 22, H / 2 + xh / 2 - 4,
                     [("arc", ORANGE), ("wind", ink)])
    write(f"arcwind-banner-{theme}.svg", svg(W, H, body + wd, bg, "arcwind"))
    write(f"arcwind-banner-{theme}-transparent.svg", svg(W, H, body + wd, None, "arcwind"))

    # app icon (512) and favicon (thicker strokes for <=48 px)
    S = 512
    k = 2.2
    body = (f'<rect width="{S}" height="{S}" rx="112" fill="{bg}"/>'
            + helix_thick(256 - 90 * k, 256 - 85 * k, k, ink))
    write(f"arcwind-icon-{theme}.svg", svg(S, S, body, None, "arcwind icon"))
    k = 2.45
    body = (f'<rect width="{S}" height="{S}" rx="112" fill="{bg}"/>'
            + helix_thick(256 - 90 * k, 256 - 85 * k, k, ink, red_w=14, blue_w=10))
    write(f"arcwind-favicon-{theme}.svg", svg(S, S, body, None, "arcwind favicon"))

    # bare mark, transparent (Sphinx html_logo)
    write(f"arcwind-mark-{theme}.svg", svg(180, 172, helix_thick(0, 2, 1.0, ink), None, "arcwind"))
print("wrote", sorted(os.listdir(OUT)))
