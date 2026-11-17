"""Plots of `hlahaplotools stats` results, written as standalone SVG (no plotting library needed)."""

from __future__ import annotations

from html import escape

# Sequential blue scale, light (0) to dark (1).
_SCALE = [(0xF7, 0xFB, 0xFF), (0xC6, 0xDB, 0xEF), (0x6B, 0xAE, 0xD6), (0x21, 0x71, 0xB5), (0x08, 0x30, 0x6B)]
_CELL, _GAP = 52, 2


def _colour(v: float) -> str:
    v = min(max(v, 0.0), 1.0) * (len(_SCALE) - 1)
    i = min(int(v), len(_SCALE) - 2)
    t = v - i
    r, g, b = (round(a + (c - a) * t) for a, c in zip(_SCALE[i], _SCALE[i + 1], strict=True))
    return f"#{r:02x}{g:02x}{b:02x}"


def _fmt(x) -> str:
    return "–" if x is None else f"{x:.3g}"


def ld_svg(ld: list[dict], loci: list[str] | None = None, alpha: float = 0.05) -> str:
    """Two lower-triangle heatmaps of linkage disequilibrium per pair of loci: D′ (Hedrick)
    and Wn (Cohen), on a 0–1 scale, with the value in each cell.

    `ld` is the list returned by `linkage_disequilibrium` (q_value, if present, comes from
    `adjust_fdr`); `loci` gives the order (default: order of first appearance). Pairs with
    q < `alpha` (or P, when there is no q) get an asterisk. Pairs not tested stay blank.
    Hovering a cell shows n, the allele counts, P and q.
    """
    if loci is None:
        loci = []
        for r in ld:
            for l in (r["locus_a"], r["locus_b"]):
                if l not in loci:
                    loci.append(l)
    pos = {l: i for i, l in enumerate(loci)}
    by_pair = {}
    for r in ld:
        a, b = sorted((r["locus_a"], r["locus_b"]), key=pos.__getitem__)
        by_pair[(a, b)] = r
    k = max(len(loci) - 1, 1)
    label_w = max(36, 8 * max((len(l) for l in loci), default=1) + 12)
    panel = label_w + k * _CELL
    top, bottom = 48, 70
    width = max(2 * panel + 40, label_w + 560)  # room for the colour bar and its note
    height = top + k * _CELL + label_w + bottom
    tested = any(r.get("p_value") is not None for r in ld)

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" font-family="Helvetica, Arial, sans-serif" font-size="12">',
        f'<rect width="{width}" height="{height}" fill="#ffffff"/>',
    ]
    for p, (key, name) in enumerate([("d_prime", "D′ (Hedrick)"), ("wn", "Wn (Cohen)")]):
        x0 = p * (panel + 40)
        out.append(f'<text x="{x0 + label_w}" y="24" font-size="14" font-weight="bold">{name}</text>')
        for i, row in enumerate(loci[1:], start=1):
            y = top + (i - 1) * _CELL
            out.append(f'<text x="{x0 + label_w - 6}" y="{y + _CELL / 2 + 4}" text-anchor="end">{escape(row)}</text>')
            for j, col in enumerate(loci[:i]):
                r = by_pair.get((col, row))
                if r is None:
                    continue
                v = r[key]
                x = x0 + label_w + j * _CELL
                q = r.get("q_value", r.get("p_value"))
                star = "*" if q is not None and q < alpha else ""
                ink = "#ffffff" if v > 0.55 else "#1a1a1a"
                tip = (
                    f"{col}~{row}: {name} = {v:.3f}; n = {r['n']}; alleles {r['alleles_a']} × {r['alleles_b']}; "
                    f"P = {_fmt(r.get('p_value'))}; q = {_fmt(r.get('q_value'))}"
                )
                out.append(
                    f'<g><title>{escape(tip)}</title><rect x="{x}" y="{y}" width="{_CELL - _GAP}" '
                    f'height="{_CELL - _GAP}" fill="{_colour(v)}"/><text x="{x + (_CELL - _GAP) / 2}" '
                    f'y="{y + _CELL / 2 + 3}" text-anchor="middle" fill="{ink}">{v:.2f}{star}</text></g>'
                )
        ybase = top + k * _CELL + 6
        for j, col in enumerate(loci[:-1]):
            x = x0 + label_w + j * _CELL + _CELL / 2
            out.append(
                f'<text x="{x}" y="{ybase}" text-anchor="end" transform="rotate(-45 {x} {ybase})">{escape(col)}</text>'
            )
    # Colour bar.
    ybar = height - bottom + 18
    out.append('<defs><linearGradient id="ldscale">')
    for s in range(len(_SCALE)):
        f = s / (len(_SCALE) - 1)
        out.append(f'<stop offset="{f:.2f}" stop-color="{_colour(f)}"/>')
    out.append("</linearGradient></defs>")
    out.append(f'<rect x="{label_w}" y="{ybar}" width="200" height="10" fill="url(#ldscale)" stroke="#999"/>')
    for f in (0, 0.5, 1):
        out.append(f'<text x="{label_w + 200 * f}" y="{ybar + 24}" text-anchor="middle">{f:g}</text>')
    note = f"* q &lt; {alpha:g} (permutation test, Benjamini–Hochberg)" if tested else "No permutation test run"
    out.append(f'<text x="{label_w + 220}" y="{ybar + 9}" fill="#555">{note}</text>')
    out.append("</svg>")
    return "\n".join(out) + "\n"
