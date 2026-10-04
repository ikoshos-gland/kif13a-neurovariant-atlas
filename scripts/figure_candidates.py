"""Draw the KIF13A candidate landscape as a standalone SVG (no plotting library needed).

x: protein position (MANE, 1805 aa), y: AlphaMissense score, colour: evidence tier.
Dashed lines: the KIF1A-calibrated cutoff 0.968 and AlphaMissense's generic 0.564.

    mcp\\atlas_server\\.venv\\Scripts\\python.exe scripts\\figure_candidates.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))

from atlas_mcp import core  # noqa: E402

W, H = 1100, 460
L, R, T, B = 70, 30, 30, 90          # margins
PLOT_W, PLOT_H = W - L - R, H - T - B
LENGTH = 1805
CUTOFF, GENERIC = 0.968, 0.564
DOMAINS = [(5, 352, "motor", "#cfe3f7"), (353, 380, "boyun", "#e6e5df"), (381, 468, "CC1", "#e6e5df"),
           (469, 519, "FHA", "#f7dfc9"), (520, 1805, "kuyruk", "#efeee9")]
COLOURS = {"A": "#e34948", "B": "#eb6834", "C": "#2a78d6", "D": "#6250d6", "-": "#c9c8c0"}
NAMES = {"A": "A: model + paralog", "B": "B: model + korunmus kalinti", "C": "C: yalniz model",
         "D": "D: paralog kanitli, model esigin altinda", "-": "esigin altinda"}


def x(pos: float) -> float:
    return L + PLOT_W * (pos - 1) / (LENGTH - 1)


def y(score: float) -> float:
    return T + PLOT_H * (1 - score)


def main() -> int:
    rows = list(csv.DictReader(open(sorted((ROOT / "data" / "processed").glob("KIF13A_candidates_*.tsv"))[-1],
                                    encoding="utf-8"), delimiter="\t"))
    pts = [r for r in rows if r["am_pathogenicity"]]
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Arial, Helvetica, sans-serif">',
           f'<rect width="{W}" height="{H}" fill="#ffffff"/>']
    # domain bands
    for a, b, label, colour in DOMAINS:
        out.append(f'<rect x="{x(a):.1f}" y="{T}" width="{x(b) - x(a):.1f}" height="{PLOT_H}" fill="{colour}" opacity="0.55"/>')
        out.append(f'<text x="{(x(a) + x(b)) / 2:.1f}" y="{T + PLOT_H + 16}" font-size="11" fill="#52514e" text-anchor="middle">{label}</text>')
    # axes and gridlines
    for value in (0, 0.25, 0.5, 0.75, 1.0):
        out.append(f'<line x1="{L}" y1="{y(value):.1f}" x2="{L + PLOT_W}" y2="{y(value):.1f}" stroke="#e1e0d9"/>')
        out.append(f'<text x="{L - 8}" y="{y(value) + 4:.1f}" font-size="11" fill="#898781" text-anchor="end">{value:.2f}</text>')
    for pos in (1, 250, 500, 750, 1000, 1250, 1500, 1805):
        out.append(f'<text x="{x(pos):.1f}" y="{T + PLOT_H + 32}" font-size="11" fill="#898781" text-anchor="middle">{pos}</text>')
    # cutoffs
    for value, label, colour in ((CUTOFF, f"KIF1A ile kalibre edilen esik {CUTOFF}", "#a32d2d"),
                                 (GENERIC, f"AlphaMissense genel esigi {GENERIC}", "#898781")):
        out.append(f'<line x1="{L}" y1="{y(value):.1f}" x2="{L + PLOT_W}" y2="{y(value):.1f}" stroke="{colour}" stroke-width="1.2" stroke-dasharray="6 4"/>')
        out.append(f'<text x="{L + PLOT_W - 4}" y="{y(value) - 5:.1f}" font-size="11" fill="{colour}" text-anchor="end">{label}</text>')
    # points: background first, tiers on top
    for tier in ("-", "C", "B", "D", "A"):
        for r in (p for p in pts if p["tier"] == tier):
            conf = r["confidently_observed"] == "True"
            radius = 3.4 if tier != "-" else 1.8
            out.append(f'<circle cx="{x(int(r["protein_pos"])):.1f}" cy="{y(float(r["am_pathogenicity"])):.1f}" r="{radius}" '
                       f'fill="{COLOURS[tier] if conf else "#ffffff"}" stroke="{COLOURS[tier]}" stroke-width="{1.1 if tier != "-" else 0.4}" '
                       f'opacity="{0.95 if tier != "-" else 0.5}"/>')
    # label the strongest confident tier A points, alternating above and below to avoid collisions
    top = sorted((r for r in pts if r["tier"] == "A" and r["confidently_observed"] == "True"),
                 key=lambda r: -float(r["am_pathogenicity"]))[:10]
    for i, r in enumerate(sorted(top, key=lambda r: int(r["protein_pos"]))):
        px, py = x(int(r["protein_pos"])), y(float(r["am_pathogenicity"]))
        dy = -10 if i % 2 == 0 else 16
        out.append(f'<text x="{px:.1f}" y="{py + dy:.1f}" font-size="11" fill="#0b0b0b" text-anchor="middle">{r["protein_change"]}</text>')
    # legend
    lx, ly = L, H - 34
    for tier in ("A", "B", "C", "D", "-"):
        n = sum(1 for p in pts if p["tier"] == tier)
        out.append(f'<circle cx="{lx + 5}" cy="{ly - 4}" r="4" fill="{COLOURS[tier]}"/>')
        out.append(f'<text x="{lx + 14}" y="{ly}" font-size="11" fill="#52514e">{NAMES[tier]} ({n})</text>')
        lx += 14 + 7.2 * len(NAMES[tier]) + 34
    out.append(f'<circle cx="{L + 5}" cy="{H - 14}" r="4" fill="#ffffff" stroke="#e34948" stroke-width="1.1"/>')
    out.append(f'<text x="{L + 14}" y="{H - 10}" font-size="11" fill="#52514e">ici bos daire: gnomAD AC0 filtresi, guvenilir tasiyici yok</text>')
    out.append(f'<text x="{L}" y="{T - 12}" font-size="13" fill="#0b0b0b">KIF13A missense varyantlari (gnomAD v4.1, n={len(pts)}): AlphaMissense skoru ve kanit katmani</text>')
    out.append(f'<text x="14" y="{T + PLOT_H / 2:.0f}" font-size="11" fill="#52514e" transform="rotate(-90 14 {T + PLOT_H / 2:.0f})" text-anchor="middle">AlphaMissense skoru</text>')
    out.append(f'<text x="{L + PLOT_W / 2:.0f}" y="{T + PLOT_H + 52}" font-size="11" fill="#52514e" text-anchor="middle">protein pozisyonu (MANE, 1805 aa)</text>')
    out.append("</svg>")

    path = ROOT / "figures" / "KIF13A_candidates_landscape.svg"
    path.write_text("\n".join(out), encoding="utf-8")
    print(core.rel(path), path.stat().st_size, "bytes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
