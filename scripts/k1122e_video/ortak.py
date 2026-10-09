"""K1122E videolari icin ortak yardimcilar.

Icerik: stil ve renkler, sahne zaman cizelgesi (her kare bastan cizilir), MP4 ve GIF yazimi
(ffmpeg, imageio-ffmpeg paketinden), altyazi kutulari ve basit bir 3B molekul cizici
(ortografik izdusum, derinlige gore siralanmis tup ve kure katmanlari).
"""
from __future__ import annotations

import subprocess
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib

matplotlib.use("Agg")
import imageio_ffmpeg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FFMpegWriter  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402

FFMPEG = imageio_ffmpeg.get_ffmpeg_exe()
GEN, YUK, DPI, FPS = 12.8, 7.2, 100, 30
ARKA, KUTU, YAZI, SOLUK, IZGARA = "#0e1117", "#1b2130", "#ececec", "#9aa3ad", "#2b3240"
RENK = {
    "WT": "#4cc9f0", "K1122R": "#a3d977", "K1122E": "#ff5d73",
    "oligo": "#f2a541", "opc": "#ffd166", "noron": "#6fa8dc", "glia": "#8d99ae",
    "rab10": "#5c7aa6", "switch": "#8ecae6", "heliks1": "#f2a541", "heliks2": "#e76f51",
    "donus": "#e9ecef", "sarmal": "#f2a541", "sarmal_degil": "#e9ecef",
    "arti": "#4c8dff", "eksi": "#ff5d73", "na": "#b388ff", "gtp": "#c5e384", "mg": "#7bd389",
}
plt.rcParams.update({
    "font.family": "DejaVu Sans", "text.color": YAZI, "axes.labelcolor": YAZI,
    "xtick.color": SOLUK, "ytick.color": SOLUK, "axes.edgecolor": IZGARA, "axes.facecolor": ARKA,
    "figure.facecolor": ARKA, "savefig.facecolor": ARKA, "axes.titlecolor": YAZI,
    "animation.ffmpeg_path": FFMPEG,
})


def yumusat(x: float) -> float:
    x = min(max(x, 0.0), 1.0)
    return x * x * (3 - 2 * x)


def aralik(p: float, a: float, b: float) -> float:
    """p, [a, b] araliginda 0'dan 1'e yumusak gecis."""
    return yumusat((p - a) / (b - a)) if b > a else float(p >= a)


@dataclass
class Sahne:
    saniye: float
    ciz: Callable  # ciz(fig, p), p in [0, 1]


def video_yaz(sahneler: list[Sahne], cikti: Path, fps: int = FPS) -> int:
    cikti = Path(cikti)
    cikti.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(GEN, YUK), dpi=DPI)
    yazici = FFMpegWriter(fps=fps, codec="libx264",
                          extra_args=["-pix_fmt", "yuv420p", "-crf", "21", "-preset", "slow", "-movflags", "+faststart"])
    n_kare = 0
    with yazici.saving(fig, str(cikti), DPI):
        for s in sahneler:
            n = max(int(round(s.saniye * fps)), 1)
            for i in range(n):
                fig.clf()
                s.ciz(fig, i / max(n - 1, 1))
                yazici.grab_frame()
                n_kare += 1
    plt.close(fig)
    return n_kare


def gif_yap(mp4: Path, gif: Path, bas: float = 0.0, sure: float | None = None, fps: int = 10, genislik: int = 640,
            renk: int = 256):
    """MP4'ten iki gecisli paletle GIF (README'de satir icinde oynasin diye)."""
    vf = f"fps={fps},scale={genislik}:-1:flags=lanczos"
    zaman = ["-ss", f"{bas}"] + (["-t", f"{sure}"] if sure else [])
    palet = Path(gif).with_suffix(".palet.png")
    subprocess.run([FFMPEG, "-y", *zaman, "-i", str(mp4), "-vf", vf + f",palettegen=max_colors={renk}:stats_mode=diff", str(palet)],
                   check=True, capture_output=True)
    subprocess.run([FFMPEG, "-y", *zaman, "-i", str(mp4), "-i", str(palet), "-lavfi",
                    vf + " [x]; [x][1:v] paletteuse=dither=bayer:bayer_scale=4", str(gif)],
                   check=True, capture_output=True)
    palet.unlink()


def baslik(fig, metin: str, alt: str | None = None, gorunurluk: float = 1.0, y: float = 0.94):
    fig.text(0.04, y, metin, fontsize=25, fontweight="bold", color=YAZI, va="top", alpha=gorunurluk)
    if alt:
        fig.text(0.04, y - 0.07, alt, fontsize=15, color=SOLUK, va="top", alpha=gorunurluk)


def altyazi(fig, metin: str, gorunurluk: float = 1.0, y: float = 0.035, genislik: int = 105, boyut: float = 16):
    if gorunurluk <= 0.01 or not metin:
        return
    satirlar = "\n".join(textwrap.fill(p, genislik) for p in metin.split("\n"))
    fig.text(0.5, y, satirlar, ha="center", va="bottom", fontsize=boyut, color=YAZI, alpha=gorunurluk,
             linespacing=1.35,
             bbox=dict(boxstyle="round,pad=0.55", facecolor=KUTU, edgecolor="#3a4356", alpha=0.94 * gorunurluk))


def altyazi_akisi(fig, p: float, akis: list[tuple[float, float, str]], **kw):
    """akis: (baslangic, bitis, metin) listesi; her metin kendi araliginda yumusak girer ve cikar."""
    for a, b, metin in akis:
        if a <= p <= b:
            g = min(aralik(p, a, a + 0.04), 1 - aralik(p, b - 0.04, b))
            altyazi(fig, metin, g, **kw)


def kart(fig, baslik_metni: str, maddeler: list[str], p: float, alt: str | None = None):
    """Ozet karti: maddeler sirayla belirir."""
    baslik(fig, baslik_metni, alt, gorunurluk=aralik(p, 0, 0.15))
    for i, m in enumerate(maddeler):
        g = aralik(p, 0.1 + i * 0.12, 0.2 + i * 0.12)
        fig.text(0.06, 0.70 - i * 0.12, "•  " + textwrap.fill(m, 95, subsequent_indent="    "),
                 fontsize=17, color=YAZI, va="top", alpha=g, linespacing=1.3)


# ---------------------------------------------------------------- 3B cizici

def pdb_oku(metin: str) -> list[dict]:
    atomlar = []
    for satir in metin.splitlines():
        if satir.startswith(("ATOM", "HETATM")):
            atomlar.append({
                "ad": satir[12:16].strip(), "kalinti": satir[17:20].strip(), "zincir": satir[21],
                "no": int(satir[22:26]), "xyz": np.array([float(satir[30:38]), float(satir[38:46]), float(satir[46:54])]),
                "element": (satir[76:78].strip() or satir[12:16].strip()[0]),
            })
    return atomlar


def donme(y_aci: float, x_aci: float = 0.0) -> np.ndarray:
    cy, sy, cx, sx = np.cos(y_aci), np.sin(y_aci), np.cos(x_aci), np.sin(x_aci)
    ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    return rx @ ry


def kamera_tabani(P: np.ndarray, yukari_nokta: np.ndarray | None = None) -> np.ndarray:
    """Temel bilesenlerden bakis tabani: 1. eksen dikey, 2. eksen yatay, 3. eksen bakis yonu."""
    X = P - P.mean(0)
    _, _, vt = np.linalg.svd(X, full_matrices=False)
    e1, e2, e3 = vt
    if yukari_nokta is not None and (yukari_nokta - P.mean(0)) @ e1 < 0:
        e1 = -e1
    e3 = np.cross(e2, e1)
    return np.array([e2, e1, e3])  # satirlar: ekran x, ekran y, derinlik


def catmull_rom(P: np.ndarray, n: int = 6) -> tuple[np.ndarray, np.ndarray]:
    """Noktalardan gecen yumusak egri; ikinci donus her ara noktanin ait oldugu kalinti indeksi."""
    if len(P) < 2:
        return P, np.zeros(len(P), int)
    Q = np.vstack([P[0] - (P[1] - P[0]), P, P[-1] + (P[-1] - P[-2])])
    egri, sahip = [], []
    for i in range(1, len(Q) - 2):
        p0, p1, p2, p3 = Q[i - 1], Q[i], Q[i + 1], Q[i + 2]
        for t in np.linspace(0, 1, n, endpoint=False):
            t2, t3 = t * t, t * t * t
            egri.append(0.5 * ((2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2
                               + (-p0 + 3 * p1 - 3 * p2 + p3) * t3))
            sahip.append(i - 1 if t < 0.5 else i)
    egri.append(P[-1])
    sahip.append(len(P) - 1)
    return np.array(egri), np.array(sahip)


class Sahne3B:
    """Derinlige gore katmanlanan tup ve kure ilkelleri. ekle_* ile doldur, ciz ile eksene bas."""

    def __init__(self, taban: np.ndarray, merkez: np.ndarray, olcek: float):
        self.T, self.m, self.s = taban, merkez, olcek
        self.cizgiler, self.kureler = [], []  # (x0, y0, x1, y1, z, renk, kalinlik), (x, y, z, renk, boyut, kenar)

    def izdusum(self, P: np.ndarray) -> np.ndarray:
        Q = (np.atleast_2d(P) - self.m) @ self.T.T
        Q[:, :2] *= self.s
        return Q

    def ekle_tup(self, P: np.ndarray, renkler: list, kalinliklar: list, n: int = 6, alfa: float = 1.0):
        egri, sahip = catmull_rom(np.asarray(P, float), n)
        Q = self.izdusum(egri)
        for k in range(len(Q) - 1):
            r = sahip[k]
            self.cizgiler.append((Q[k, 0], Q[k, 1], Q[k + 1, 0], Q[k + 1, 1], 0.5 * (Q[k, 2] + Q[k + 1, 2]),
                                  renkler[r], kalinliklar[r], alfa))

    def ekle_kure(self, P: np.ndarray, renk, boyut: float, kenar: str = "#0b0d12", alfa: float = 1.0):
        for q in self.izdusum(P):
            self.kureler.append((q[0], q[1], q[2], renk, boyut, kenar, alfa))

    def ekle_bag(self, a: np.ndarray, b: np.ndarray, renk, kalinlik: float = 2.5, alfa: float = 1.0):
        qa, qb = self.izdusum(a)[0], self.izdusum(b)[0]
        self.cizgiler.append((qa[0], qa[1], qb[0], qb[1], 0.5 * (qa[2] + qb[2]), renk, kalinlik, alfa))

    def ciz(self, ax, katman: int = 14, sis: float = 0.55):
        z = [c[4] for c in self.cizgiler] + [k[2] for k in self.kureler]
        if not z:
            return
        zmin, zmax = min(z), max(z)
        sinirlar = np.linspace(zmin, zmax + 1e-6, katman + 1)

        def sisle(renk, zz):
            f = sis * (zz - zmin) / max(zmax - zmin, 1e-6)  # uzak (buyuk z) koyulasir
            return tuple((1 - f) * np.array(to_rgb(renk)) + f * np.array(to_rgb(ARKA)))

        for i in range(katman):
            a, b = sinirlar[katman - i - 1], sinirlar[katman - i]  # uzaktan yakina
            cs = [c for c in self.cizgiler if a <= c[4] < b]
            if cs:
                seg = [[(c[0], c[1]), (c[2], c[3])] for c in cs]
                ax.add_collection(LineCollection(seg, colors="#07090d", linewidths=[c[6] + 2.2 for c in cs],
                                                 capstyle="round", zorder=2 * i, alpha=0.9))
                ax.add_collection(LineCollection(seg, colors=[(*sisle(c[5], c[4]), c[7]) for c in cs],
                                                 linewidths=[c[6] for c in cs], capstyle="round", zorder=2 * i + 0.5))
            ks = [k for k in self.kureler if a <= k[2] < b]
            if ks:
                ax.scatter([k[0] for k in ks], [k[1] for k in ks], s=[k[4] for k in ks],
                           c=[(*sisle(k[3], k[2]), k[6]) for k in ks], edgecolors=[k[5] for k in ks],
                           linewidths=0.8, zorder=2 * i + 1)


def bos_eksen(ax, xlim=(-1, 1), ylim=(-1, 1)):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.axis("off")
