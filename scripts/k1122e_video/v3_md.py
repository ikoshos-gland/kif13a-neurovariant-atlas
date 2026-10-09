"""Video 3: Ayni baslangictan WT, K1122R ve K1122E, 100 ns kompleks MD (Rab10-GTP-Mg + KIF13A 1085-1160).

Girdi: export_viewer2.py ciktisi (s3 baslangic yapisi, 3 kosu, 0,2 ns arayla 450 kare; hepsi ilk sistemin
Rab10 CA'larina hizali). Uretim komutu (TRUBA, debug kuyrugu, yaklasik 4 dk):
  python export_viewer2.py --run WT_wt_s3=runs/WT_wt_s3 --run K1122R_paired_wt_s3=runs/K1122R_paired_wt_s3 \
      --run K1122E_paired_wt_s3=runs/K1122E_paired_wt_s3 --stride3d 20 --stride-chart 10 --out viewer_s3_frames.json
Cikti: figures/videos/03_md.mp4 ve 03_md.gif
Kullanim: python scripts/k1122e_video/v3_md.py --kareler viewer_s3_frames.json [--hizli]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ortak import (ARKA, IZGARA, RENK, SOLUK, YAZI, Sahne, Sahne3B, altyazi_akisi, aralik, baslik, bos_eksen,
                   donme, gif_yap, kart, pdb_oku, video_yaz)

KOK = Path(__file__).resolve().parents[2]
ap = argparse.ArgumentParser()
ap.add_argument("--kareler", type=Path,
                default=KOK.parent / "so3lr-local" / "truba_k1122e" / "complex_md" / "viewer" / "viewer_s3_frames.json")
ap.add_argument("--hizli", action="store_true")
a = ap.parse_args()
veri = json.loads(a.kareler.read_text())["systems"]
AD = {"WT_wt_s3": ("WT", "K1122", RENK["WT"]), "K1122R_paired_wt_s3": ("K1122R", "R1122", RENK["K1122R"]),
      "K1122E_paired_wt_s3": ("K1122E", "E1122", RENK["K1122E"])}
OMURGA = {"N", "CA", "C", "O"}
ELEMENT = {"C": "#c8ccd2", "N": RENK["arti"], "O": RENK["eksi"]}
YUKLU = {"NZ", "NH1", "NH2", "NE", "OE1", "OE2"}


class Sistem:
    def __init__(self, d):
        self.id = d["system_id"]
        self.ad, self.k1122, self.renk = AD[self.id]
        self.atomlar = pdb_oku(d["pdb"])
        self.X = np.asarray(d["coords"], dtype=np.float32).reshape(d["n_frames"], -1, 3) / 10.0
        self.t = np.asarray(d["times_ns"])
        self.heliks = {int(k): np.asarray(v) for k, v in d["frame_helix"].items()}
        ch = d["chart"]
        self.tc = np.asarray(ch["t_ns"])
        H = np.array([ch["helix"][str(r)] for r in range(1116, 1127)], dtype=float).mean(0)
        k = 20  # 2 ns kayan ortalama (grafik noktalari 0,1 ns arayla)
        self.donus = np.convolve(H, np.ones(k) / k, mode="same")
        self.ca_rab = [i for i, x in enumerate(self.atomlar) if x["zincir"] == "A" and x["ad"] == "CA"]
        self.rab_no = [self.atomlar[i]["no"] for i in self.ca_rab]
        self.ca_kif = [i for i, x in enumerate(self.atomlar) if x["zincir"] == "B" and x["ad"] == "CA"]
        self.kif_no = [self.atomlar[i]["no"] for i in self.ca_kif]
        self.na = [i for i, x in enumerate(self.atomlar) if x["zincir"] == "E"]
        self.yan = {key: [i for i, x in enumerate(self.atomlar) if x["zincir"] == z and x["no"] == n
                          and x["ad"] not in OMURGA] for key, (z, n) in
                    {"1122": ("B", 1122), "1130": ("B", 1130), "70": ("A", 70)}.items()}

    def kare(self, f: float) -> np.ndarray:
        i = int(np.floor(f))
        j = min(i + 1, len(self.X) - 1)
        u = f - i
        return (1 - u) * self.X[i] + u * self.X[j]


S = [Sistem(d) for d in veri]
N = len(S[0].t)

# kamera: WT ilk kare; KIF13A ustte, sac tokasinin uzun ekseni yatay, dönüs ucgeni en genis gorunecek sekilde
X0 = S[0].X[0]
rab0, kif0 = X0[S[0].ca_rab], X0[S[0].ca_kif]
yukari = kif0.mean(0) - rab0.mean(0)
yukari /= np.linalg.norm(yukari)
_, _, vt = np.linalg.svd(kif0 - kif0.mean(0))
uzun = vt[0] - (vt[0] @ yukari) * yukari
uzun /= np.linalg.norm(uzun)
TABAN = np.array([uzun, yukari, np.cross(uzun, yukari)])
a0 = S[0].atomlar


def nokta(sist, X, key, adlar):
    return np.mean([X[i] for i in sist.yan[key] if sist.atomlar[i]["ad"] in adlar], axis=0)


UCLU = np.array([nokta(S[0], X0, "1122", {"NZ"}), nokta(S[0], X0, "1130", {"OE1", "OE2"}),
                 nokta(S[0], X0, "70", {"NH1", "NH2", "NE"})])
en, ACI, EGIM = -1, 0.0, 0.0
for aci in np.linspace(-np.pi, np.pi, 145):
    for egim in np.linspace(-0.4, 0.5, 10):
        Q = (UCLU - UCLU.mean(0)) @ (donme(aci, egim) @ TABAN).T
        alan = abs((Q[1, 0] - Q[0, 0]) * (Q[2, 1] - Q[0, 1]) - (Q[1, 1] - Q[0, 1]) * (Q[2, 0] - Q[0, 0]))
        if alan > en:
            en, ACI, EGIM = alan, aci, egim
T_KAM = donme(ACI, EGIM) @ TABAN
# cerceve: dönüs, heliks 2'nin basi ve Rab10 switch bolgeleri izdusumde ortalanir ve sigdirilir
P = np.vstack([kif0[[S[0].kif_no.index(n) for n in range(1108, 1138)]],
               rab0[[S[0].rab_no.index(n) for n in range(60, 76)]], UCLU])
Q = (P - P.mean(0)) @ T_KAM.T
cx, cy = (Q[:, 0].max() + Q[:, 0].min()) / 2, (Q[:, 1].max() + Q[:, 1].min()) / 2
MERKEZ = P.mean(0) + T_KAM[0] * cx + T_KAM[1] * cy
OLCEK = min(0.92 / ((Q[:, 0].max() - Q[:, 0].min()) / 2 + 2.0), 0.85 / ((Q[:, 1].max() - Q[:, 1].min()) / 2 + 2.0))


def kif_renk(n, sarmal):
    if not sarmal:
        return RENK["sarmal_degil"]
    return RENK["heliks1"] if n <= 1120 else RENK["heliks2"]


def sarmal_mi(sist, n, i):
    if n in sist.heliks:
        return bool(sist.heliks[n][i])
    return 1088 <= n <= 1113 or 1129 <= n <= 1148


def panel(fig, x0, sist, f):
    ax = fig.add_axes([x0, 0.33, 0.32, 0.53])
    bos_eksen(ax, (-1, 1), (-0.932, 0.932))
    ax.set_facecolor(ARKA)
    i = int(round(f))
    X = sist.kare(f)
    s = Sahne3B(T_KAM, MERKEZ, OLCEK)
    s.ekle_tup(X[sist.ca_rab], [RENK["switch"] if (36 <= n <= 44 or 66 <= n <= 80) else RENK["rab10"]
                                for n in sist.rab_no], [4.0] * len(sist.ca_rab), n=5)
    sar = [sarmal_mi(sist, n, i) for n in sist.kif_no]
    s.ekle_tup(X[sist.ca_kif], [kif_renk(n, h) for n, h in zip(sist.kif_no, sar)],
               [7.5 if h else 5.0 for h in sar], n=5)
    for key in ("70", "1130", "1122"):
        ids = sist.yan[key]
        for u, ia in enumerate(ids):
            for ib in ids[u + 1:]:
                if np.linalg.norm(X[ia] - X[ib]) < 1.9:
                    s.ekle_bag(X[ia], X[ib], "#aab2bd", 2.6)
        for ia in ids:
            x = sist.atomlar[ia]
            s.ekle_kure(X[ia], ELEMENT.get(x["element"][0], "#c8ccd2"), 120 if x["ad"] in YUKLU else 45)
    oks = [X[ia] for key in ("1122", "1130") for ia in sist.yan[key] if sist.atomlar[ia]["ad"] in {"OE1", "OE2"}]
    for ia in sist.na:
        if oks and min(np.linalg.norm(X[ia] - o) for o in oks) < 6.0:
            s.ekle_kure(X[ia], RENK["na"], 190)
    s.ciz(ax, katman=12)
    ax.text(0, 0.86, sist.ad, ha="center", fontsize=19, fontweight="bold", color=sist.renk)
    h1124 = sarmal_mi(sist, 1124, i)
    fig.text(x0 + 0.16, 0.305, "1124: " + ("sarmal" if h1124 else "çözülmüş"), ha="center", fontsize=13.5,
             color=RENK["heliks2"] if h1124 else YAZI, fontweight="bold")


def grafik(fig, t_simdi, g=1.0):
    ax = fig.add_axes([0.07, 0.135, 0.88, 0.135])
    ax.set_facecolor(ARKA)
    for k in ("top", "right"):
        ax.spines[k].set_visible(False)
    ax.grid(color=IZGARA, lw=0.7)
    for sist in S:
        m = sist.tc <= t_simdi
        ax.plot(sist.tc[m], sist.donus[m], color=sist.renk, lw=2.0, alpha=g)
    ax.axvline(t_simdi, color=SOLUK, lw=1, ls=":")
    ax.set_xlim(10, 100)
    ax.set_ylim(0, 0.62)
    ax.set_yticks([0, 0.3, 0.6])
    ax.set_yticklabels(["0", "0,3", "0,6"], fontsize=10)
    ax.tick_params(axis="x", labelsize=10)
    ax.set_title("dönüş bölgesi 1116-1126 sarmal oranı (2 ns kayan ortalama); yatay eksen zaman (ns)", fontsize=11,
                 color=SOLUK, pad=3, loc="left")


def s_giris(fig, p):
    baslik(fig, "MD: aynı başlangıç, üç dizi", "Rab10-GTP-Mg + KIF13A 1085-1160, 0,15 M NaCl, 310 K, 100 ns",
           gorunurluk=aralik(p, 0, 0.2), y=0.66)
    fig.text(0.04, 0.46, "Üç koşu birebir aynı yapıdan başlıyor; tek fark 1122. kalıntı.\n"
             "WT: lizin (+).  K1122R: arjinin (+), yükü koruyan kontrol.  K1122E: glutamat (−).\n"
             "Bütün kareler Rab10'a hizalı; kamera sabit, dönüşe ve Rab10 R70'e bakıyor.",
             fontsize=15, color=SOLUK, alpha=aralik(p, 0.2, 0.45), linespacing=1.6)


def s_md(fig, p):
    f = p * (N - 1)
    i = int(round(f))
    t = S[0].t[i]
    baslik(fig, "Aynı başlangıçtan üç koşu", None, y=0.965)
    fig.text(0.96, 0.935, f"t = {t:5.1f} ns".replace(".", ","), ha="right", fontsize=20, color=YAZI,
             fontweight="bold", family="DejaVu Sans Mono")
    for k, sist in enumerate(S):
        panel(fig, 0.01 + 0.33 * k, sist, f)
    grafik(fig, t)
    fig.text(0.04, 0.885, "■ sarmal (heliks 1)", color=RENK["heliks1"], fontsize=11.5)
    fig.text(0.17, 0.885, "■ sarmal (heliks 2)", color=RENK["heliks2"], fontsize=11.5)
    fig.text(0.30, 0.885, "■ sarmal değil", color=RENK["sarmal_degil"], fontsize=11.5)
    fig.text(0.40, 0.885, "● Na+", color=RENK["na"], fontsize=11.5)
    fig.text(0.46, 0.885, "● N (+)", color=RENK["arti"], fontsize=11.5)
    fig.text(0.53, 0.885, "● O (−)", color=RENK["eksi"], fontsize=11.5)
    altyazi_akisi(fig, p, [
        (0.0, 0.22, "Üç koşu aynı başlangıçtan. Ortada dönüş, altta Rab10 (mavi) ve onun R70'i; turuncu ve kiremit "
                    "sarmal, beyaz sarmal değil."),
        (0.22, 0.46, "İlk 40 ns'de üçünde de dönüş benzer. Alttaki eğriler dönüş bölgesinin sarmal oranı."),
        (0.46, 0.76, "K1122E'de yaklaşık 50 ns'den sonra heliks 2'nin ilk dönüşü (1124-1125) çözülüyor. WT ve "
                     "K1122R'de dönüş 100 ns boyunca aynı kalıyor."),
        (0.76, 1.0, "Rab10'a değen heliks 2'nin ikinci yarısı yerinde kalıyor; ara yüz kırılmıyor. Mor küreler "
                    "E1122'nin çevresinde tutulan Na+ iyonları."),
    ], y=0.012, boyut=14.5, genislik=118)


def v(x):
    return f"{x:.2f}".replace(".", ",")


def s_ozet(fig, p):
    sonlar = {sist.ad: np.mean(sist.heliks[1124]) for sist in S}
    kart(fig, "Bu başlangıçta ve bütün koşularda", [
        "Bu başlangıçta 1124 sarmal oranı: WT " + v(sonlar["WT"]) + ", K1122R " + v(sonlar["K1122R"])
        + ", K1122E " + v(sonlar["K1122E"]) + ".",
        "12 başlangıç yapısının hepsinde: kontroller 0,89-0,96, mutantlar 0,20-0,82 (tam ayrışma, p=0,002).",
        "Bazı koşularda çözülme 30-50 ns'den sonra başlıyor; 100 ns olayın son boyunu göstermiyor.",
        "Rab10 ara yüzünün çekirdeği, toplam temas ve gömülü yüzey değişmiyor.",
    ], p)


h = 0.6 if a.hizli else 1.0
sahneler = [Sahne(6 * h, s_giris), Sahne(30 * h, s_md), Sahne(8 * h, s_ozet)]
cikti = KOK / "figures" / "videos" / "03_md.mp4"
n = video_yaz(sahneler, cikti, fps=12 if a.hizli else 24)
print("kare", n, cikti)
gif_yap(cikti, cikti.with_suffix(".gif"), bas=12 * h, sure=24 * h, fps=6, genislik=640, renk=96)
print("gif", cikti.with_suffix(".gif"))
