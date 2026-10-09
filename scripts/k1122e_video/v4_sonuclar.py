"""Video 4: Butun kosular: K1122E neyi degistiriyor, neyi degistirmiyor?

Girdi: results/k1122e_md_20261009/ (kompleks: tum_kosular_20261009.tsv, tum_kosular_bloklar_20261009.tsv,
       kosular/*/*/{analysis.json, interface_residues.tsv}; monomer: monomer_tam_100ns_0910.tsv)
Istatistik birimi baslangic yapisi: ayni kol ve baslangicin rig ve TRUBA kopyalarinin ortalamasi.
Cikti: figures/videos/04_sonuclar.mp4 ve 04_sonuclar.gif
Kullanim: python scripts/k1122e_video/v4_sonuclar.py [--hizli]
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from matplotlib.ticker import FuncFormatter
from scipy import stats

from ortak import (ARKA, IZGARA, RENK, SOLUK, YAZI, Sahne, altyazi_akisi, aralik, baslik, gif_yap, kart,
                   video_yaz)

KOK = Path(__file__).resolve().parents[2]
ap = argparse.ArgumentParser()
ap.add_argument("--hizli", action="store_true")
a = ap.parse_args()
R = KOK / "results" / "k1122e_md_20261009"
KOLLAR = ["WT", "K1122R_paired", "K1122E_paired", "K1122E_AF3"]
KOL_AD = {"WT": "WT", "K1122R_paired": "K1122R", "K1122E_paired": "K1122E\n(WT yapısından)",
          "K1122E_AF3": "K1122E\n(AF3 modeli)"}
KOL_RENK = {"WT": RENK["WT"], "K1122R_paired": RENK["K1122R"], "K1122E_paired": RENK["K1122E"],
            "K1122E_AF3": "#ff9f9f"}
KONTROL = {"WT", "K1122R_paired"}


def tsv(p):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))


def kol_bas(kosu):
    k, s = kosu.rsplit("_s", 1)
    return k.replace("_wt", "").replace("_k1122e", ""), "s" + s


# birim duzeyi degerler
kos = tsv(R / "kompleks" / "tum_kosular_20261009.tsv")
birim = defaultdict(lambda: defaultdict(list))
for r in kos:
    u = kol_bas(r["kosu"])
    for k in ("H1124", "temas cifti", "gomulu SASA nm2", "R70-E1130"):
        birim[u][k].append(float(r[k]))
    kd = R / "kompleks" / "kosular" / r["kaynak"] / r["kosu"]
    aj = json.loads((kd / "analysis.json").read_text())
    if isinstance(aj.get("gbsa"), dict):
        birim[u]["GBSA"].append(aj["gbsa"]["mean_kj_mol"])
    occ = {}
    for x in tsv(kd / "interface_residues.tsv"):
        m = re.match(r"KIF13A [A-Z](\d+)", x["residue"])
        if m:
            occ[int(m.group(1))] = float(x["occupancy"])
    for n in (1089, 1092, 1106, 1110, 1118):
        birim[u][f"r{n}"].append(occ.get(n, 0.0))
UNITS = sorted(birim, key=lambda u: (KOLLAR.index(u[0]), u[1]))


def deger(k):
    return {u: float(np.mean(birim[u][k])) for u in UNITS if birim[u][k]}


def pdeg(k):
    v = deger(k)
    c = [v[u] for u in v if u[0] in KONTROL]
    m = [v[u] for u in v if u[0] not in KONTROL]
    return stats.mannwhitneyu(m, c, alternative="two-sided", method="exact").pvalue


def virgul(x, _=None):
    return f"{x:g}".replace(".", ",").replace("-", "−")


def eksen_sus(ax):
    ax.set_facecolor(ARKA)
    for k in ("top", "right"):
        ax.spines[k].set_visible(False)
    ax.grid(axis="y", color=IZGARA, lw=0.8)
    ax.set_axisbelow(True)
    ax.yaxis.set_major_formatter(FuncFormatter(virgul))


def serit(ax, k, g_bas, g_bit, p, boy=70, etiket=True):
    v = deger(k)
    sira = [u for u in UNITS if u in v]
    for j, u in enumerate(sira):
        g = aralik(p, g_bas + (g_bit - g_bas) * j / len(sira), g_bas + (g_bit - g_bas) * (j + 1) / len(sira))
        if g <= 0:
            continue
        x = KOLLAR.index(u[0]) + {"s1": -0.18, "s2": -0.06, "s3": 0.06, "s5": 0.18}[u[1]]
        ax.scatter([x], [v[u]], s=boy, color=KOL_RENK[u[0]], alpha=g, edgecolors="#07090d", zorder=3)
        if etiket and g > 0.9:
            ax.text(x + 0.05, v[u], u[1], fontsize=9, color=SOLUK, va="center")
    ax.set_xticks(range(4))
    ax.set_xticklabels([KOL_AD[k] for k in KOLLAR], fontsize=11)
    ax.set_xlim(-0.6, 3.6)


def s_baslik(fig, p):
    baslik(fig, "Bütün koşular: ne değişti, ne değişmedi?",
           "Kompleks MD: 12 başlangıç yapısı, 15 yörünge, her biri 100 ns. Monomer MD: 8 koşu.",
           gorunurluk=aralik(p, 0, 0.25), y=0.62)
    fig.text(0.04, 0.42, "İstatistik birimi başlangıç yapısı: aynı yapıdan koşulan kopyaların ortalaması alınır.\n"
             "Kontrol = WT + K1122R (6 birim), mutant = iki K1122E kolu (6 birim). Test: tam Mann-Whitney.",
             fontsize=15, color=SOLUK, alpha=aralik(p, 0.3, 0.55), linespacing=1.6)


def s_1124(fig, p):
    baslik(fig, "1124: heliks 2'nin ilk kalıntısı", "Sarmal oranı (DSSP), her nokta bir başlangıç yapısı")
    ax = fig.add_axes([0.08, 0.25, 0.86, 0.53])
    eksen_sus(ax)
    v = deger("H1124")
    cmin = min(x for u, x in v.items() if u[0] in KONTROL)
    mmax = max(x for u, x in v.items() if u[0] not in KONTROL)
    g = aralik(p, 0.55, 0.7)
    ax.axhspan(mmax, cmin, color="#ffffff", alpha=0.06 * g, lw=0)
    ax.text(3.55, (mmax + cmin) / 2, "boşluk", ha="right", va="center", fontsize=11, color=SOLUK, alpha=g)
    serit(ax, "H1124", 0.05, 0.5, p, boy=110)
    ax.set_ylim(0, 1.05)
    ax.text(1.5, 0.08, "tam ayrışma, p=0,002", ha="center", fontsize=16, fontweight="bold", color=YAZI,
            alpha=aralik(p, 0.6, 0.72))
    altyazi_akisi(fig, p, [
        (0.0, 0.5, "WT ve K1122R'de 1124 hep sarmal (0,89-0,96). K1122E'nin iki kolunda da düşük (0,20-0,82)."),
        (0.5, 1.0, "Altı mutant başlangıcın hepsi altı kontrolün altında. Bu düzende en küçük olası p, 0,002. "
                   "Yükü koruyan K1122R etkiyi göstermiyor; etki yükün ters dönmesine bağlı."),
    ])


def s_monomer(fig, p):
    baslik(fig, "Rab10 olmadan: yalnız saç tokası", "Monomer MD, 1092-1146, her nokta bir koşu (10-100 ns)")
    ax = fig.add_axes([0.08, 0.25, 0.86, 0.53])
    eksen_sus(ax)
    mon = tsv(R / "monomer" / "monomer_tam_100ns_0910.tsv")
    sirala = {"WT": 0, "K1122R": 1, "K1122E": 2}
    renk = {"WT": RENK["WT"], "K1122R": RENK["K1122R"], "K1122E": RENK["K1122E"]}
    mon.sort(key=lambda r: sirala[r["gorev"].split("_")[1]])
    for j, r in enumerate(mon):
        var = r["gorev"].split("_")[1]
        g = aralik(p, 0.05 + 0.05 * j, 0.1 + 0.05 * j)
        x = sirala[var] + (int(r["gorev"][-1]) - 2) * 0.12
        ax.scatter([x], [float(r["DSSP_H1124"])], s=110, color=renk[var], alpha=g, edgecolors="#07090d", zorder=3)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["WT (3)", "K1122R (2)", "K1122E (3)"], fontsize=12)
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("1124 sarmal oranı", fontsize=12)
    ax.text(1, 0.12, "5 kontrole karşı 3 mutant: tam ayrışma, p=0,036", ha="center", fontsize=15,
            fontweight="bold", color=YAZI, alpha=aralik(p, 0.5, 0.62))
    altyazi_akisi(fig, p, [
        (0.0, 0.5, "Aynı ölçüm Rab10 olmadan: WT 0,89-0,94, K1122R 0,92-0,95, K1122E 0,60-0,67."),
        (0.5, 1.0, "Etki Rab10'a bağlı değil; saç tokasının kendi özelliği. Saç tokasının toplam sarmallığı "
                   "(önceden seçilmiş birincil ölçüm) ise K1122R ile örtüşüyor."),
    ])


def s_zaman(fig, p):
    baslik(fig, "Zaman içinde dönüş", "Dönüş bölgesi 1116-1126 sarmal oranı, 10 ns'lik bloklar")
    ax = fig.add_axes([0.08, 0.25, 0.86, 0.53])
    eksen_sus(ax)
    ax.grid(axis="x", color=IZGARA, lw=0.8)
    blok = [r for r in tsv(R / "kompleks" / "tum_kosular_bloklar_20261009.tsv")
            if r["pencere"] == "tam" and r["kaynak"] == "truba_tam"]
    kesme = 10 + 90 * aralik(p, 0.05, 0.6)
    for r in blok:
        kol, bas = kol_bas(r["kosu"])
        y = [float(x) for x in r["blok_donus_heliks(10ns)"].split()]
        t = 15 + 10 * np.arange(len(y))
        m = t <= kesme
        mutant = kol not in KONTROL
        ax.plot(t[m], np.array(y)[m], color=KOL_RENK[kol], lw=2.6 if mutant else 1.6,
                ls="--" if kol == "K1122E_AF3" else "-", marker="o", ms=4, alpha=0.95 if mutant else 0.7)
        if mutant and bas == "s3" and m.all() and p > 0.62:
            ax.text(t[-1] + 1.5, y[-1], {"K1122E_paired": "WT yapısından", "K1122E_AF3": "AF3 modeli"}[kol] + " s3",
                    fontsize=10.5, color=KOL_RENK[kol],
                    va="center", alpha=aralik(p, 0.62, 0.7))
    ax.set_xlim(10, 112)
    ax.set_ylim(0, 0.55)
    ax.set_xlabel("zaman (ns)", fontsize=12)
    ax.xaxis.set_major_formatter(FuncFormatter(virgul))
    for i, (ad, rr, ls) in enumerate((("WT", RENK["WT"], "-"), ("K1122R", RENK["K1122R"], "-"),
                                      ("K1122E, WT yapısından", RENK["K1122E"], "-"),
                                      ("K1122E, AF3 modeli", "#ff9f9f", "--"))):
        ax.plot([], [], color=rr, ls=ls, lw=2.2, label=ad)
    ax.legend(loc="lower left", fontsize=10.5, frameon=False, ncol=4, labelcolor=YAZI)
    altyazi_akisi(fig, p, [
        (0.0, 0.5, "12 TRUBA koşusunun 10 ns'lik ortalamaları. Kontroller (mavi, yeşil) baştan sona düz."),
        (0.5, 1.0, "İki s3 koşusunda dönüş 30-50 ns'den sonra çözülüyor ve 100 ns'de hâlâ düşüyor. 100 ns olayın "
                   "son boyunu yakalamıyor; 60 ns'lik ön analiz etkiyi küçük gösteriyordu."),
    ])


def s_arayuz(fig, p):
    baslik(fig, "Rab10 ara yüzü kırılıyor mu?", "Her nokta bir başlangıç yapısı")
    for j, (k, ad) in enumerate((("temas cifti", "ara yüz temas çifti"), ("gomulu SASA nm2", "gömülü yüzey (nm²)"),
                                 ("R70-E1130", "R70-E1130 köprüsü"))):
        ax = fig.add_axes([0.06 + 0.32 * j, 0.27, 0.27, 0.48])
        eksen_sus(ax)
        g0 = 0.05 + 0.15 * j
        serit(ax, k, g0, g0 + 0.25, p, boy=60, etiket=False)
        ax.set_xticklabels(["WT", "K1122R", "K1122E\nWT yap.", "K1122E\nAF3"], fontsize=9.5)
        ax.set_title(ad, fontsize=13, pad=8)
        ax.text(0.5, 0.02, f"p={pdeg(k):.2f}".replace(".", ","), transform=ax.transAxes, ha="center",
                fontsize=13, color=YAZI, fontweight="bold", alpha=aralik(p, g0 + 0.25, g0 + 0.35))
        v = list(deger(k).values())
        ax.set_ylim(min(v) - 0.25 * (max(v) - min(v)) - 0.02, max(v) + 0.1 * (max(v) - min(v)) + 0.02)
    altyazi_akisi(fig, p, [
        (0.0, 0.55, "Temas sayısı, gömülü yüzey ve R70-E1130 köprüsü kollar arasında örtüşüyor."),
        (0.55, 1.0, "Ara yüzün çekirdeği (heliks 2'nin ikinci yarısı) her kolda %95-100 dolu. K1122E, Rab10 ara "
                    "yüzünü kırmıyor."),
    ])


def s_heliks1(fig, p):
    baslik(fig, "Yeni ve temkinli: heliks 1 tarafı", "Rab10'a temas oranı, kalıntı başına")
    for j, n in enumerate((1089, 1092, 1106, 1118, 1110)):
        ax = fig.add_axes([0.05 + 0.19 * j, 0.27, 0.16, 0.48])
        eksen_sus(ax)
        v = deger(f"r{n}")
        g = aralik(p, 0.05 + 0.08 * j, 0.25 + 0.08 * j)
        for u, x in v.items():
            xx = (0 if u[0] in KONTROL else 1) + {"s1": -0.15, "s2": -0.05, "s3": 0.05, "s5": 0.15}[u[1]]
            ax.scatter([xx], [x], s=45, color=KOL_RENK[u[0]], alpha=g, edgecolors="#07090d", zorder=3)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["kontrol", "mutant"], fontsize=10)
        ax.set_xlim(-0.5, 1.5)
        ax.set_ylim(-0.05, 1.05)
        ad = {1089: "E1089", 1092: "N1092", 1106: "R1106", 1118: "N1118", 1110: "D1110"}[n]
        ax.set_title(ad, fontsize=13, pad=6)
        ax.text(0.5, 0.93, f"p={pdeg(f'r{n}'):.3f}".replace(".", ","), transform=ax.transAxes, ha="center",
                fontsize=11, color=YAZI, alpha=g)
    altyazi_akisi(fig, p, [
        (0.0, 0.5, "Heliks 1'in Rab10'a değen kalıntılarının çoğu mutantta daha az değiyor; D1110 tersine daha çok. "
                   "E1089 tam ayrışıyor."),
        (0.5, 1.0, "Ama bunlar düzeltilmemiş p değerleri (43 kalıntı tarandı) ve en güçlü sinyaller parçanın yapay "
                   "kesik uçlarına yakın. Bu yüzden bir hipotez: tam boy ya da daha uzun bir parçayla sınanmalı."),
    ])


def s_gbsa(fig, p):
    baslik(fig, "Kaba bağlanma enerjisi (GBSA)", "Aynı başlangıçta K1122E eksi WT (kJ/mol); pozitif = daha zayıf")
    ax = fig.add_axes([0.2, 0.27, 0.6, 0.48])
    eksen_sus(ax)
    v = deger("GBSA")
    farklar = [(b, v[("K1122E_paired", b)] - v[("WT", b)]) for b in ("s2", "s3", "s5")]
    for j, (b, d) in enumerate(farklar):
        g = aralik(p, 0.05 + 0.12 * j, 0.25 + 0.12 * j)
        ax.bar(j, d * g, color=RENK["K1122E"] if d > 0 else RENK["WT"], width=0.55)
        if g > 0.95:
            ax.text(j, d + (4 if d > 0 else -4), f"{d:+.0f}".replace("-", "−"), ha="center",
                    va="bottom" if d > 0 else "top", fontsize=14, color=YAZI, fontweight="bold")
    ax.axhline(0, color=SOLUK, lw=1)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["başlangıç s2", "başlangıç s3", "başlangıç s5"], fontsize=12)
    ax.set_ylim(-65, 55)
    altyazi_akisi(fig, p, [
        (0.0, 0.5, "Yalnız s2 varken mutant daha zayıf bağlanıyor gibi görünüyordu."),
        (0.5, 1.0, "Üç başlangıçta işaret değişiyor: biri daha zayıf, ikisi daha güçlü. Tutarlı bir yön yok; yük "
                   "değiştiren mutasyonda bu yöntem zaten güvenilir değil."),
    ])


def s_ozet(fig, p):
    kart(fig, "Kısım 2'nin özeti", [
        "Desteklenen: K1122E dönüşü ve heliks 2'nin başını (1124-1125) yerel olarak çözüyor; Rab10'lu ve "
        "Rab10'suz sistemde, K1122R'de değil.",
        "Desteklenmeyen: Rab10 ara yüzünün kırılması, daha sıkı bağlanma, GBSA'da tutarlı zayıflama.",
        "Hipotez: heliks 1 tarafında temasların yeniden dağılması (kesik uçlara yakın, düzeltilmemiş).",
        "Açık: K1104E kontrolü, daha uzun koşular, tam boy bağlam, varyantın kalıtımı.",
    ], p)


h = 0.6 if a.hizli else 1.0
sahneler = [Sahne(5 * h, s_baslik), Sahne(11 * h, s_1124), Sahne(9 * h, s_monomer), Sahne(11 * h, s_zaman),
            Sahne(11 * h, s_arayuz), Sahne(11 * h, s_heliks1), Sahne(9 * h, s_gbsa), Sahne(10 * h, s_ozet)]
cikti = KOK / "figures" / "videos" / "04_sonuclar.mp4"
n = video_yaz(sahneler, cikti, fps=12 if a.hizli else 30)
print("kare", n, cikti)
gif_yap(cikti, cikti.with_suffix(".gif"), bas=5 * h, sure=31 * h, fps=8, genislik=640)
print("gif", cikti.with_suffix(".gif"))
