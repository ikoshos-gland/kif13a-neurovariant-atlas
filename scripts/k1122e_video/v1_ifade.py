"""Video 1: KIF13A beyinde nerede ifade ediliyor? (Allen Brain Cell Atlas ve insan yasam boyu verisi)

Girdi: kif13a-allenbrain deposunun results/ tablolari (varsayilan: bu deponun yanindaki klon).
Cikti: figures/videos/01_ifade.mp4 ve 01_ifade.gif
Kullanim: python scripts/k1122e_video/v1_ifade.py [--allenbrain ../kif13a-allenbrain] [--hizli]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from ortak import (ARKA, IZGARA, RENK, SOLUK, YAZI, Sahne, altyazi_akisi, aralik, baslik, gif_yap, kart,
                   video_yaz)

KOK = Path(__file__).resolve().parents[2]
ap = argparse.ArgumentParser()
ap.add_argument("--allenbrain", type=Path, default=KOK.parent / "kif13a-allenbrain")
ap.add_argument("--hizli", action="store_true", help="dusuk kare hizi ile deneme")
a = ap.parse_args()
R = a.allenbrain / "results"

sira = pd.read_csv(R / "derived_wholebrain_subclass_ranking.csv").sort_values("mean_Kif13a", ascending=False)
sira = sira.reset_index(drop=True)
sira["sira"] = np.arange(1, len(sira) + 1)
tespit = pd.read_csv(R / "derived_detection_vs_level.csv").set_index("subclass")
donor = pd.read_csv(R / "Kif13a_AgingMouse_OPC_Oligo_donor_summary.csv")
yas = pd.read_csv(R / "derived_donor_level_age_test.csv").set_index("subclass")

TR = {"327 Oligo NN": "olgun oligodendrosit", "326 OPC NN": "OPC (öncül)", "145 MH Tac2 Glut": "medial habenula, kolinerjik",
      "334 Microglia NN": "mikroglia", "319 Astro-TE NN": "astrosit", "331 Peri NN": "perisit"}


def renk(sc: str) -> str:
    if sc.startswith("327 Oligo"):
        return RENK["oligo"]
    if sc.startswith("326 OPC"):
        return RENK["opc"]
    return RENK["glia"] if sc.endswith(" NN") else RENK["noron"]


def virgul(v, _=None):
    return f"{v:g}".replace(".", ",").replace("-", "−")


def eksen_sus(ax):
    from matplotlib.ticker import FuncFormatter
    ax.xaxis.set_major_formatter(FuncFormatter(virgul))
    ax.yaxis.set_major_formatter(FuncFormatter(virgul))
    ax.set_facecolor(ARKA)
    for k in ("top", "right"):
        ax.spines[k].set_visible(False)
    ax.grid(axis="x", color=IZGARA, lw=0.8)
    ax.set_axisbelow(True)


def s_baslik(fig, p):
    g = aralik(p, 0, 0.3)
    baslik(fig, "KIF13A beyinde nerede ifade ediliyor?",
           "Fare beyni: Allen Brain Cell Atlas, 2,34 milyon hücre, 265 hücre alt sınıfı", gorunurluk=g, y=0.62)
    fig.text(0.04, 0.42, "Ölçü: log2(CPM+1). CPM, bir hücredeki mRNA okumalarının milyonda kaçının\n"
             "KIF13A'dan geldiğini söyler. Bu bir mRNA ölçüsüdür, protein değil.",
             fontsize=15, color=SOLUK, alpha=aralik(p, 0.3, 0.6), linespacing=1.5)


def s_ilk15(fig, p):
    baslik(fig, "En yüksek 15 hücre alt sınıfı", "Çubuk: o hücre tipindeki ortalama KIF13A mRNA'sı")
    ax = fig.add_axes([0.33, 0.24, 0.60, 0.56])
    eksen_sus(ax)
    ust = sira.head(15).iloc[::-1]
    for j, (_, r) in enumerate(ust.iterrows()):
        i = 14 - j
        g = aralik(p, 0.03 + i * 0.025, 0.2 + i * 0.025)
        ax.barh(j, r.mean_Kif13a * g, color=renk(r.subclass), height=0.72)
        if g > 0.95:
            ax.text(r.mean_Kif13a + 0.08, j, f"{r.mean_Kif13a:.2f}".replace(".", ","), va="center", fontsize=10.5,
                    color=YAZI)
        if r.subclass in TR and g > 0.95:
            ax.text(r.mean_Kif13a + 0.65, j, TR[r.subclass], va="center", fontsize=12, color=renk(r.subclass),
                    fontweight="bold")
    ax.set_yticks(range(15))
    ax.set_yticklabels(ust.subclass, fontsize=10.5)
    ax.set_xlim(0, 10.5)
    ax.set_xlabel("ortalama log2(CPM+1)", fontsize=12)
    fig.text(0.70, 0.87, "■ oligodendrosit soyu", color=RENK["oligo"], fontsize=12)
    fig.text(0.70, 0.83, "■ nöron", color=RENK["noron"], fontsize=12)
    fig.text(0.85, 0.87, "■ diğer glia/damar", color=RENK["glia"], fontsize=12)
    altyazi_akisi(fig, p, [
        (0.0, 0.48, "Her çubuk bir hücre alt sınıfı. Uzunluk, o hücrelerdeki ortalama KIF13A mRNA'sı (log2 ölçeği)."),
        (0.48, 1.0, "En üstte oligodendrosit soyu var: olgun oligodendrosit 1., öncülü OPC 2. Glia dışında en yüksek "
                    "grup medial habenulanın kolinerjik nöronları (3.)."),
    ])


def s_hepsi(fig, p):
    baslik(fig, "265 alt sınıfın tamamı", "Soldan sağa: en yüksekten en düşüğe")
    ax = fig.add_axes([0.08, 0.26, 0.86, 0.52])
    eksen_sus(ax)
    ax.grid(axis="y", color=IZGARA, lw=0.8)
    n = int(len(sira) * aralik(p, 0.0, 0.45)) or 1
    gor = sira.head(n)
    ax.scatter(gor.sira, gor.mean_Kif13a, s=16, c=[renk(s) for s in gor.subclass], alpha=0.9, lw=0)
    ax.set_xlim(-3, 270)
    ax.set_ylim(0, 9.3)
    ax.set_xlabel("sıra", fontsize=12)
    ax.set_ylabel("ortalama log2(CPM+1)", fontsize=12)
    g = aralik(p, 0.45, 0.6)
    for sc, (dx, dy) in {"327 Oligo NN": (12, 0.7), "326 OPC NN": (13, 0.3), "145 MH Tac2 Glut": (40, -0.6),
                         "334 Microglia NN": (-62, 1.5), "319 Astro-TE NN": (16, 1.9), "331 Peri NN": (-6, 1.9)}.items():
        r = sira[sira.subclass == sc].iloc[0]
        ax.annotate(f"{TR[sc]} ({int(r.sira)}.)", (r.sira, r.mean_Kif13a), (r.sira + dx, r.mean_Kif13a + dy),
                    color=YAZI, fontsize=12, alpha=g, arrowprops=dict(arrowstyle="-", color=SOLUK, alpha=g))
        ax.scatter([r.sira], [r.mean_Kif13a], s=70, facecolors="none", edgecolors=YAZI, alpha=g, lw=1.4)
    altyazi_akisi(fig, p, [
        (0.0, 0.55, "Glia hücrelerinin çoğu sıralamanın ortasında ya da altında: mikroglia 168., astrosit 178., "
                    "perisit 240."),
        (0.55, 1.0, "Ama bu sıralama iki ayrı şeyi tek bir sayıda birleştiriyor: genin kaç hücrede yakalandığını ve "
                    "yakalandığı hücrede ne kadar olduğunu."),
    ])


def s_tespit(fig, p):
    baslik(fig, "Yakalanma oranı ve düzey ayrı ayrı", "Sıralamadaki farkın iki bileşeni")
    tipler = [("327 Oligo NN", "olgun\noligodendrosit"), ("326 OPC NN", "OPC"), ("334 Microglia NN", "mikroglia"),
              ("319 Astro-TE NN", "astrosit"), ("331 Peri NN", "perisit")]
    renkler = [RENK["oligo"], RENK["opc"], RENK["glia"], RENK["glia"], RENK["glia"]]
    for k, (sutun, ad, x0, ust) in enumerate((("fraction_expressing", "Geni ifade eden hücre oranı", 0.07, 1.0),
                                              ("cpm_among_expressing", "İfade eden hücrede düzey (CPM)", 0.55, 380))):
        ax = fig.add_axes([x0, 0.26, 0.40, 0.50])
        eksen_sus(ax)
        ax.grid(axis="x", visible=False)
        ax.grid(axis="y", color=IZGARA, lw=0.8)
        g = aralik(p, 0.05 + 0.15 * k, 0.35 + 0.15 * k)
        deg = [tespit.loc[t, sutun] for t, _ in tipler]
        ax.bar(range(5), [v * g for v in deg], color=renkler, width=0.66)
        for i, v in enumerate(deg):
            if g > 0.95:
                ax.text(i, v * 1.02, (f"{v:.2f}" if k == 0 else f"{v:.0f}").replace(".", ","), ha="center",
                        va="bottom", fontsize=11, color=YAZI)
        ax.set_xticks(range(5))
        ax.set_xticklabels([ad for _, ad in tipler], fontsize=11)
        ax.set_ylim(0, ust * 1.22)
        ax.set_title(ad, fontsize=14, pad=10)
        oran = tespit.loc["327 Oligo NN", sutun] / tespit.loc["319 Astro-TE NN", sutun]
        go = aralik(p, 0.45 + 0.1 * k, 0.6 + 0.1 * k)
        ust_y = max(deg[0], deg[3]) + ust * 0.10
        ax.plot([0, 0, 3, 3], [deg[0] + ust * 0.06, ust_y, ust_y, deg[3] + ust * 0.06], color=YAZI, alpha=go, lw=1.3)
        ax.text(1.5, ust_y + ust * 0.015, f"yaklaşık {oran:.1f} kat".replace(".", ","),
                ha="center", va="bottom", fontsize=13, color=YAZI, alpha=go, fontweight="bold")
    altyazi_akisi(fig, p, [
        (0.0, 0.42, "Sol: genin yakalandığı hücre oranı. Sağ: yakalandığı hücrede ne kadar olduğu."),
        (0.42, 0.74, "Oligodendrosit, astrositin yaklaşık 1,8 katı hücrede yakalanıyor ve yakalandığında yaklaşık "
                     "2 katı düzeyde. Eski '16 kat' ifadesi, log ortalamalarının farkını kat değişimi sanmaktan "
                     "geliyordu."),
        (0.74, 1.0, "Perisit sıralamada 240. ama ifade ettiğinde astrositten yüksek. Yakalanma oranı hücrenin RNA "
                    "miktarına bağlı; bu yüzden 'en üstte' kesin, 'oligodendrosite özgü' kanıtlanmadı."),
    ])


def s_motor(fig, p):
    baslik(fig, "Kontrol: motor nöron hipotezi", "Motor nöronlar 84 hücre tipi arasında hangi sırada?")
    ax = fig.add_axes([0.12, 0.27, 0.76, 0.48])
    ax.set_facecolor(ARKA)
    ax.axis("off")
    genler = [("Kif5a", 1, RENK["noron"]), ("Kif1a", 6, RENK["noron"]), ("Kif13a", 34, RENK["K1122E"])]
    for j, (gen, rk, rr) in enumerate(genler):
        y = 2 - j
        ax.plot([1, 84], [y, y], color=IZGARA, lw=6, solid_capstyle="round")
        g = aralik(p, 0.08 + 0.12 * j, 0.25 + 0.12 * j)
        x = 84 - (84 - rk) * g
        ax.scatter([x], [y], s=260, color=rr, zorder=3, edgecolors="#07090d")
        ax.text(-3, y, gen, ha="right", va="center", fontsize=16, color=YAZI, fontweight="bold")
        if g > 0.95:
            ax.text(x, y + 0.3, f"{rk}.", ha="center", fontsize=14, color=rr, fontweight="bold")
    ax.text(1, -0.75, "1. sıra (en yüksek)", fontsize=12, color=SOLUK)
    ax.text(84, -0.75, "84. sıra", fontsize=12, color=SOLUK, ha="right")
    ax.set_xlim(-12, 90)
    ax.set_ylim(-1.1, 2.8)
    altyazi_akisi(fig, p, [
        (0.0, 0.5, "Aynı motor nöronlar Kif5a'da 1., Kif1a'da 6. sırada; yani kinezinleri bol taşıyorlar. "
                   "Kif13a'da ise 34. sıradalar."),
        (0.5, 1.0, "Omurilik motor nöronları oligodendrositin log2 ölçeğinde 2,07 altında (p=0,023). Omurilikte de "
                   "en yüksek tip oligodendrosit. Motor nöron hipotezi negatif."),
    ])


def s_yas(fig, p):
    baslik(fig, "Yaşlanma: farede artış, insanda yok", "Her nokta bir fare (donör düzeyinde test)")
    ax = fig.add_axes([0.07, 0.27, 0.46, 0.48])
    eksen_sus(ax)
    ax.grid(axis="x", visible=False)
    ax.grid(axis="y", color=IZGARA, lw=0.8)
    rng = np.random.default_rng(3)
    for k, (sc, ad) in enumerate((("327 Oligo NN", "olgun oligodendrosit"), ("326 OPC NN", "OPC"))):
        for m, (kat, etiket, rr) in enumerate((("adult", "genç", "#8fb3d9"), ("aged", "yaşlı", RENK["oligo"]))):
            x0 = k * 2.6 + m
            v = donor[(donor.subclass_name == sc) & (donor.donor_age_category == kat)].mean_Kif13a.dropna().values
            g = aralik(p, 0.05 + 0.08 * (2 * k + m), 0.25 + 0.08 * (2 * k + m))
            n = int(len(v) * g)
            ax.scatter(x0 + rng.uniform(-0.22, 0.22, len(v))[:n], v[:n], s=20, color=rr, alpha=0.85, lw=0)
            if g > 0.95:
                ax.plot([x0 - 0.3, x0 + 0.3], [v.mean()] * 2, color=YAZI, lw=2.4)
                ax.text(x0, 4.55, etiket, ha="center", fontsize=11.5, color=SOLUK)
        r = yas.loc[sc]
        gp = aralik(p, 0.5, 0.6)
        fark = f"+{r.delta_log2:.3f}".replace(".", ",")
        pm = f"{r.mannwhitney_p:.4f}".replace(".", ",") if r.mannwhitney_p < 0.01 else f"{r.mannwhitney_p:.2f}".replace(".", ",")
        ax.text(k * 2.6 + 0.5, 9.25, f"{ad}\n{fark} log2, p={pm}", ha="center", fontsize=12, color=YAZI, alpha=gp)
    ax.set_xlim(-0.6, 4.2)
    ax.set_ylim(4.3, 10.0)
    ax.set_xticks([])
    ax.set_ylabel("fare başına ortalama log2(CPM+1)", fontsize=11.5)

    bx = fig.add_axes([0.62, 0.27, 0.33, 0.48])
    eksen_sus(bx)
    gb = aralik(p, 0.55, 0.75)
    etki = [("fare, olgun\noligodendrosit", 0.250, None, RENK["oligo"]), ("fare, OPC", 0.131, None, RENK["opc"]),
            ("insan, olgun\noligodendrosit", -0.057, (-0.096, -0.017), RENK["WT"])]
    for j, (ad, d, ga, rr) in enumerate(etki):
        y = 2 - j
        bx.barh(y, d * gb, color=rr, height=0.55)
        if ga and gb > 0.95:
            bx.plot(ga, [y, y], color=YAZI, lw=2)
        bx.text(-0.16, y, ad, ha="right", va="center", fontsize=11.5, color=YAZI)
    bx.axvline(0, color=SOLUK, lw=1)
    bx.set_xlim(-0.15, 0.32)
    bx.set_yticks([])
    bx.set_xlabel("yaşlı eksi genç (log2)", fontsize=11.5)
    altyazi_akisi(fig, p, [
        (0.0, 0.5, "Fare: 53 genç ve 44 yaşlı; her fare tek değer. Olgun oligodendrositte ifade yaşla yaklaşık %19 "
                   "artıyor (p=0,0002). OPC'de artış anlamlı değil."),
        (0.5, 1.0, "İnsan: 335 donörde aynı karşılaştırma artış göstermiyor (-0,057 log2). Farenin etkisini görme gücü "
                   "0,99'un üzerindeydi; yani bu gerçek bir negatif."),
    ])


def s_ozet(fig, p):
    kart(fig, "Kısım 1'in özeti", [
        "Beyinde KIF13A en çok oligodendrosit soyunda: farede olgun oligodendrosit 1., insanda OPC 1.",
        "Glia dışında medial habenulanın kolinerjik nöronları öne çıkıyor; motor nöronlarda düşük.",
        "'16 kat' yanlıştı: fark yakalanmada yaklaşık 1,8, düzeyde yaklaşık 2 kat.",
        "Yaşla artış yalnız farede; insanda yok.",
        "Ölçüm mRNA. Proteinin nerede olduğu için genişletme mikroskopisi (ExM) gibi deneyler gerekir.",
    ], p)


hiz = 0.6 if a.hizli else 1.0
sahneler = [Sahne(4.5 * hiz, s_baslik), Sahne(10 * hiz, s_ilk15), Sahne(9 * hiz, s_hepsi), Sahne(12 * hiz, s_tespit),
            Sahne(9 * hiz, s_motor), Sahne(12 * hiz, s_yas), Sahne(8 * hiz, s_ozet)]
cikti = KOK / "figures" / "videos" / "01_ifade.mp4"
n = video_yaz(sahneler, cikti, fps=12 if a.hizli else 30)
print("kare", n, cikti)
gif_yap(cikti, cikti.with_suffix(".gif"), bas=4.5, sure=19, fps=8, genislik=640)
print("gif", cikti.with_suffix(".gif"))
