"""Video 2: K1122E nerede duruyor? Rab10-GTP-Mg + KIF13A 1085-1160 modelinde mutasyon yeri ve ara yuz.

Girdi: results/k1122e_md_20261009/yapi/{WT_s3_baslangic.pdb, K1122E_s3_baslangic.pdb}
       (AF3 WT modelinden kurulan MD baslangic yapilari, agir atomlar; K1122E ayni yapidan kuruldu)
Cikti: figures/videos/02_yapi.mp4 ve 02_yapi.gif
Kullanim: python scripts/k1122e_video/v2_yapi.py [--hizli]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from ortak import (ARKA, RENK, SOLUK, YAZI, Sahne, Sahne3B, altyazi_akisi, aralik, baslik, bos_eksen, donme,
                   gif_yap, kart, pdb_oku, video_yaz)

KOK = Path(__file__).resolve().parents[2]
ap = argparse.ArgumentParser()
ap.add_argument("--hizli", action="store_true")
a = ap.parse_args()
Y = KOK / "results" / "k1122e_md_20261009" / "yapi"
wt = pdb_oku((Y / "WT_s3_baslangic.pdb").read_text())
mt = pdb_oku((Y / "K1122E_s3_baslangic.pdb").read_text())
OMURGA = {"N", "CA", "C", "O", "OXT"}
ELEMENT = {"C": "#c8ccd2", "N": RENK["arti"], "O": RENK["eksi"], "S": "#ffd166"}
YUKLU = {"NZ", "NH1", "NH2", "NE", "OE1", "OE2", "OD1", "OD2"}


def ca(atomlar, zincir):
    d = {x["no"]: x["xyz"] for x in atomlar if x["zincir"] == zincir and x["ad"] == "CA"}
    nos = sorted(d)
    return nos, np.array([d[n] for n in nos])


def yan(atomlar, zincir, no):
    return [x for x in atomlar if x["zincir"] == zincir and x["no"] == no and x["ad"] not in OMURGA]


def atom(atomlar, zincir, no, adlar):
    return [x["xyz"] for x in atomlar if x["zincir"] == zincir and x["no"] == no and x["ad"] in adlar]


rab_no, rab_ca = ca(wt, "A")
kif_no, kif_ca = ca(wt, "B")
gtp = [x for x in wt if x["zincir"] == "C"]
mg = np.array([x["xyz"] for x in wt if x["zincir"] == "D"])


def rab_renk(n):
    return RENK["switch"] if (36 <= n <= 44 or 66 <= n <= 80) else RENK["rab10"]


def kif_renk(n):
    if n <= 1117:
        return RENK["heliks1"]
    if n <= 1123:
        return RENK["donus"]
    return RENK["heliks2"] if n <= 1148 else "#b8a99a"


# kamera: KIF13A ustte, saç tokasinin uzun ekseni yatay
merkez_tum = np.vstack([rab_ca, kif_ca]).mean(0)
yukari = kif_ca.mean(0) - rab_ca.mean(0)
yukari /= np.linalg.norm(yukari)
_, _, vt = np.linalg.svd(kif_ca - kif_ca.mean(0))
uzun = vt[0] - (vt[0] @ yukari) * yukari
uzun /= np.linalg.norm(uzun)
TABAN = np.array([uzun, yukari, np.cross(uzun, yukari)])
R_TUM = np.max(np.linalg.norm(np.vstack([rab_ca, kif_ca]) - merkez_tum, axis=1))
X_YARI, Y_YARI = 1.185, 0.5
donus_merkez = np.vstack([kif_ca[[kif_no.index(n) for n in range(1118, 1132)]],
                          rab_ca[[rab_no.index(70)]]]).mean(0)
cekirdek = [1134, 1137, 1138, 1141, 1144, 1145, 1148]
rab_yuz = [42, 44, 46, 63, 70, 71, 74, 78]
ara_merkez = np.vstack([kif_ca[[kif_no.index(n) for n in cekirdek]],
                        rab_ca[[rab_no.index(n) for n in rab_yuz]]]).mean(0)


def sahne_kur(aci, merkez, olcek, egim=0.18):
    T = donme(aci, egim) @ TABAN
    s = Sahne3B(T, merkez, olcek)
    s.ekle_tup(rab_ca, [rab_renk(n) for n in rab_no], [5.0 if 36 <= n <= 44 or 66 <= n <= 80 else 3.8 for n in rab_no])
    s.ekle_tup(kif_ca, [kif_renk(n) for n in kif_no], [7.5 if kif_renk(n) != RENK["donus"] else 5.5 for n in kif_no])
    for x in gtp:
        s.ekle_kure(x["xyz"], RENK["gtp"], 22 * olcek / 0.012)
    s.ekle_kure(mg, RENK["mg"], 70 * olcek / 0.012)
    return s


UCLU = np.array([atom(wt, "B", 1122, {"NZ"})[0], np.mean(atom(wt, "B", 1130, {"OE1", "OE2"}), axis=0),
                 np.mean(atom(wt, "A", 70, {"NH1", "NH2", "NE"}), axis=0)])
UCLU_MERKEZ = UCLU.mean(0)


def en_iyi_aci():
    """K1122, E1130 ve R70'in izdusumde en genis ucgeni yaptigi bakis (dikey egim sinirli)."""
    en, sec = -1, (0.0, 0.18)
    for aci in np.linspace(-np.pi, np.pi, 145):
        for egim in np.linspace(-0.5, 0.6, 12):
            Q = (UCLU - UCLU_MERKEZ) @ (donme(aci, egim) @ TABAN).T
            alan = abs(np.cross(Q[1, :2] - Q[0, :2], Q[2, :2] - Q[0, :2]))
            if alan > en:
                en, sec = alan, (aci, egim)
    return sec


ACI_YAKIN, EGIM_YAKIN = en_iyi_aci()


def yan_ekle(s, atomlar, zincir, no, olcek, alfa=1.0, tek_renk=None, buyuk=1.0):
    ys = yan(atomlar, zincir, no)
    for i, p in enumerate(ys):
        for q in ys[i + 1:]:
            if np.linalg.norm(p["xyz"] - q["xyz"]) < 1.9:
                s.ekle_bag(p["xyz"], q["xyz"], "#aab2bd", 3.2 * buyuk, alfa)
    for p in ys:
        renk = tek_renk or ELEMENT.get(p["element"][0], "#c8ccd2")
        boy = (190 if p["ad"] in YUKLU else 80) * buyuk * (olcek / 0.039) ** 1.5
        s.ekle_kure(p["xyz"], renk, boy, alfa=alfa)


def etiket(ax, s, p, metin, renk, dx=0.0, dy=0.06, g=1.0, boyut=13, disari=None):
    q = s.izdusum(p)[0]
    if disari is not None:
        o = s.izdusum(disari)[0]
        v = q[:2] - o[:2]
        v = v / (np.linalg.norm(v) + 1e-9)
        dx, dy = 0.36 * v[0], 0.24 * v[1]
    ax.annotate(metin, (q[0], q[1]), (q[0] + dx, q[1] + dy), color=renk, fontsize=boyut, fontweight="bold",
                ha="center", alpha=g, zorder=100,
                arrowprops=dict(arrowstyle="-", color=renk, alpha=0.7 * g, lw=1.0))


def kesikli(ax, s, p1, p2, renk, metin, g=1.0, boyut=11.5, kaydir=(0.02, 0.0)):
    q1, q2 = s.izdusum(p1)[0], s.izdusum(p2)[0]
    ax.plot([q1[0], q2[0]], [q1[1], q2[1]], ls=(0, (3, 2)), color=renk, lw=2.0, alpha=g, zorder=99)
    ax.text((q1[0] + q2[0]) / 2 + kaydir[0], (q1[1] + q2[1]) / 2 + kaydir[1], metin, color=renk, fontsize=boyut, alpha=g, zorder=100,
            fontweight="bold")


def en_yakin(A, B):
    A, B = np.array(A), np.array(B)
    d = np.linalg.norm(A[:, None] - B[None], axis=2)
    i, j = np.unravel_index(d.argmin(), d.shape)
    return A[i], B[j], d[i, j]


def eksen(fig):
    ax = fig.add_axes([0.0, 0.14, 1.0, 0.75])
    bos_eksen(ax, (-X_YARI, X_YARI), (-Y_YARI, Y_YARI))
    ax.set_facecolor(ARKA)
    return ax


def s_genel(fig, p):
    baslik(fig, "K1122E nerede duruyor?", "Rab10 (GTP bağlı) ve KIF13A'nın Rab10 bağlama bölgesi, 1085-1160")
    ax = eksen(fig)
    olcek = 0.47 / R_TUM
    s = sahne_kur(2 * np.pi * p, merkez_tum, olcek)
    s.ciz(ax)
    g = aralik(p, 0.08, 0.2)
    etiket(ax, s, rab_ca.mean(0), "Rab10", RENK["switch"], dx=-0.28, dy=-0.12, g=g, boyut=15)
    etiket(ax, s, kif_ca.mean(0), "KIF13A saç tokası", RENK["heliks1"], dx=0.30, dy=0.17, g=g, boyut=15)
    etiket(ax, s, np.array([x["xyz"] for x in gtp]).mean(0), "GTP + Mg", RENK["gtp"], dx=0.25, dy=-0.17, g=g)
    altyazi_akisi(fig, p, [
        (0.0, 0.5, "Model: AlphaFold 3 ile kurulan Rab10-KIF13A kompleksi, MD'nin başlangıç yapısı. Mavi: Rab10; "
                   "açık mavi: switch I ve switch II, GTP'nin durumunu okuyan bölgeler."),
        (0.5, 1.0, "KIF13A parçası bir saç tokası: heliks 1 (turuncu), dönüş (beyaz) ve heliks 2 (kiremit). "
                   "Rab10'a heliks 2 ile oturuyor."),
    ])


def s_yakin(fig, p, mutant=False):
    ax = eksen(fig)
    olcek0, olcek1 = 0.47 / R_TUM, 0.47 / 11.0
    if not mutant:
        baslik(fig, "Dönüşe yaklaşalım", "K1122, E1130 ve Rab10 R70")
        z = aralik(p, 0.0, 0.35)
    else:
        baslik(fig, "K1122E: artı yük eksiye dönüyor", "Aynı yapıda 1122'nin yan zinciri değişiyor")
        z = 1.0
    merkez = (1 - z) * merkez_tum + z * UCLU_MERKEZ
    olcek = olcek0 * (olcek1 / olcek0) ** z
    salin = 0.12 * np.sin(2 * np.pi * p) * z
    s = sahne_kur(z * ACI_YAKIN + salin, merkez, olcek, egim=(1 - z) * 0.18 + z * EGIM_YAKIN)
    gk = aralik(p, 0.3, 0.45) if not mutant else 1.0
    m = aralik(p, 0.1, 0.35) if mutant else 0.0
    if gk > 0:
        yan_ekle(s, wt, "A", 70, olcek, alfa=gk)
        yan_ekle(s, wt, "B", 1130, olcek, alfa=gk)
        if 1 - m > 0.02:
            yan_ekle(s, wt, "B", 1122, olcek, alfa=gk * (1 - m))
        if m > 0.02:
            yan_ekle(s, mt, "B", 1122, olcek, alfa=m)
    s.ciz(ax)
    ge = aralik(p, 0.42, 0.55) if not mutant else 1.0
    k_nz, e_oe, d1 = en_yakin(atom(wt, "B", 1122, {"NZ"}), atom(wt, "B", 1130, {"OE1", "OE2"}))
    e2, r_nh, d2 = en_yakin(atom(wt, "B", 1130, {"OE1", "OE2"}), atom(wt, "A", 70, {"NE", "NH1", "NH2"}))
    if not mutant:
        etiket(ax, s, k_nz, "K1122 (+)", RENK["arti"], g=ge, disari=UCLU_MERKEZ)
        etiket(ax, s, e_oe, "E1130 (\u2212)", RENK["eksi"], g=ge, disari=UCLU_MERKEZ)
        etiket(ax, s, r_nh, "R70 (+), Rab10", RENK["switch"], g=ge, disari=UCLU_MERKEZ)
        gb = aralik(p, 0.55, 0.68)
        kesikli(ax, s, k_nz, e_oe, YAZI, f"{d1:.1f} Å".replace(".", ","), g=gb)
        kesikli(ax, s, e2, r_nh, YAZI, f"{d2:.1f} Å".replace(".", ","), g=gb)
        altyazi_akisi(fig, p, [
            (0.0, 0.4, "Saç tokasının dönüşüne yaklaşıyoruz. Mutasyonun yeri, 1122, bu dönüşte."),
            (0.4, 0.72, "WT'de K1122'nin artı yükü, heliks 2'deki E1130'un eksi yükünü tutuyor. E1130 da Rab10'un "
                        "switch II bölgesindeki R70 ile tuz köprüsü kuruyor."),
            (0.72, 1.0, "Yani K1122 ara yüzde değil, ikinci kabukta: Rab10'a doğrudan değil, E1130 üzerinden "
                        "bağlı. K1122 ile R70 arasında yaklaşık 4,5 Å var ve ikisi de artı yüklü."),
        ])
    else:
        e_oe_m = atom(mt, "B", 1122, {"OE1", "OE2"})
        a1, b1, d3 = en_yakin(e_oe_m, atom(mt, "B", 1130, {"OE1", "OE2"}))
        etiket(ax, s, k_nz if m < 0.5 else a1, "K1122 (+)" if m < 0.5 else "E1122 (\u2212)",
               RENK["arti"] if m < 0.5 else RENK["eksi"], disari=UCLU_MERKEZ)
        etiket(ax, s, e_oe, "E1130 (\u2212)", RENK["eksi"], disari=UCLU_MERKEZ)
        etiket(ax, s, r_nh, "R70 (+), Rab10", RENK["switch"], disari=UCLU_MERKEZ)
        gi = aralik(p, 0.4, 0.55)
        kesikli(ax, s, a1, b1, RENK["eksi"], "aynı yük: itme", g=gi, kaydir=(-0.36, -0.07))
        altyazi_akisi(fig, p, [
            (0.0, 0.42, "K1122E'de lizin (+) yerine glutamat (\u2212) geliyor."),
            (0.42, 1.0, "E1122 artık E1130'u tutamaz; iki eksi yük birbirini iter. Suda bu çift için ölçtüğümüz "
                        "serbest enerji itici: +3,6 ile +5 kJ/mol. Soru, bunun dönüşü ve Rab10 ara yüzünü ne kadar "
                        "etkilediği."),
        ])


def yz_merkez(atomlar, zincir, no):
    ys = yan(atomlar, zincir, no)
    if ys:
        return np.mean([x["xyz"] for x in ys], axis=0)
    nos, P = ca(atomlar, zincir)
    return P[nos.index(no)]


def s_arayuz(fig, p):
    baslik(fig, "Asıl ara yüz 1122'den uzakta", "Heliks 2'nin ikinci yarısı ve Rab10'un korunmuş yüzeyi")
    ax = eksen(fig)
    z = aralik(p, 0.0, 0.3)
    merkez = (1 - z) * UCLU_MERKEZ + z * ara_merkez
    olcek1, olcek2 = 0.47 / 11.0, 0.47 / 26.0
    olcek = olcek1 * (olcek2 / olcek1) ** z
    s = sahne_kur((1 - z) * ACI_YAKIN + 0.2 * np.sin(np.pi * p) * z, merkez, olcek,
                  egim=(1 - z) * EGIM_YAKIN + z * 0.25)
    g = aralik(p, 0.25, 0.4)
    boy = 330 * (olcek / 0.018) ** 1.2
    if g > 0:
        for n in rab_yuz:
            s.ekle_kure(yz_merkez(wt, "A", n), RENK["switch"], boy, alfa=g)
        for n in cekirdek:
            s.ekle_kure(yz_merkez(wt, "B", n), RENK["heliks2"], boy, alfa=g)
        s.ekle_kure(yz_merkez(wt, "B", 1122), RENK["arti"], boy, alfa=g)
    s.ciz(ax)
    ge = aralik(p, 0.4, 0.52)
    o = np.vstack([yz_merkez(wt, "B", n) for n in cekirdek] + [yz_merkez(wt, "A", n) for n in rab_yuz]).mean(0)
    for zn, n, ad, rr in (("A", 46, "F46", RENK["switch"]), ("A", 63, "W63", RENK["switch"]),
                          ("A", 78, "Y78", RENK["switch"]), ("B", 1137, "W1137", RENK["heliks2"]),
                          ("B", 1144, "R1144", RENK["heliks2"]), ("B", 1122, "K1122", RENK["arti"])):
        etiket(ax, s, yz_merkez(wt, zn, n), ad, rr, g=ge, boyut=12, disari=o)
    altyazi_akisi(fig, p, [
        (0.0, 0.45, "Rab10'a sürekli değen kalıntılar heliks 2'nin ikinci yarısında: V1134, W1137, V1138, T1141, "
                    "R1144, N1145, L1148 (kiremit küreler). Bunlar 1122'den 15-36 Å uzakta."),
        (0.45, 1.0, "Rab10 tarafında (açık mavi küreler) switch I, interswitch ve switch II var; korunmuş "
                    "F46/W63/Y78 üçlüsü Rab efektörlerinin çoğunun tutunduğu yüzey. KIF13A tam bu yüzeye oturuyor."),
    ])


def s_ozet(fig, p):
    kart(fig, "Yapıdan çıkan beklenti", [
        "K1122 dönüşte, ikinci kabukta: E1130'u tutuyor, E1130 da Rab10 R70'i.",
        "K1122E bu tutamağı kaldırıyor ve dönüşe ikinci bir eksi yük getiriyor.",
        "Asıl ara yüz 15-36 Å uzakta; doğrudan bir kopma beklenmez.",
        "Beklenen etki dönüşün kendisinde: MD bunu sınıyor (bir sonraki video).",
    ], p)


h = 0.6 if a.hizli else 1.0
sahneler = [Sahne(10 * h, s_genel), Sahne(12 * h, s_yakin), Sahne(9 * h, lambda f, p: s_yakin(f, p, True)),
            Sahne(11 * h, s_arayuz), Sahne(7 * h, s_ozet)]
cikti = KOK / "figures" / "videos" / "02_yapi.mp4"
n = video_yaz(sahneler, cikti, fps=12 if a.hizli else 24)
print("kare", n, cikti)
gif_yap(cikti, cikti.with_suffix(".gif"), bas=10 * h, sure=21 * h, fps=8, genislik=640)
print("gif", cikti.with_suffix(".gif"))
