"""Rab10 ile KIF13A arasindaki butun ara yuz: hangi kalintilar degiyor, K1122E bunlardan birini degistiriyor mu?

Girdi: her kosunun analyze_complex.py ciktilari (interface_contacts.tsv, interface_residues.tsv,
intermolecular_salt_bridges.tsv). Temas: en yakin agir atomlar <= 0.45 nm, 100 ps arayla, ilk 10 ns atilmis.
Tuz koprusu: yuklu grup atomlari <= 0.40 nm, her 10 ps.
Istatistik birimi baslangic yapisi: ayni kol + ayni baslangicin rig ve TRUBA kopyalarinin ortalamasi alinir
(WT 3, K1122R 3, K1122E paired 3, K1122E AF3 3). Kontrol = WT + K1122R (6), mutant = iki K1122E kolu (6).
Bir kalinti bir koside hic listelenmemisse dolulugu 0 sayilir (tablo yalniz en az bir kez temas edenleri yazar).
"""
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

kok = Path(sys.argv[1])
cikti = Path(sys.argv[2])
cikti.mkdir(parents=True, exist_ok=True)

KOLLAR = ["WT", "K1122R_paired", "K1122E_paired", "K1122E_AF3"]
KONTROL, MUTANT = {"WT", "K1122R_paired"}, {"K1122E_paired", "K1122E_AF3"}
APOLAR = set("AVLIMFWYCP")
POZ, NEG = set("KR"), set("DE")


def kol_baslangic(ad):
    kol, s = ad.rsplit("_s", 1)
    return kol.replace("_wt", "").replace("_k1122e", ""), "s" + s


def ayristir(etiket):
    """'Rab10 M1 N-terminus' -> ('Rab10', 'M', 1); numarasiz uc gruplar (NME, ACE) -> None."""
    mol, geri = etiket.split(" ", 1)
    m = re.match(r"([A-Z])(\d+)", geri)
    if not m:
        return None
    return mol, m.group(1), int(m.group(2))


def bolge_kif(n):
    if n < 1088:
        return "N ucu"
    if n <= 1117:
        return "heliks 1"
    if n <= 1123:
        return "donus"
    if n <= 1148:
        return "heliks 2"
    return "C ucu"


def bolge_rab(n):
    if n <= 12:
        return "N ucu"
    if 36 <= n <= 44:
        return "switch I"
    if 45 <= n <= 65:
        return "interswitch"
    if 66 <= n <= 80:
        return "switch II"
    return "diger"


def tsv(p):
    return list(csv.DictReader(p.open(), delimiter="\t")) if p.exists() else []


# kopya duzeyinde oku
kopyalar = defaultdict(list)  # (kol, baslangic) -> [ {kalinti: occ}, {cift: occ}, {kopru: occ} ]
adlar = {}  # (mol, numara) -> {kol: tek harf}
for f in sorted(kok.glob("*/*/interface_contacts.tsv")):
    d = f.parent
    kol, bas = kol_baslangic(d.name)
    res, cift, kopru = {}, {}, {}
    for r in tsv(d / "interface_residues.tsv"):
        a = ayristir(r["residue"])
        if a and a[0] in {"Rab10", "KIF13A"}:
            res[(a[0], a[2])] = float(r["occupancy"])
            adlar.setdefault((a[0], a[2]), {})[kol] = a[1]
    for r in tsv(f):
        a, b = ayristir(r["receptor"]), ayristir(r["kif13a"])
        if a and b and a[0] == "Rab10":
            cift[(a[2], b[2])] = float(r["occupancy"])
    for r in tsv(d / "intermolecular_salt_bridges.tsv"):
        a, b = ayristir(r["receptor"]), ayristir(r["kif13a"])
        if a and b and a[0] == "Rab10":
            kopru[(a[2], b[2])] = max(kopru.get((a[2], b[2]), 0.0), float(r["occupancy"]))
    kopyalar[(kol, bas)].append((res, cift, kopru, f"{d.parent.name}/{d.name}"))

birimler = sorted(kopyalar, key=lambda k: (KOLLAR.index(k[0]), k[1]))


def ortala(i):
    """Her birim icin kopyalarin ortalamasi (eksik = 0)."""
    out = {}
    for u in birimler:
        sozluk = [k[i] for k in kopyalar[u]]
        anahtar = set().union(*sozluk)
        out[u] = {k: float(np.mean([s.get(k, 0.0) for s in sozluk])) for k in anahtar}
    return out


RES, CIFT, KOPRU = ortala(0), ortala(1), ortala(2)


def ad(mol, n):
    harfler = adlar.get((mol, n), {})
    if mol == "KIF13A" and n == 1122:
        return "K/R/E1122"
    return f"{harfler.get('WT') or next(iter(harfler.values()), '?')}{n}"


def kol_ozet(deger):
    """deger: birim -> sayi. Kol basina 'ort (min-max)'."""
    out = {}
    for kol in KOLLAR:
        v = [deger[u] for u in birimler if u[0] == kol]
        out[kol] = (np.mean(v), min(v), max(v)) if v else (np.nan, np.nan, np.nan)
    return out


def karsilastir(deger):
    k = np.array([deger[u] for u in birimler if u[0] in KONTROL])
    m = np.array([deger[u] for u in birimler if u[0] in MUTANT])
    ayrisma = "mutant hep alti" if m.max() < k.min() else ("mutant hep ustu" if m.min() > k.max() else "")
    p = stats.mannwhitneyu(m, k, alternative="two-sided", method="exact").pvalue if (k.std() + m.std()) > 0 else 1.0
    return k, m, ayrisma, p


def yaz(ad_, satirlar):
    with (cikti / ad_).open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(satirlar[0]), delimiter="\t")
        w.writeheader()
        w.writerows(satirlar)


# 1) kalinti duzeyi, iki taraf
for mol, bolge in (("KIF13A", bolge_kif), ("Rab10", bolge_rab)):
    numaralar = sorted({k[1] for u in birimler for k in RES[u] if k[0] == mol})
    satirlar = []
    for n in numaralar:
        deger = {u: RES[u].get((mol, n), 0.0) for u in birimler}
        if max(deger.values()) < 0.10:
            continue
        oz = kol_ozet(deger)
        k, m, ayr, p = karsilastir(deger)
        # baslica karsi taraf ortaklari (12 birimin ortalamasi)
        ort = defaultdict(float)
        for u in birimler:
            for (rn, kn), o in CIFT[u].items():
                if (mol == "KIF13A" and kn == n) or (mol == "Rab10" and rn == n):
                    ort[(rn, kn)] += o / len(birimler)
        enler = sorted(ort.items(), key=lambda x: -x[1])[:3]
        diger = "Rab10" if mol == "KIF13A" else "KIF13A"
        ortak = ", ".join(f"{ad(diger, rn if mol == 'KIF13A' else kn)} {o:.2f}" for (rn, kn), o in enler if o >= 0.05)
        kopru = defaultdict(float)
        for u in birimler:
            for (rn, kn), o in KOPRU[u].items():
                if (mol == "KIF13A" and kn == n) or (mol == "Rab10" and rn == n):
                    kopru[(rn, kn)] += o / len(birimler)
        kopru_s = ", ".join(f"{ad(diger, rn if mol == 'KIF13A' else kn)} {o:.2f}"
                            for (rn, kn), o in sorted(kopru.items(), key=lambda x: -x[1]) if o >= 0.05)
        satirlar.append({
            "kalinti": ad(mol, n), "bolge": bolge(n),
            **{kol: f"{oz[kol][0]:.2f} ({oz[kol][1]:.2f}-{oz[kol][2]:.2f})" for kol in KOLLAR},
            "kontrol_aralik": f"{k.min():.2f}-{k.max():.2f}", "mutant_aralik": f"{m.min():.2f}-{m.max():.2f}",
            "fark_ort": f"{m.mean() - k.mean():+.2f}", "ayrisma": ayr, "MWU_p": f"{p:.3f}",
            "baslica_ortaklar": ortak, "tuz_koprusu": kopru_s, "_ort": float(np.mean(list(deger.values())))})
    yaz(f"ara_yuz_kalinti_{mol}.tsv", [{k: v for k, v in s.items() if k != "_ort"} for s in satirlar])
    print(f"\n=== {mol} ara yuz kalintilari (12 birimden en az birinde doluluk >= 0.10): {len(satirlar)} ===")
    print(f"{'kalinti':<10} {'bolge':<12} " + " ".join(f"{k[:13]:>17}" for k in KOLLAR) + f" {'fark':>6} {'p':>6}  ayrisma")
    for s in satirlar:
        print(f"{s['kalinti']:<10} {s['bolge']:<12} " + " ".join(f"{s[k]:>17}" for k in KOLLAR)
              + f" {s['fark_ort']:>6} {s['MWU_p']:>6}  {s['ayrisma']}")
    print("  baslica ortaklar:")
    for s in sorted(satirlar, key=lambda x: -x["_ort"]):
        if s["_ort"] >= 0.3:
            print(f"    {s['kalinti']:<10} {s['baslica_ortaklar']}" + (f"  | kopru: {s['tuz_koprusu']}" if s["tuz_koprusu"] else ""))

# 2) cift duzeyi
anahtarlar = sorted({k for u in birimler for k in CIFT[u]})
satirlar = []
for c in anahtarlar:
    deger = {u: CIFT[u].get(c, 0.0) for u in birimler}
    if max(deger.values()) < 0.30:
        continue
    oz = kol_ozet(deger)
    k, m, ayr, p = karsilastir(deger)
    satirlar.append({"Rab10": ad("Rab10", c[0]), "KIF13A": ad("KIF13A", c[1]),
                     "Rab10_bolge": bolge_rab(c[0]), "KIF13A_bolge": bolge_kif(c[1]),
                     **{kol: f"{oz[kol][0]:.2f} ({oz[kol][1]:.2f}-{oz[kol][2]:.2f})" for kol in KOLLAR},
                     "fark_ort": f"{m.mean() - k.mean():+.2f}", "ayrisma": ayr, "MWU_p": f"{p:.3f}"})
yaz("ara_yuz_ciftler.tsv", satirlar)
print(f"\n=== temas ciftleri (12 birimden en az birinde >= 0.30): {len(satirlar)}; ayrisanlar: ===")
for s in satirlar:
    if s["ayrisma"]:
        print(f"  Rab10 {s['Rab10']:<5} - KIF13A {s['KIF13A']:<9} " + " ".join(f"{s[k]:>17}" for k in KOLLAR)
              + f" {s['fark_ort']:>6} p={s['MWU_p']} {s['ayrisma']}")

# 3) 1122 konumunun Rab10 ile temaslari, kopya duzeyinde
print("\n=== 1122 konumunun Rab10 temaslari (kopya duzeyinde, >= 0.02) ===")
for u in birimler:
    for res, cift, kopru, kaynak in kopyalar[u]:
        t = sorted(((rn, o) for (rn, kn), o in cift.items() if kn == 1122 and o >= 0.02), key=lambda x: -x[1])
        b = [(rn, o) for (rn, kn), o in kopru.items() if kn == 1122 and o >= 0.02]
        print(f"  {kaynak:<38} toplam {res.get(('KIF13A', 1122), 0.0):.2f} | "
              + ", ".join(f"{ad('Rab10', rn)} {o:.2f}" for rn, o in t)
              + (f" | kopru: " + ", ".join(f"{ad('Rab10', rn)} {o:.2f}" for rn, o in b) if b else ""))
print("\nbirimler:", ", ".join(f"{u[0]} {u[1]} ({len(kopyalar[u])} kopya)" for u in birimler))
