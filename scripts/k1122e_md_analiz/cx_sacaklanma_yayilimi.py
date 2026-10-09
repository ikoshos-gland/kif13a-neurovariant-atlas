"""Saçaklanma ara yüze yansıyor mu? Mevcut zaman serilerinde kosullu analiz (yeni hesap gerekmez).

Soru: bir kosunun icinde donus saçaklandigi karelerde ara yuz olculeri degisiyor mu?
Kareler 10 ps arayla ve guclu otokorelasyonlu, o yuzden iki duzeyde bakiyoruz:
  (a) kare duzeyi Spearman: yalniz betimleyici, p degeri verilmiyor
  (b) 10 ns blok ortalamalari, kosu ici ortalamasi cikarilmis (sabit etki), kol icinde havuzlanmis
"""
import csv
import sys
from pathlib import Path
from statistics import mean

import numpy as np
from scipy import stats

kok = Path(sys.argv[1])
ARA_YUZ = ["e1130_r70_nm", "intermolecular_sb_count", "kif_helix2_rmsd_nm"]
kollar = {}
print("kosu                       n     donus_heliks   " + "  ".join(f"rho({a.split('_')[0]})" for a in ARA_YUZ))
for f in sorted(kok.glob("*/*/timeseries.tsv")):
    rows = list(csv.DictReader(f.open(), delimiter="\t"))
    t = np.array([float(r["production_ns"]) for r in rows])
    m = t > 10.0
    t = t[m]
    th = np.array([float(r["turn_helix_fraction"]) for r in rows])[m]
    ad = f.parent.name
    kol = ad.rsplit("_s", 1)[0].replace("_wt", "").replace("_k1122e", "")
    rhos = []
    for a in ARA_YUZ:
        v = np.array([float(r[a]) for r in rows])[m]
        rhos.append(stats.spearmanr(th, v).statistic)
        # 10 ns bloklari, kosu ici merkezlenmis
        bi = np.floor((t - 10.0) / 10.0).astype(int)
        bth, bv = [], []
        for b in np.unique(bi):
            s = bi == b
            if s.sum() > 500:
                bth.append(th[s].mean())
                bv.append(v[s].mean())
        if len(bth) >= 3:
            bth = np.array(bth) - np.mean(bth)
            bv = np.array(bv) - np.mean(bv)
            kollar.setdefault((kol, a), []).extend(zip(bth, bv))
    print(f"{ad:<26} {len(t):>5}  {th.mean():>12.3f}   " + "  ".join(f"{r:>+10.2f}" for r in rhos))

print("\n10 ns bloklari, kosu ici merkezlenmis, kol icinde havuzlanmis (Pearson):")
print(f"{'kol':<16} {'olcum':<26} {'n_blok':>6} {'r':>7} {'p':>9}")
for (kol, a), xy in sorted(kollar.items()):
    x = np.array([p[0] for p in xy])
    y = np.array([p[1] for p in xy])
    if len(x) < 6:
        continue
    r = stats.pearsonr(x, y)
    print(f"{kol:<16} {a:<26} {len(x):>6} {r.statistic:>+7.2f} {r.pvalue:>9.3f}")
