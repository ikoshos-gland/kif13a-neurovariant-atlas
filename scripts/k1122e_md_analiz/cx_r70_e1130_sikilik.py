"""Rig ciftinden cikan 'mutant R70-E1130'u sikilastiriyor' iddiasinin olcumlerini butun kosularda hesaplar.

Olcumler (10 ns sonrasi, 10 ps kareler): R70-E1130 dolulugu (<=0.40 nm), uzaklik dagiliminin 90. yuzdeligi,
en uzun kesintisiz kapali sure (ns), acilma sayisi; interface_contacts.tsv'den F71-E1130 ve R70-V1127 temaslari.
"""
import csv
import sys
from pathlib import Path

import numpy as np

kok = Path(sys.argv[1])
print("kaynak\tkosu\tns\tR70E1130_dolu\tp90_nm\ten_uzun_kapali_ns\tacilma_sayisi\tF71-E1130\tR70-V1127")
for f in sorted(kok.glob("*/*/timeseries.tsv")):
    rows = list(csv.DictReader(f.open(), delimiter="\t"))
    t = np.array([float(r["production_ns"]) for r in rows])
    d = np.array([float(r["e1130_r70_nm"]) for r in rows])
    m = t > 10.0
    t, d = t[m], d[m]
    kapali = d <= 0.40
    en_uzun, cur = 0, 0
    for k in kapali:
        cur = cur + 1 if k else 0
        en_uzun = max(en_uzun, cur)
    acilma = int(np.sum(kapali[:-1] & ~kapali[1:]))
    temas = {}
    ic = f.parent / "interface_contacts.tsv"
    if ic.exists():
        for r in csv.DictReader(ic.open(), delimiter="\t"):
            temas[(r["receptor"], r["kif13a"])] = float(r["occupancy"])
    f71 = temas.get(("Rab10 F71", "KIF13A E1130"), 0.0)
    v1127 = temas.get(("Rab10 R70", "KIF13A V1127"), 0.0)
    print(f"{f.parent.parent.name}\t{f.parent.name}\t{t.max():.0f}\t{kapali.mean():.3f}\t{np.percentile(d, 90):.3f}"
          f"\t{en_uzun * 0.01:.1f}\t{acilma}\t{f71:.2f}\t{v1127:.2f}")
