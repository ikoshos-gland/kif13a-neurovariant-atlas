"""timeseries.tsv'den kosu basina 10 ns blok ortalamalari ve kosu ortalamalari.

E1130-R70 dolulugu: yuklu atomlar arasi en kisa uzaklik <= 0.40 nm (analyze_complex.py tanimi).
Yalniz 10-50 ns penceresi de ayrica verilir (onizlemeler ~48 ns'de bitiyor; adil kiyas icin).
"""
import csv
import sys
from pathlib import Path
from statistics import mean

kok = Path(sys.argv[1])
print("kaynak\tkosu\tns_son\tpencere\tdonus_heliks\tR70E1130_dolu\t1122E1130_nm\tsactokasi_rmsd_nm\tblok_donus_heliks(10ns)")
for f in sorted(kok.glob("*/*/timeseries.tsv")):
    rows = list(csv.DictReader(f.open(), delimiter="\t"))
    t = [float(r["production_ns"]) for r in rows]
    for ad, lo, hi in (("tam", 10.0, 1e9), ("10-50", 10.0, 50.0)):
        sel = [r for r, x in zip(rows, t) if lo <= x <= hi]
        if not sel:
            continue
        th = [float(r["turn_helix_fraction"]) for r in sel]
        sb = [float(r["e1130_r70_nm"]) <= 0.40 for r in sel]
        d = [float(r["res1122_e1130_nm"]) for r in sel]
        rm = [float(r["kif_hairpin_ca_rmsd_nm"]) for r in sel]
        bl = []
        b = lo
        while b < min(hi, max(t)):
            x = [v for v, tt in zip(th, [float(r["production_ns"]) for r in sel]) if b <= tt < b + 10]
            if len(x) > 100:
                bl.append(f"{mean(x):.2f}")
            b += 10
        print(f"{f.parent.parent.name}\t{f.parent.name}\t{max(t):.1f}\t{ad}\t{mean(th):.3f}\t{mean(sb):.3f}\t{mean(d):.2f}\t{mean(rm):.3f}\t{' '.join(bl)}")
