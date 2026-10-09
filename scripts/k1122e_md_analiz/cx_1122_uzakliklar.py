"""Baslangic yapilarinda (complex_heavy.pdb) 1122 konumunun ara yuz kalintilarina ve Rab10'a uzakligi.

Zincir A = Rab10, zincir B = KIF13A 1085-1160 (MANE numaralari). Cikti: CA(1122)-CA(x) uzakliklari (A),
1122 yan zincirinin ve 1124'un Rab10'a en yakin agir atom uzakligi.
"""
import sys
from pathlib import Path

import numpy as np

HEDEF = [1089, 1092, 1106, 1113, 1117, 1118, 1124, 1130, 1134, 1137, 1141, 1144, 1148, 1158]


def oku(p):
    at = []
    for l in open(p):
        if l.startswith(("ATOM", "HETATM")):
            at.append((l[21], int(l[22:26]), l[17:20], l[12:16].strip(),
                       np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return at


print("sistem".ljust(24), " ".join(f"{n:>5}" for n in HEDEF), "| 1122->Rab10  1124->Rab10")
for d in sorted(Path(sys.argv[1]).glob("*/complex_heavy.pdb")):
    at = oku(d)
    ca = {(c, r): x for c, r, n, a, x in at if a == "CA"}
    sc1122 = np.array([x for c, r, n, a, x in at if c == "B" and r == 1122 and a not in ("N", "C", "O", "CA")])
    r1124 = np.array([x for c, r, n, a, x in at if c == "B" and r == 1124])
    rab = [(r, n, x) for c, r, n, a, x in at if c == "A"]
    rx = np.array([x for _, _, x in rab])

    def en_yakin(s):
        dd = np.linalg.norm(s[:, None, :] - rx[None, :, :], axis=2)
        i = np.unravel_index(dd.argmin(), dd.shape)
        return dd.min(), f"{rab[i[1]][1]}{rab[i[1]][0]}"

    a, b = en_yakin(sc1122)
    c, e = en_yakin(r1124)
    print(f"{d.parent.name:<24}", " ".join(f"{np.linalg.norm(ca[('B', 1122)] - ca[('B', n)]):>5.1f}" for n in HEDEF),
          f"| {a:4.1f} A {b:<7} {c:4.1f} A {e}")
