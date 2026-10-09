"""Kompleks MD (Rab10-GTP-Mg + KIF13A 1085-1160) uc noktalarini tum kosulardan tek tabloda toplar.

Girdi: <kok>/<kaynak>/<kosu>/{analysis.json, kif13a_dssp.tsv, intermolecular_salt_bridges.tsv}
kaynak: truba_tam (bitmis 100 ns), truba_oniz (yarim kosu kopyasi), rig (rig'de kosulmus).
"""
import csv
import json
import sys
from pathlib import Path

kok = Path(sys.argv[1])
satirlar = []
for kaynak_dir in sorted(p for p in kok.iterdir() if p.is_dir()):
    for kosu in sorted(p for p in kaynak_dir.iterdir() if p.is_dir()):
        aj = kosu / "analysis.json"
        if not aj.exists():
            continue
        a = json.loads(aj.read_text())
        ep = a.get("endpoints", {})
        dssp = {}
        f = kosu / "kif13a_dssp.tsv"
        if f.exists():
            for r in csv.DictReader(f.open(), delimiter="\t"):
                dssp[int(r["kif13a"])] = float(r["helix_H"])
        sb = {}
        f = kosu / "intermolecular_salt_bridges.tsv"
        if f.exists():
            for r in csv.DictReader(f.open(), delimiter="\t"):
                sb[f'{r["receptor"].replace("Rab10 ", "")}-{r["kif13a"].replace("KIF13A ", "")}'] = float(r["occupancy"])
        ad = kosu.name
        kol = ad.rsplit("_s", 1)[0].replace("_wt", "").replace("_k1122e", "")
        k = {
            "kaynak": kaynak_dir.name,
            "kosu": ad,
            "kol": kol,
            "kare": a.get("frames"),
            "R70-E1130": ep.get("e1130_r70_sb_occupancy"),
            "1122-E1130 kopru": ep.get("k1122_e1130_sb_occupancy"),
            "1122-E1130 temas": ep.get("res1122_e1130_contact_occupancy"),
            "1122-E1130 nm": ep.get("res1122_e1130_dist_nm_mean"),
            "D45-R1144": sb.get("D45-R1144"),
            "D45-R1106": sb.get("D45-R1106"),
            "temas cifti": ep.get("interface_contact_pairs_mean"),
            "gomulu SASA nm2": ep.get("buried_sasa_nm2_mean"),
            "donus heliks": ep.get("turn_helix_fraction_1116_1126"),
            "H1117": dssp.get(1117),
            "H1118": dssp.get(1118),
            "H1124": dssp.get(1124),
            "H1125": dssp.get(1125),
            "sac tokasi RMSD nm": ep.get("kif_hairpin_ca_rmsd_nm_mean"),
            "heliksler arasi aci": ep.get("interhelix_angle_deg_mean"),
            "Na+ 1122": ep.get("na_occupancy_1122"),
            "molekul arasi kopru": ep.get("intermolecular_sb_count_mean"),
        }
        for g in ("gbsa", "mmgbsa", "gbsa_relative"):
            if g in a:
                k["gbsa"] = json.dumps(a[g])[:120]
        satirlar.append(k)

alanlar = list(dict.fromkeys(x for s in satirlar for x in s))
w = csv.DictWriter(sys.stdout, fieldnames=alanlar, delimiter="\t")
w.writeheader()
for s in satirlar:
    w.writerow({x: (f"{v:.3f}" if isinstance(v, float) else v) for x, v in s.items()})
