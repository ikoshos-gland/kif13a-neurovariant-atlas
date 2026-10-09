"""Monomer v3 onizleme sonuclarini tek tabloda toplar (10-50 ns penceresi, kopyalar 50 ns'de kirpildi)."""
import csv
import json
import sys
from pathlib import Path

kok = Path(sys.argv[1])
satirlar = []
for d in sorted(p for p in kok.iterdir() if p.is_dir()):
    a = json.loads((d / "analysis.json").read_text())
    ek = json.loads((d / "ek_olcumler.json").read_text()) if (d / "ek_olcumler.json").exists() else {}
    ce = a["comparison_endpoints"]
    alpha = {int(r["mane_position"]): float(r["alpha_torsion_occupancy"])
             for r in csv.DictReader((d / "alpha_residue_occupancies.tsv").open(), delimiter="\t")}
    kopru = {}
    for r in csv.DictReader((d / "contact_occupancies.tsv").open(), delimiter="\t"):
        if r["metric"] == "salt_bridge":
            kopru[f'{r["name_1"]}{r["residue_1"]}-{r["name_2"]}{r["residue_2"]}'] = float(r["occupancy"])
    hk = ek.get("heliks_H", {})
    uz = {u["ad"]: u for u in ek.get("uzakliklar", [])}
    k1122 = next((v for k, v in uz.items() if "1122-" in k and k.endswith("1130")), {})
    s = {
        "gorev": d.name,
        "kare": a["frames"],
        "pencere_ns": "%.2f-%.2f" % tuple(a["analysis_time_ns"]),
        "cekirdek_alfa (birincil)": ce["core_alpha_torsion_fraction"]["mean"],
        "cekirdek_i,i+4_Hbag (birincil)": ce["core_i_i4_hbond_fraction"]["mean"],
        "bolge_alfa_1118-1126": ce["site_alpha_torsion_fraction"]["mean"],
        "1122_temas": ce["site_nonlocal_heavy_contact_count"]["mean"],
        "cekirdek_tuz_koprusu": ce["core_salt_bridge_count"]["mean"],
        "1122_tuz_koprusu": ce["site_salt_bridge_count"]["mean"],
        "cekirdek_CA_RMSD_nm": ce["core_ca_rmsd_nm"]["mean"],
        "oz_goruntu_temas": a["periodic_self_contact_frame_fraction"],
        "alfa1117": alpha.get(1117), "alfa1118": alpha.get(1118),
        "alfa1123": alpha.get(1123), "alfa1124": alpha.get(1124), "alfa1125": alpha.get(1125),
        "DSSP_H1117": hk.get("1117"), "DSSP_H1118": hk.get("1118"),
        "DSSP_H1124": hk.get("1124"), "DSSP_H1125": hk.get("1125"),
        "1122-E1130_nm": k1122.get("ort_nm"), "1122-E1130_dolu": k1122.get("dolu_0.40nm"),
        "K1119-D1126_dolu": uz.get("LYS1119-ASP1126", {}).get("dolu_0.40nm"),
        "K1119-E1122_dolu": uz.get("LYS1119-GLU1122", {}).get("dolu_0.40nm"),
    }
    for k, v in sorted(kopru.items(), key=lambda x: -x[1])[:6]:
        s["kopru " + k] = v
    satirlar.append(s)
alanlar = list(dict.fromkeys(x for s in satirlar for x in s))
w = csv.DictWriter(sys.stdout, fieldnames=alanlar, delimiter="\t")
w.writeheader()
for s in satirlar:
    w.writerow({x: (f"{v:.3f}" if isinstance(v, float) else v) for x, v in s.items()})
