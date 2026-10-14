"""Monomer MD v3 (30 tasks): one row per run over the full analysis window, plus K1122E vs controls (WT + K1122R).

    python monomer_v3_tablo.py --runs k1122e_md_production_v3 --out monomer_v3_tam_pencere_20261014.tsv

Inputs per task: analysis.json, metrics_timeseries.tsv (10-100 ns primary 1092-1146, 10-50 ns padded 1084-1154),
alpha_residue_occupancies.tsv. Exact two-sided permutation p on the difference of means, 3 vs 6 (minimum 1/84).
"""
import argparse, csv, glob, itertools, json, os
import numpy as np

ARMS = ["WT", "K1122R", "K1122E", "K1104E", "E1128K"]
KEYS = ["cekirdek_alfa", "cekirdek_hbag", "bolge_1118_1126_alfa", "alfa_1117", "alfa_1120", "alfa_1124", "alfa_1125",
        "cekirdek_ca_rmsd_nm"]


def perm(a, b):
    p, n = np.array(a + b), len(a)
    obs, hits, tot = abs(np.mean(a) - np.mean(b)), 0, 0
    for idx in itertools.combinations(range(len(p)), n):
        m = np.zeros(len(p), bool); m[list(idx)] = True
        hits += abs(p[m].mean() - p[~m].mean()) >= obs - 1e-12; tot += 1
    return hits / tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = []
    for d in sorted(glob.glob(os.path.join(a.runs, "*_r0*"))):
        part, arm, rep = os.path.basename(d).split("_")
        ts = list(csv.DictReader(open(f"{d}/metrics_timeseries.tsv"), delimiter="\t"))
        t = np.array([float(r["production_time_ns"]) for r in ts])
        m = lambda k: float(np.mean([float(r[k]) for r in ts]))
        al = {int(r["mane_position"]): float(r["alpha_torsion_occupancy"])
              for r in csv.DictReader(open(f"{d}/alpha_residue_occupancies.tsv"), delimiter="\t")}
        rows.append(dict(task=os.path.basename(d), parca=part, kol=arm, tekrar=rep, pencere_ns=f"{t.min():.2f}-{t.max():.2f}",
                         kare=len(t), cekirdek_alfa=m("core_alpha_torsion_fraction"), cekirdek_hbag=m("core_i_i4_hbond_fraction"),
                         bolge_1118_1126_alfa=m("site_alpha_torsion_fraction"), cekirdek_ca_rmsd_nm=m("core_ca_rmsd_nm"),
                         alfa_1117=al[1117], alfa_1120=al[1120], alfa_1124=al[1124], alfa_1125=al[1125],
                         oz_temas=m("periodic_self_contact_pair_count")))
    cols = list(rows[0])
    with open(a.out, "w", newline="\n") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(f"{r[k]:.4f}" if isinstance(r[k], float) else str(r[k]) for k in cols) + "\n")
    for part in ("primary", "padded"):
        print("==", part)
        for k in KEYS:
            g = {arm: [r[k] for r in rows if r["parca"] == part and r["kol"] == arm] for arm in ARMS}
            print(f" {k:22s}", "  ".join(f"{arm} {min(v):.3f}-{max(v):.3f}" for arm, v in g.items()),
                  f"| p(K1122E vs WT+K1122R)={perm(g['K1122E'], g['WT'] + g['K1122R']):.3f}")


if __name__ == "__main__":
    main()
