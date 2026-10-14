"""Rab10 catalytic Q68 and the gamma-phosphate water in the finished complex runs (10-100 ns, 2026-10-14).

    python q68_analiz.py run --run-dir ../runs/WT_wt_s2 --system-dir ../systems/WT_wt_s2 --out ../results_q68/WT_wt_s2.json
    python q68_analiz.py ozet --in-dir ../results_q68 --out ../results_q68/ozet

Question: does K1122E (or K1122R) change the pre-catalytic geometry of Rab10 switch II, i.e. where Q68 sits relative
to the GTP gamma-phosphate and whether a water is placed for in-line attack? Classical MD cannot hydrolyse GTP;
these are geometric proxies only.

Solute trajectory (10 ps): Q68 amide (OE1, NE2) to PG and to the gamma oxygens; Q68 CA to PG; R70-E1130 and
1122-E1130 for context. Full-system trajectory (100 ps, runs_solvent): water oxygens near PG, excluding the Mg
shell (O-Mg < 2.6 A):
  attack water   O-PG < 4.0 A and angle O3B-PG-O >= 150 deg (in line with the leaving beta-gamma oxygen)
  bridged water  attack-region water (O-PG < 4.5 A) that is also within 3.4 A of Q68 OE1/NE2
Summary: one value per run (mean over frames); arms pooled as controls (WT + K1122R) vs mutants (K1122E paired +
K1122E AF3), 6 vs 6, exact two-sided permutation p on the difference of means (minimum attainable 2/924 = 0.0022).
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def min_image(vec, box):
    frac = vec @ np.linalg.inv(box)
    return (frac - np.round(frac)) @ box


def run(args):
    import mdtraj as md
    from analyze_complex import Complex
    from complex_common import load_protocol
    protocol = load_protocol()
    cx = Complex(args.run_dir, args.system_dir, protocol, allow_smoke=False)
    t = cx.traj
    q = cx.res("Rab10", 68)
    qa = {a.name: a.index for a in q.atoms}
    gtp = cx.res("GTP", 0)
    ga = {a.name: a.index for a in gtp.atoms}
    amide = [qa["OE1"], qa["NE2"]]
    gam_o = [ga["O1G"], ga["O2G"], ga["O3G"]]
    dist = lambda A, B: md.compute_distances(t, np.array([(i, j) for i in A for j in B])).min(axis=1) * 10
    chg = lambda mol, num, names: [a.index for a in cx.res(mol, num).atoms if a.name in names]
    sol = {
        "q68_amid_PG": dist(amide, [ga["PG"]]),
        "q68_amid_Ogama": dist(amide, gam_o),
        "q68_CA_PG": dist([qa["CA"]], [ga["PG"]]),
        "r70_e1130": dist(chg("Rab10", 70, {"NE", "NH1", "NH2"}), chg("KIF13A", 1130, {"OE1", "OE2"})),
    }
    out = {"system_id": args.system_dir.name, "solute_frames": len(t),
           "solute": {k: {"mean": float(v.mean()), "p10": float(np.percentile(v, 10)), "p50": float(np.median(v)),
                          "p90": float(np.percentile(v, 90))} for k, v in sol.items()},
           "q68_amid_PG_lt5": float((sol["q68_amid_PG"] < 5.0).mean()),
           "q68_amid_Ogama_lt35": float((sol["q68_amid_Ogama"] < 3.5).mean())}
    # full system, 100 ps
    run_json = json.loads((args.run_dir / "run.json").read_text())
    sdir = args.run_dir.parent.parent / f"{args.run_dir.parent.name}_solvent" / args.run_dir.name
    sdcd = sdir / "system_solvent.dcd"
    full = md.load_prmtop(str(args.system_dir / "complex.prmtop"))
    # indices in the full topology
    rmap = json.loads((args.system_dir / "residue_map.json").read_text())["residues"]
    ridx = lambda mol, num: [r["prmtop_residue_index"] for r in rmap if r["molecule"] == mol and r["uniprot_number"] == num][0]
    gres = [r["prmtop_residue_index"] for r in rmap if r["molecule"] == "GTP"][0]
    mres = [r["prmtop_residue_index"] for r in rmap if r["molecule"] == "MG"][0]
    fq = {a.name: a.index for a in full.residue(ridx("Rab10", 68)).atoms}
    fg = {a.name: a.index for a in full.residue(gres).atoms}
    fmg = [a.index for a in full.residue(mres).atoms][0]
    wo = np.array([a.index for a in full.atoms if a.residue.name in {"WAT", "HOH"} and a.name == "O"])
    sel = np.array([fg["PG"], fg["O3B"], fmg, fq["OE1"], fq["NE2"]] + wo.tolist())
    siv = run_json["solvent_report_interval_steps"]
    start, dt = run_json["production_start_step"], run_json["timestep_fs"]
    first = (start // siv + 1) * siv
    n_att, n_br, frames = [], [], 0
    k = 0
    for chunk in md.iterload(str(sdcd), top=full, chunk=50, atom_indices=sel):
        tns = (first + (k + np.arange(len(chunk))) * siv - start) * dt * 1e-6
        k += len(chunk)
        for f in np.where(tns > cx.discard_ns)[0]:
            x = chunk.xyz[f] * 10
            box = chunk.unitcell_vectors[f] * 10
            pg, o3b, mg, oe1, ne2 = x[0], x[1], x[2], x[3], x[4]
            w = pg + min_image(x[5:] - pg, box)
            dpg = np.linalg.norm(w - pg, axis=1)
            near = dpg < 4.5
            dmg = np.linalg.norm(w[near] - mg, axis=1)
            wn, dn = w[near][dmg >= 2.6], dpg[near][dmg >= 2.6]
            v1, v2 = o3b - pg, wn - pg
            ang = np.degrees(np.arccos(np.clip((v2 @ v1) / (np.linalg.norm(v2, axis=1) * np.linalg.norm(v1) + 1e-9), -1, 1)))
            att = (dn < 4.0) & (ang >= 150)
            dq = np.minimum(np.linalg.norm(wn - oe1, axis=1), np.linalg.norm(wn - ne2, axis=1))
            n_att.append(int(att.sum()))
            n_br.append(int((dq < 3.4).sum()))
            frames += 1
    n_att, n_br = np.array(n_att), np.array(n_br)
    out.update({"solvent_frames": frames, "attack_water_frac": float((n_att > 0).mean()),
                "attack_water_mean": float(n_att.mean()), "bridged_water_frac": float((n_br > 0).mean()),
                "bridged_water_mean": float(n_br.mean())})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1))
    print(json.dumps({"system": out["system_id"], "q68_PG_mean": round(out["solute"]["q68_amid_PG"]["mean"], 2),
                      "attack": round(out["attack_water_frac"], 3), "bridged": round(out["bridged_water_frac"], 3)}))


def perm_p(a, b):
    a, b = np.asarray(a), np.asarray(b)
    pool, n, obs = np.concatenate([a, b]), len(a), abs(a.mean() - b.mean())
    hits = tot = 0
    for idx in itertools.combinations(range(len(pool)), n):
        m = np.zeros(len(pool), bool); m[list(idx)] = True
        hits += abs(pool[m].mean() - pool[~m].mean()) >= obs - 1e-12
        tot += 1
    return hits / tot


def ozet(args):
    rows = [json.loads(p.read_text()) for p in sorted(args.in_dir.glob("*.json"))]
    arm = lambda s: "K1122E_AF3" if s.startswith("K1122E_AF3") else s.split("_wt_")[0].replace("_paired", "")
    keys = [("q68_amid_PG", lambda r: r["solute"]["q68_amid_PG"]["mean"]),
            ("q68_amid_Ogama", lambda r: r["solute"]["q68_amid_Ogama"]["mean"]),
            ("q68_CA_PG", lambda r: r["solute"]["q68_CA_PG"]["mean"]),
            ("q68_amid_PG_lt5", lambda r: r["q68_amid_PG_lt5"]),
            ("attack_water_frac", lambda r: r["attack_water_frac"]),
            ("bridged_water_frac", lambda r: r["bridged_water_frac"]),
            ("r70_e1130", lambda r: r["solute"]["r70_e1130"]["mean"])]
    lines = ["system\tarm\t" + "\t".join(k for k, _ in keys)]
    for r in rows:
        lines.append(r["system_id"] + "\t" + arm(r["system_id"]) + "\t" + "\t".join(f"{f(r):.4f}" for _, f in keys))
    ctrl = [r for r in rows if arm(r["system_id"]) in {"WT", "K1122R"}]
    mut = [r for r in rows if arm(r["system_id"]) in {"K1122E", "K1122E_AF3"}]
    summ = {}
    for k, f in keys:
        a, b = [f(r) for r in ctrl], [f(r) for r in mut]
        summ[k] = {"controls": [round(v, 4) for v in a], "mutants": [round(v, 4) for v in b],
                   "controls_mean": float(np.mean(a)), "mutants_mean": float(np.mean(b)),
                   "p_perm": perm_p(a, b) if len(a) == 6 and len(b) == 6 else None}
    Path(str(args.out) + ".tsv").write_text("\n".join(lines) + "\n")
    Path(str(args.out) + ".json").write_text(json.dumps({"n_controls": len(ctrl), "n_mutants": len(mut), "endpoints": summ}, indent=1))
    for k, v in summ.items():
        print(f"{k:22s} kontrol {v['controls_mean']:.3f}  mutant {v['mutants_mean']:.3f}  p={v['p_perm']}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run"); r.add_argument("--run-dir", type=Path, required=True)
    r.add_argument("--system-dir", type=Path, required=True); r.add_argument("--out", type=Path, required=True)
    o = sub.add_parser("ozet"); o.add_argument("--in-dir", type=Path, required=True); o.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    run(a) if a.cmd == "run" else ozet(a)


if __name__ == "__main__":
    main()
