"""Frames + per-frame measurements for the 100 ns comparison viewer (WT vs K1122E from the same start).

    python export_viewer2.py --run WT_wt_s2=runs_rig/WT_wt_s2 --run K1122E_paired_wt_s2=runs_rig/K1122E_paired_wt_s2 \
        --stride3d 30 --stride-chart 10 --out viewer2_frames.json

3D frames: backbone N/CA/C/O of Rab10 and KIF13A, heavy side chains of the residues in SIDE, GTP, Mg and the
Na+ ions that come closest to the KIF13A acidic turn; coordinates in 0.1 A integers, all systems superposed on
the Rab10 CA of the first system's first frame. Distances are computed exactly from the full trajectory (charged
atoms, closest pair, periodic minimum image for Na+) and stored per frame, so labels do not inherit the 0.1 A
rounding. Charts: every --stride-chart frames, plus DSSP helix flags for KIF13A 1114-1128.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import mdtraj as md

from analyze_complex import Complex
from complex_common import load_protocol

HERE = Path(__file__).resolve().parent
CHAIN = {"Rab10": "A", "KIF13A": "B", "GTP": "C", "MG": "D"}
RENAME = {"HID": "HIS", "HIE": "HIS", "HIP": "HIS", "CYX": "CYS", "gtp": "GTP"}
SIDE = {"Rab10": {45, 70, 71}, "KIF13A": {1106, 1117, 1118, 1121, 1122, 1124, 1125, 1126, 1127, 1130, 1144}}
CHG = {"LYS": ["NZ"], "ARG": ["NE", "NH1", "NH2"], "GLU": ["OE1", "OE2"], "ASP": ["OD1", "OD2"]}
BACKBONE = {"N", "CA", "C", "O"}
PAIRS = [("e1130_r70", ("KIF13A", 1130), ("Rab10", 70)),
         ("r1122_e1130", ("KIF13A", 1122), ("KIF13A", 1130)),
         ("r1122_r70", ("KIF13A", 1122), ("Rab10", 70)),
         ("d45_r1144", ("Rab10", 45), ("KIF13A", 1144))]
HELIX_RES = list(range(1114, 1129))
N_NA = 4


def pdb_line(serial, name, resname, chain, resnum, xyz, element, het):
    nm = name if len(name) == 4 else " " + name.ljust(3)
    rec = "HETATM" if het else "ATOM  "
    return (f"{rec}{serial:5d} {nm} {resname:>3s} {chain}{resnum:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2s}")


def charged(cx, mol, num):
    res = cx.res(mol, num)
    return [a.index for a in res.atoms if a.name in CHG.get(cx.chem(res), [])]


def min_image(vec, box):
    frac = vec @ np.linalg.inv(box)
    return (frac - np.round(frac)) @ box


def export(system_id, run_dir, protocol, stride3d, stride_chart, ref_xyz=None):
    cx = Complex(run_dir, HERE.parent / "systems" / system_id, protocol, allow_smoke=False)
    traj, top = cx.traj, cx.traj.topology
    n = len(traj)
    rab_ca = [a.index for a in top.atoms if a.name == "CA"
              and (cx.rows_by_index.get(a.residue.index) or {}).get("molecule") == "Rab10"]
    # exact distances, every frame (nm -> A)
    dist = {}
    for key, ga, gb in PAIRS:
        a, b = charged(cx, *ga), charged(cx, *gb)
        pairs = np.array([(i, j) for i in a for j in b])
        dist[key] = md.compute_distances(traj, pairs, periodic=False).min(axis=1) * 10
    # Na+: closest ion to residue 1122 / 1130 charged atoms (periodic), and the ions to draw
    na = cx.ions["Na+"]
    g1122, g1130 = charged(cx, "KIF13A", 1122), charged(cx, "KIF13A", 1130)
    d_na = {}
    for key, g in (("na_1122", g1122), ("na_1130", g1130)):
        pairs = np.array([(i, j) for i in g for j in na])
        d = md.compute_distances(traj, pairs, periodic=True).reshape(n, len(g), len(na)).min(axis=1) * 10
        d_na[key] = d
    near = (np.minimum(d_na["na_1122"], d_na["na_1130"]) < 6.0).mean(axis=0)
    chosen = [na[k] for k in np.argsort(-near)[:N_NA] if near[k] > 0]
    dist["na_1122"] = d_na["na_1122"].min(axis=1)
    dist["na_count_1122"] = (d_na["na_1122"] < 3.2).sum(axis=1).astype(float)
    # DSSP helix flags
    kif_res = {cx.rows_by_index[r.index]["uniprot_number"]: r.index for r in top.residues
               if r.index in cx.rows_by_index and cx.rows_by_index[r.index]["molecule"] == "KIF13A"}
    ch_idx = np.arange(0, n, stride_chart)
    dssp = md.compute_dssp(traj[ch_idx], simplified=True)
    helix = {u: (dssp[:, kif_res[u]] == "H").astype(int) for u in HELIX_RES}
    # atoms for 3D
    keep, meta = [], []
    for atom in top.atoms:
        row = cx.rows_by_index.get(atom.residue.index)
        if row is None or row["cap"] or row["molecule"] not in CHAIN or atom.element.symbol == "H":
            continue
        mol = row["molecule"]
        if mol in SIDE and atom.name not in BACKBONE and row["uniprot_number"] not in SIDE[mol]:
            continue
        keep.append(atom.index)
        num = row["uniprot_number"] if row["uniprot_number"] is not None else 1
        meta.append((atom.name, RENAME.get(row["resname_prmtop"], row["resname_prmtop"]), CHAIN[mol], num,
                     atom.element.symbol, mol in ("GTP", "MG")))
    for k, i in enumerate(chosen):
        keep.append(i)
        meta.append(("NA", "NA", "E", k + 1, "NA", True))
    idx3 = np.arange(0, n, stride3d)
    sub = traj[idx3]
    xyz = sub.xyz.copy()
    # Na+ to the periodic copy nearest the acidic turn
    anchor = xyz[:, g1122 + g1130].mean(axis=1)
    for i in chosen:
        for f in range(len(sub)):
            xyz[f, i] = anchor[f] + min_image(xyz[f, i] - anchor[f], sub.unitcell_vectors[f])
    sub.xyz = xyz
    if ref_xyz is None:
        ref_xyz = sub.xyz[0, rab_ca].copy()
    # every frame onto the reference Rab10 CA (same atom order in all systems: identical Rab10 topology)
    full_ref = sub[0].xyz.copy()
    full_ref[0, rab_ca] = ref_xyz
    sub.superpose(md.Trajectory(full_ref, sub.topology), 0, atom_indices=rab_ca)
    out_xyz = sub.xyz[:, keep, :] * 10.0
    center = ref_xyz.mean(axis=0) * 10.0
    out_xyz -= center
    lines = [pdb_line(i + 1, nm, r, c, num, out_xyz[0, i], e, h) for i, (nm, r, c, num, e, h) in enumerate(meta)]
    times = cx.times
    r = lambda a, d=2: [round(float(v), d) for v in a]
    return ref_xyz, {
        "system_id": system_id, "n_atoms": len(keep), "n_frames": len(sub), "n_na": len(chosen),
        "times_ns": r(times[idx3], 2),
        "pdb": "\n".join(lines) + "\nEND\n",
        "coords": np.round(out_xyz.reshape(len(sub), -1) * 10).astype(int).tolist(),
        "frame_dist": {k: r(v[idx3], 2) for k, v in dist.items()},
        "chart": {"t_ns": r(times[ch_idx], 2), **{k: r(v[ch_idx], 2) for k, v in dist.items()},
                  "helix": {str(u): helix[u].tolist() for u in HELIX_RES}},
        "frame_helix": {str(u): [int(helix[u][min(len(ch_idx) - 1, round(i / stride_chart))]) for i in idx3]
                        for u in HELIX_RES},
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="append", required=True, help="SYSTEM_ID=RUN_DIR")
    ap.add_argument("--stride3d", type=int, default=30)
    ap.add_argument("--stride-chart", type=int, default=10)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    protocol = load_protocol()
    out, ref = [], None
    for s in a.run:
        sid, rd = s.split("=", 1)
        ref, d = export(sid, Path(rd), protocol, a.stride3d, a.stride_chart, ref)
        out.append(d)
    a.out.write_text(json.dumps({"systems": out}, separators=(",", ":")))
    print({o["system_id"]: (o["n_atoms"], o["n_frames"], o["n_na"]) for o in out}, a.out.stat().st_size // 1024, "kB")


if __name__ == "__main__":
    main()
