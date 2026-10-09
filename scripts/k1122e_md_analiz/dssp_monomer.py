"""Monomer onizleme ek olcumleri: kalinti bazinda DSSP heliks orani (kompleks analiziyle ayni tanim,
mdtraj simplified DSSP 'H') ve 1122 cevresindeki yuklu grup uzakliklari. Ilk 10 ns atilir.
Kullanim: python dssp_monomer.py <gorev_dizini>   ->  <gorev_dizini>/ek_olcumler.json
"""
import json
import sys
from pathlib import Path

import mdtraj as md
import numpy as np

d = Path(sys.argv[1])
build = json.loads((d / "build.json").read_text())
run = json.loads((d / "run.json").read_text())
traj = md.load_dcd(str(d / "trajectory.dcd"), top=str(d / "protein.pdb"))
start, iv, dt = run["production_start_step"], run["report_interval_steps"], run["timestep_fs"]
first = ((start // iv) + 1) * iv
t = (first + np.arange(len(traj)) * iv - start) * dt / 1e6
traj = traj[t > 10.0]
res = {}
for entry, r in zip(build["protein_residue_map"], traj.topology.residues):
    if entry["mane_position"] is not None:
        res[entry["mane_position"]] = r
dssp = md.compute_dssp(traj, simplified=True)
idx = {r.index: j for j, r in enumerate(traj.topology.residues)}
heliks = {p: float((dssp[:, idx[r.index]] == "H").mean()) for p, r in sorted(res.items()) if 1100 <= p <= 1140}
YUK = {"LYS": ["NZ"], "ARG": ["NE", "NH1", "NH2"], "GLU": ["OE1", "OE2"], "ASP": ["OD1", "OD2"]}


def grup(p):
    r = res[p]
    return [a.index for a in r.atoms if a.name in YUK.get(r.name, [])]


def uzaklik(p, q):
    a, b = grup(p), grup(q)
    pairs = np.array([[i, j] for i in a for j in b])
    dmin = md.compute_distances(traj, pairs).min(axis=1)
    return {"ad": f"{res[p].name}{p}-{res[q].name}{q}", "ort_nm": float(dmin.mean()),
            "dolu_0.40nm": float((dmin <= 0.40).mean())}


ciftler = [(1122, 1130), (1119, 1126), (1119, 1122), (1122, 1125), (1122, 1126), (1129, 1122)]
out = {"gorev": build["task"]["task_id"], "kare": int(len(traj)), "pencere_ns": [10.0, float(t.max())],
       "heliks_H": heliks, "uzakliklar": [uzaklik(p, q) for p, q in ciftler]}
(d / "ek_olcumler.json").write_text(json.dumps(out, indent=1))
print(json.dumps({"gorev": out["gorev"], "kare": out["kare"], "H1117": heliks.get(1117), "H1124": heliks.get(1124)}))
