"""Endpoints of the Rab10(GTP, Mg):KIF13A complex runs, keyed by UniProt/MANE numbers, and replica aggregation.

    python analyze_complex.py run --run-dir runs/WT_wt_s5 --system-dir systems/WT_wt_s5 [--gbsa] [--solvent]
    python analyze_complex.py aggregate --runs-root runs --out-dir results
    python analyze_complex.py selftest --out-dir /tmp/cx_selftest --system-dir systems/WT_wt_s5

run: one trajectory (trajectory.dcd = solute + 2 Mg-shell waters + ions, 10 ps). Production frames after the
burn-in (protocol analysis.discard_production_ns) are used; smoke runs only with --allow-smoke (no burn-in,
diagnostic only). KIF13A, GTP, Mg and the shell waters are re-imaged next to Rab10 in every frame (whole
complex, triclinic minimum image), so intra-complex geometry is computed without periodic artefacts;
ion distances use the periodic minimum image.

Endpoints (complex_protocol.json analysis block): E1130-R70 salt bridge (charged atoms <= 0.40 nm) occupancy
and distance distribution; residue 1122 charged atoms to R70 and to E1130; K/R1122-E1130 occupancy; every
Rab10(+GTP):KIF13A salt bridge with occupancy; interface residue contacts (closest heavy atoms <= 0.45 nm);
buried SASA; KIF13A turn 1116-1126 DSSP and helix1/helix2 angle; KIF13A helix-2 and Rab10 switch I/II RMSD
after superposition on the Rab10 core; Na+ within 0.32 nm of E1122/E1124/D1125/D1126/E1130 carboxylate O;
Rab10 H72 environment; Mg coordination; periodic image distance; optional GBn2 interaction energy every
1 ns (relative only: no entropy, single trajectory, ff19SB not trained for GB).

aggregate: per arm per-replica values, mean, t-based 95% CI; Welch differences vs WT (and the other
protocol comparisons) with exact permutation p and the minimum p the design can attain; paired-by-start
differences with exact sign-flip p for the paired arms.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from complex_common import DEFAULT_PROTOCOL, load_protocol, provenance, read_json, sha256, utc_now, write_json

POS = {"ARG": ["NE", "NH1", "NH2"], "LYS": ["NZ"], "HIP": ["ND1", "NE2"]}
NEG = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}
GTP_NEG = ["O1G", "O2G", "O3G", "O1B", "O2B", "O1A", "O2A"]
ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIE": "H",
       "HID": "H", "HIP": "H", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
       "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("empty\n")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        fields = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", restval="")
        writer.writeheader()
        writer.writerows(rows)


def autocorrelation_quality(values) -> dict:
    values = np.asarray(values, dtype=float)
    if len(values) < 4 or float(np.var(values)) == 0:
        return {"integrated_autocorrelation_frames": 1.0, "effective_frames": float(len(values))}
    centered = values - values.mean()
    size = 1 << (2 * len(values) - 1).bit_length()
    spectrum = np.fft.rfft(centered, n=size)
    cov = np.fft.irfft(spectrum * np.conjugate(spectrum), n=size)[:len(values)] / np.arange(len(values), 0, -1)
    corr = cov / cov[0]
    paired = corr[1:len(values) - (len(values) % 2)]
    paired = paired[:len(paired) // 2 * 2].reshape(-1, 2).sum(axis=1)
    stop = next((i for i, item in enumerate(paired) if item <= 0), len(paired))
    tau = max(1.0, float(1 + 2 * paired[:stop].sum()))
    return {"integrated_autocorrelation_frames": tau, "effective_frames": float(len(values) / tau)}


def describe(values, cutoff: float | None = None) -> dict:
    values = np.asarray(values, dtype=float)
    out = {"mean": float(values.mean()), "sd": float(values.std()), "median": float(np.median(values)),
           "q05": float(np.quantile(values, 0.05)), "q95": float(np.quantile(values, 0.95)),
           "first_vs_second_half": float(values[len(values) // 2:].mean() - values[:max(1, len(values) // 2)].mean()),
           **autocorrelation_quality(values)}
    if cutoff is not None:
        out["occupancy"] = float((values <= cutoff).mean())
    return out


# ----------------------------------------------------------------------------- loading and imaging

class Complex:
    """Trajectory + topology + UniProt-keyed residue lookup for one run."""

    def __init__(self, run_dir: Path, system_dir: Path, protocol: dict, allow_smoke: bool, discard_ns: float | None = None):
        import mdtraj as md
        self.run_dir, self.system_dir, self.protocol = run_dir, system_dir, protocol
        self.run = read_json(run_dir / "run.json")
        if not self.run.get("completed"):
            raise ValueError(f"{run_dir.name}: run not complete")
        if self.run["smoke"] and not allow_smoke:
            raise ValueError("Smoke trajectory excluded from scientific analysis (use --allow-smoke for diagnostics)")
        self.qc = read_json(system_dir / "qc.json")
        self.rmap = read_json(system_dir / "residue_map.json")
        self.mg = read_json(system_dir / "mg_restraints.json")
        selection = read_json(run_dir / "trajectory_atoms.json")["atom_indices"]
        if self.run["prmtop_sha256"] != sha256(system_dir / "complex.prmtop"):
            raise ValueError("complex.prmtop differs from the one used by the run")
        full_top = md.load_prmtop(str(system_dir / "complex.prmtop"))
        self.new_index = {old: new for new, old in enumerate(selection)}
        top = full_top.subset(selection)
        traj = md.load_dcd(str(run_dir / "trajectory.dcd"), top=top)
        start, interval, dt_fs = self.run["production_start_step"], self.run["report_interval_steps"], self.run["timestep_fs"]
        expected = (start + self.run["completed_production_steps"]) // interval - start // interval
        if len(traj) != expected:
            raise ValueError(f"trajectory has {len(traj)} frames, run.json implies {expected}")
        first = (start // interval + 1) * interval
        times = (first + np.arange(len(traj)) * interval - start) * dt_fs * 1e-6
        discard = 0.0 if self.run["smoke"] else protocol["analysis"]["discard_production_ns"]
        if discard_ns is not None:
            discard = discard_ns
        self.discard_ns = discard
        keep = times > discard
        self.traj, self.times = traj[keep], times[keep]
        if len(self.traj) == 0:
            raise ValueError("No frames after burn-in")
        self.frame_ps = interval * dt_fs * 1e-3
        ref = md.load(str(system_dir / "complex.inpcrd"), top=full_top).atom_slice(selection)
        self.reference = ref
        self.residues: dict[tuple[str, int], Any] = {}
        self.rows_by_index = {}
        top_res = list(self.traj.topology.residues)
        for row in self.rmap["residues"]:
            residue = top_res[row["prmtop_residue_index"]]
            if residue.name != row["resname_prmtop"] and not (residue.name == "HIS" and row["resname_prmtop"] in {"HID", "HIE", "HIP"}):
                raise ValueError(f"Residue mismatch at {row['prmtop_residue_index']}: {residue.name} != {row['resname_prmtop']}")
            self.rows_by_index[residue.index] = row
            if row["uniprot_number"] is not None:
                self.residues[(row["molecule"], row["uniprot_number"])] = residue
            elif row["molecule"] in {"GTP", "MG"}:
                self.residues[(row["molecule"], 0)] = residue
        self.groups = {name: [] for name in ["Rab10", "KIF13A", "GTP", "MG"]}
        for residue in top_res[:len(self.rmap["residues"])]:
            self.groups[self.rows_by_index[residue.index]["molecule"]] += [a.index for a in residue.atoms]
        self.shell_waters = [[self.new_index[i] for i in w.values()] for w in self.mg["shell_waters_unrestrained"]]
        self.ions = {"Na+": [a.index for a in self.traj.topology.atoms if a.residue.name in {"Na+", "NA"}],
                     "Cl-": [a.index for a in self.traj.topology.atoms if a.residue.name in {"Cl-", "CL"}]}
        self.image_complex()

    def chem(self, residue) -> str:
        """Amber residue name from the residue map (mdtraj turns HID/HIE/HIP into HIS)."""
        row = self.rows_by_index.get(residue.index)
        return row["resname_prmtop"] if row else residue.name

    def res(self, molecule: str, number: int):
        return self.residues[(molecule, number)]

    def atoms(self, molecule: str, number: int, names) -> list[int]:
        residue = self.res(molecule, number)
        found = [a.index for a in residue.atoms if a.name in names]
        if not found:
            raise KeyError(f"{molecule} {number} has none of {names}")
        return found

    def label(self, residue) -> str:
        row = self.rows_by_index.get(residue.index)
        if row is None:
            return f"{residue.name}{residue.index}"
        if row["molecule"] in {"GTP", "MG"}:
            return row["molecule"]
        if row["cap"]:
            return f"{row['molecule']} {residue.name}"
        return f"{row['molecule']} {ONE.get(row['resname_prmtop'], row['resname_prmtop'])}{row['uniprot_number']}"

    def image_complex(self) -> None:
        """Shift KIF13A, GTP, Mg and each shell water to the periodic copy nearest the Rab10 centre."""
        xyz = self.traj.xyz
        boxes = self.traj.unitcell_vectors
        if boxes is None:
            raise ValueError("Trajectory has no box vectors")
        anchor = xyz[:, self.groups["Rab10"]].mean(axis=1)
        cells = np.array(list(itertools.product((-1, 0, 1), repeat=3)), dtype=float)
        movers = [self.groups["KIF13A"], self.groups["GTP"], self.groups["MG"]] + self.shell_waters
        self.image_shifts = 0
        for group in movers:
            com = xyz[:, group].mean(axis=1)
            for f in range(len(xyz)):
                box = boxes[f]
                delta = com[f] - anchor[f]
                base = np.round(delta @ np.linalg.inv(box))
                candidates = (base + cells) @ box
                k = int(np.argmin(np.linalg.norm(delta - candidates, axis=1)))
                shift = candidates[k]
                if np.any(np.abs(shift) > 1e-6):
                    xyz[f, group] -= shift
                    self.image_shifts += 1


def min_dist(traj, left: list[int], right: list[int], periodic: bool = False) -> np.ndarray:
    import mdtraj as md
    pairs = np.array([(a, b) for a in left for b in right])
    return md.compute_distances(traj, pairs, periodic=periodic).min(axis=1)


def charged_groups(cx: Complex) -> list[dict]:
    groups = []
    for (molecule, number), residue in cx.residues.items():
        chem = cx.chem(residue)
        names = POS.get(chem) or NEG.get(chem)
        if molecule in {"Rab10", "KIF13A"} and names:
            groups.append({"molecule": molecule, "number": number, "label": cx.label(residue),
                           "sign": 1 if chem in POS else -1,
                           "atoms": [a.index for a in residue.atoms if a.name in names]})
    m1 = cx.res("Rab10", 1)
    groups.append({"molecule": "Rab10", "number": 1, "label": "Rab10 M1 N-terminus", "sign": 1,
                   "atoms": [a.index for a in m1.atoms if a.name == "N"]})
    gtp = cx.residues[("GTP", 0)]
    groups.append({"molecule": "GTP", "number": 0, "label": "GTP phosphates", "sign": -1,
                   "atoms": [a.index for a in gtp.atoms if a.name in GTP_NEG]})
    return groups


def helix_axis(xyz: np.ndarray) -> np.ndarray:
    centered = xyz - xyz.mean(axis=0)
    axis = np.linalg.svd(centered, full_matrices=False)[2][0]
    return axis if np.dot(axis, xyz[-1] - xyz[0]) > 0 else -axis


# ----------------------------------------------------------------------------- per-run analysis

def analyze_run(run_dir: Path, system_dir: Path, protocol: dict, *, allow_smoke: bool, gbsa: bool,
                solvent: bool, out_dir: Path | None = None, discard_ns: float | None = None) -> dict:
    import mdtraj as md
    cfg = protocol["analysis"]
    cx = Complex(run_dir, system_dir, protocol, allow_smoke, discard_ns)
    out_dir = out_dir or run_dir / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
    traj, times = cx.traj, cx.times
    sb_cut = cfg["salt_bridge_cutoff_nm"]
    res1122 = cx.res("KIF13A", 1122)
    name1122 = cx.chem(res1122)
    charged1122 = POS.get(name1122) or NEG.get(name1122)
    r70 = cx.atoms("Rab10", 70, POS["ARG"])
    e1130 = cx.atoms("KIF13A", 1130, NEG["GLU"])
    a1122 = cx.atoms("KIF13A", 1122, charged1122)
    series = {
        "e1130_r70_nm": min_dist(traj, e1130, r70),
        "res1122_r70_nm": min_dist(traj, a1122, r70),
        "res1122_e1130_nm": min_dist(traj, a1122, e1130),
    }
    endpoints: dict[str, Any] = {}
    dist = {key: describe(value, sb_cut) for key, value in series.items()}
    low, high, step = cfg["histogram_nm"]
    edges = np.arange(low, high + step / 2, step)
    histograms = {key: {"edges_nm": edges.round(4).tolist(),
                        "counts": np.histogram(np.clip(value, low, high), bins=edges)[0].tolist()}
                  for key, value in series.items()}
    endpoints["e1130_r70_sb_occupancy"] = dist["e1130_r70_nm"]["occupancy"]
    endpoints["e1130_r70_dist_nm_mean"] = dist["e1130_r70_nm"]["mean"]
    endpoints["res1122_r70_dist_nm_mean"] = dist["res1122_r70_nm"]["mean"]
    endpoints["res1122_r70_contact_occupancy"] = dist["res1122_r70_nm"]["occupancy"]
    endpoints["res1122_e1130_dist_nm_mean"] = dist["res1122_e1130_nm"]["mean"]
    endpoints["res1122_e1130_contact_occupancy"] = dist["res1122_e1130_nm"]["occupancy"]
    endpoints["k1122_e1130_sb_occupancy"] = dist["res1122_e1130_nm"]["occupancy"] if name1122 in POS else None
    endpoints["e1122_r70_sb_occupancy"] = dist["res1122_r70_nm"]["occupancy"] if name1122 in NEG else None

    # all intermolecular salt bridges: receptor (Rab10 + GTP) vs KIF13A
    groups = charged_groups(cx)
    receptor = [g for g in groups if g["molecule"] in {"Rab10", "GTP"}]
    ligand = [g for g in groups if g["molecule"] == "KIF13A"]
    sb_rows, sb_count = [], np.zeros(len(traj))
    for g1 in receptor:
        for g2 in ligand:
            if g1["sign"] == g2["sign"]:
                continue
            d = min_dist(traj, g1["atoms"], g2["atoms"])
            occ = d <= sb_cut
            sb_count += occ
            if occ.any() or d.min() < 0.6:
                sb_rows.append({"receptor": g1["label"], "kif13a": g2["label"], "occupancy": float(occ.mean()),
                                "min_dist_mean_nm": float(d.mean()), "min_dist_min_nm": float(d.min())})
    sb_rows.sort(key=lambda r: -r["occupancy"])
    write_tsv(out_dir / "intermolecular_salt_bridges.tsv", sb_rows)
    endpoints["intermolecular_sb_count_mean"] = float(sb_count.mean())
    series["intermolecular_sb_count"] = sb_count

    # interface residue contacts (closest heavy atom), strided
    stride = max(1, int(round(cfg["contact_stride_ps"] / cx.frame_ps)))
    sub = traj[::stride]
    rec_res = [r for r in traj.topology.residues if r.index in cx.rows_by_index and cx.rows_by_index[r.index]["molecule"] in {"Rab10", "GTP", "MG"}]
    lig_res = [r for r in traj.topology.residues if r.index in cx.rows_by_index and cx.rows_by_index[r.index]["molecule"] == "KIF13A"]

    def ref_atom(residue) -> int:
        names = {a.name: a.index for a in residue.atoms}
        return names.get("CA", names.get("PB", names.get("MG", names.get("C", residue.atom(0).index))))

    ca_pairs = np.array([(ref_atom(a), ref_atom(b)) for a in rec_res for b in lig_res])
    ca_d = md.compute_distances(sub, ca_pairs, periodic=False)
    candidate = ca_d.min(axis=0) < cfg["contact_candidate_ca_nm"] + 0.6 * np.array(
        [1.0 if cx.rows_by_index[a.index]["molecule"] in {"GTP"} else 0.0 for a in rec_res for _ in lig_res])
    res_pairs = [(a.index, b.index) for (a, b), keep in zip(itertools.product(rec_res, lig_res), candidate) if keep]
    contact_rows, per_res = [], defaultdict(lambda: np.zeros(len(sub), dtype=bool))
    pair_count = np.zeros(len(sub))
    if res_pairs:
        cd, pairs_used = md.compute_contacts(sub, res_pairs, scheme="closest-heavy", periodic=False)
        contact = cd <= cfg["contact_cutoff_nm"]
        pair_count = contact.sum(axis=1)
        top_res = list(traj.topology.residues)
        for k, (i, j) in enumerate(pairs_used):
            if contact[:, k].any():
                contact_rows.append({"receptor": cx.label(top_res[i]), "kif13a": cx.label(top_res[j]),
                                     "occupancy": float(contact[:, k].mean()), "min_dist_mean_nm": float(cd[:, k].mean())})
                per_res[cx.label(top_res[i])] |= contact[:, k]
                per_res[cx.label(top_res[j])] |= contact[:, k]
    contact_rows.sort(key=lambda r: -r["occupancy"])
    write_tsv(out_dir / "interface_contacts.tsv", contact_rows)
    interface = sorted(({"residue": k, "occupancy": float(v.mean())} for k, v in per_res.items()), key=lambda r: -r["occupancy"])
    write_tsv(out_dir / "interface_residues.tsv", interface)
    endpoints["interface_contact_pairs_mean"] = float(pair_count.mean())
    endpoints["interface_residues_occ50"] = int(sum(1 for r in interface if r["occupancy"] >= 0.5))

    # buried SASA and H72 side-chain SASA
    sstride = max(1, int(round(cfg["sasa_stride_ps"] / cx.frame_ps)))
    ss = traj[::sstride]
    rec_atoms = cx.groups["Rab10"] + cx.groups["GTP"] + cx.groups["MG"]
    lig_atoms = cx.groups["KIF13A"]
    both = sorted(rec_atoms + lig_atoms)
    sasa_c = md.shrake_rupley(ss.atom_slice(both), mode="atom")
    sasa_r = md.shrake_rupley(ss.atom_slice(sorted(rec_atoms)), mode="atom").sum(axis=1)
    sasa_l = md.shrake_rupley(ss.atom_slice(sorted(lig_atoms)), mode="atom").sum(axis=1)
    buried = sasa_r + sasa_l - sasa_c.sum(axis=1)
    endpoints["buried_sasa_nm2_mean"] = float(buried.mean())
    strided = {"contact_pairs": (stride, pair_count), "buried_sasa_nm2": (sstride, buried)}
    position_in_both = {old: k for k, old in enumerate(both)}
    h72 = cx.res("Rab10", 72)
    h72_side = [position_in_both[a.index] for a in h72.atoms if a.name not in {"N", "CA", "C", "O", "H", "HA"}]
    h72_sasa = sasa_c[:, h72_side].sum(axis=1)

    # KIF13A turn DSSP and interhelix angle
    kif = traj.atom_slice(cx.groups["KIF13A"])
    dssp = md.compute_dssp(kif, simplified=True)
    kif_numbers = [cx.rows_by_index[r.index]["uniprot_number"] for r in traj.topology.residues if r.index in cx.rows_by_index
                   and cx.rows_by_index[r.index]["molecule"] == "KIF13A"]
    t0, t1 = cfg["turn_dssp"]
    cols = [k for k, n in enumerate(kif_numbers) if n is not None and t0 <= n <= t1]
    turn_h = (dssp[:, cols] == "H")
    series["turn_helix_fraction"] = turn_h.mean(axis=1)
    endpoints["turn_helix_fraction_1116_1126"] = float(turn_h.mean())
    dssp_rows = [{"kif13a": kif_numbers[k], "helix_H": float((dssp[:, k] == "H").mean()),
                  "strand_E": float((dssp[:, k] == "E").mean()), "coil_C": float((dssp[:, k] == "C").mean())}
                 for k in range(len(kif_numbers)) if kif_numbers[k] is not None]
    write_tsv(out_dir / "kif13a_dssp.tsv", dssp_rows)
    h1 = [cx.atoms("KIF13A", n, ["CA"])[0] for n in range(cfg["helix1"][0], cfg["helix1"][1] + 1)]
    h2 = [cx.atoms("KIF13A", n, ["CA"])[0] for n in range(cfg["helix2"][0], cfg["helix2"][1] + 1)]
    angles = np.array([math.degrees(math.acos(np.clip(np.dot(helix_axis(f[h1]), helix_axis(f[h2])), -1, 1)))
                       for f in traj.xyz])
    series["interhelix_angle_deg"] = angles
    endpoints["interhelix_angle_deg_mean"] = float(angles.mean())

    # RMSD after superposition on the Rab10 core (reference = prepared start)
    c0, c1 = cfg["rab10_core"]
    excluded = {n for lo, hi in cfg["rab10_core_exclude"] for n in range(lo, hi + 1)}
    core = [cx.atoms("Rab10", n, ["CA"])[0] for n in range(c0, c1 + 1) if n not in excluded]
    fitted = traj.slice(range(len(traj)), copy=True)
    fitted.superpose(cx.reference, atom_indices=core, ref_atom_indices=core)

    def rmsd(indices: list[int]) -> np.ndarray:
        diff = fitted.xyz[:, indices] - cx.reference.xyz[0, indices]
        return np.sqrt((diff ** 2).sum(axis=2).mean(axis=1))

    ca = lambda mol, lo, hi: [cx.atoms(mol, n, ["CA"])[0] for n in range(lo, hi + 1)]
    series["rab10_core_ca_rmsd_nm"] = rmsd(core)
    series["kif_helix2_rmsd_nm"] = rmsd(ca("KIF13A", *cfg["helix2"]))
    series["rab10_switch1_rmsd_nm"] = rmsd(ca("Rab10", *cfg["rab10_switch1"]))
    series["rab10_switch2_rmsd_nm"] = rmsd(ca("Rab10", *cfg["rab10_switch2"]))
    series["kif_hairpin_ca_rmsd_nm"] = rmsd(ca("KIF13A", cfg["helix1"][0], cfg["helix2"][1]))
    for key in ["rab10_core_ca_rmsd_nm", "kif_helix2_rmsd_nm", "rab10_switch1_rmsd_nm", "rab10_switch2_rmsd_nm",
                "kif_hairpin_ca_rmsd_nm"]:
        endpoints[f"{key}_mean"] = float(series[key].mean())

    # Na+ near acidic KIF13A residues (periodic minimum image)
    na_rows = []
    for number in cfg["na_residues"]:
        residue = cx.res("KIF13A", number)
        if cx.chem(residue) not in NEG or not cx.ions["Na+"]:
            endpoints[f"na_occupancy_{number}"] = None
            continue
        oxy = [a.index for a in residue.atoms if a.name in NEG[cx.chem(residue)]]
        d = md.compute_distances(traj, np.array([(o, n) for o in oxy for n in cx.ions["Na+"]]), periodic=True)
        d = d.reshape(len(traj), len(oxy), len(cx.ions["Na+"])).min(axis=1)
        bound = d <= cfg["na_cutoff_nm"]
        occ = float(bound.any(axis=1).mean())
        endpoints[f"na_occupancy_{number}"] = occ
        na_rows.append({"kif13a": f"{ONE[cx.chem(residue)]}{number}", "occupancy_any_na": occ,
                        "mean_na_count": float(bound.sum(axis=1).mean())})
    write_tsv(out_dir / "na_contacts.tsv", na_rows)

    # Rab10 H72 environment
    ring = [a.index for a in h72.atoms if a.name in {"ND1", "NE2"}]
    gtp = cx.residues[("GTP", 0)]
    gamma = [a.index for a in gtp.atoms if a.name in {"O1G", "O2G", "O3G"}]
    h72_env = {"state": cx.chem(h72),
               "ring_to_R70_nm": describe(min_dist(traj, ring, r70)),
               "ring_to_E1130_nm": describe(min_dist(traj, ring, e1130)),
               "ring_to_res1122_nm": describe(min_dist(traj, ring, a1122)),
               "ring_to_GTP_gamma_O_nm": describe(min_dist(traj, ring, gamma)),
               "sidechain_sasa_nm2": describe(h72_sasa)}
    polar = [a.index for a in traj.topology.atoms if a.element is not None and a.element.symbol in {"N", "O"}
             and a.residue.index != h72.index and a.residue.name not in {"Na+", "Cl-", "NA", "CL"}]
    hb = defaultdict(int)
    hbond_any = np.zeros(len(traj), dtype=bool)
    for ring_atom in ring:
        pairs = np.array([(ring_atom, p) for p in polar])
        d = md.compute_distances(traj, pairs, periodic=False)
        close = d <= cfg["hbond_heavy_cutoff_nm"]
        hbond_any |= close.any(axis=1)
        atoms_list = list(traj.topology.atoms)
        for k in np.where(close.any(axis=0))[0]:
            partner = atoms_list[polar[k]]
            hb[f"{traj.topology.atom(ring_atom).name}...{cx.label(partner.residue)}:{partner.name}"] += int(close[:, k].sum())
    h72_env["ring_polar_partner_fraction"] = float(hbond_any.mean())
    h72_env["top_partners"] = [{"pair": k, "fraction": v / len(traj)} for k, v in Counter(hb).most_common(6)]
    endpoints["h72_ring_to_res1122_nm_mean"] = h72_env["ring_to_res1122_nm"]["mean"]
    endpoints["h72_ring_to_e1130_nm_mean"] = h72_env["ring_to_E1130_nm"]["mean"]
    endpoints["h72_ring_to_r70_nm_mean"] = h72_env["ring_to_R70_nm"]["mean"]
    endpoints["h72_polar_partner_fraction"] = h72_env["ring_polar_partner_fraction"]
    endpoints["h72_sidechain_sasa_nm2_mean"] = h72_env["sidechain_sasa_nm2"]["mean"]

    # Mg coordination
    mg_atom = cx.groups["MG"][0]
    mg_rows = []
    for pair in cx.mg["pairs"]:
        lig = cx.new_index[pair["ligand_index"]]
        d = md.compute_distances(traj, np.array([[mg_atom, lig]]), periodic=True)[:, 0]
        mg_rows.append({"ligand": f"{pair['residue']}:{pair['atom']}", "flat_bottom_nm": pair["flat_bottom_nm"],
                        "mean_nm": float(d.mean()), "max_nm": float(d.max()),
                        "fraction_within_flat_bottom": float((d <= pair["flat_bottom_nm"]).mean())})
    oxygens = [a.index for a in traj.topology.atoms if a.element is not None and a.element.symbol == "O"
               and a.residue.name not in {"Na+", "Cl-"}]
    d = md.compute_distances(traj, np.array([[mg_atom, o] for o in oxygens]), periodic=True)
    cn = (d <= cfg["mg_cutoff_nm"]).sum(axis=1)
    shell_o = [w[[traj.topology.atom(i).name for i in w].index("O")] for w in cx.shell_waters]
    shell_d = md.compute_distances(traj, np.array([[mg_atom, o] for o in shell_o]), periodic=True)
    mg_info = {"restrained_pairs": mg_rows, "cn_tracked_O": describe(cn),
               "cn_tracked_values": {int(k): int(v) for k, v in zip(*np.unique(cn, return_counts=True))},
               "shell_water_O_nm": [describe(shell_d[:, k]) for k in range(len(shell_o))],
               "note": "tracked oxygens = protein, caps, GTP and the two placed shell waters; other waters are only in the solvent DCD (--solvent)"}
    endpoints["mg_cn_tracked_mean"] = float(cn.mean())
    endpoints["mg_restrained_max_distance_nm"] = float(max(r["max_nm"] for r in mg_rows))
    if solvent:
        mg_info["cn_all_water_solvent_dcd"] = solvent_mg_cn(cx, cfg)

    # periodic self-image distance of the re-imaged complex
    istride = max(1, int(round(cfg["image_stride_ps"] / cx.frame_ps)))
    probe = [a.index for a in traj.topology.atoms if a.name == "CA" or a.index == mg_atom
             or (a.residue.index == gtp.index and a.name in {"PA", "PB", "PG"})]
    image_min = []
    for f in range(0, len(traj), istride):
        x = traj.xyz[f, probe]
        box = traj.unitcell_vectors[f]
        best = np.inf
        for cell in itertools.product((-1, 0, 1), repeat=3):
            if cell == (0, 0, 0):
                continue
            shifted = x + np.array(cell, dtype=float) @ box
            best = min(best, float(np.min(np.linalg.norm(x[:, None] - shifted[None], axis=-1))))
        image_min.append(best)
    image_min = np.array(image_min)
    endpoints["image_min_distance_nm"] = float(image_min.min())

    gbsa_info = None
    if gbsa:
        gbsa_info = gbsa_energies(cx, protocol)
        endpoints["gbsa_interaction_kj_mol_mean"] = gbsa_info["mean_kj_mol"]

    rows = [{"production_ns": float(t), **{k: float(v[i]) for k, v in series.items() if len(v) == len(traj)}}
            for i, t in enumerate(times)]
    write_tsv(out_dir / "timeseries.tsv", rows)
    strided_rows = []
    for key, (every, values) in strided.items():
        strided_rows += [{"series": key, "production_ns": float(times[k * every]), "value": float(v)} for k, v in enumerate(values)]
    write_tsv(out_dir / "timeseries_strided.tsv", strided_rows)
    eligible = (not cx.run["smoke"] and cx.run["duration_override_ns"] is None and discard_ns is None
                and cx.run["completed_production_steps"] == cx.run["target_production_steps"])
    summary = {
        "system_id": cx.run["system_id"], "arm": cx.run["arm"], "start": cx.run["start"], "replica": cx.run["replica"],
        "pair_key": cx.run["pair_key"], "residue_1122": name1122, "seed": cx.run["seed"],
        "eligible_for_inference": eligible, "smoke": cx.run["smoke"], "utc": utc_now(),
        "frames": len(traj), "analysis_window_ns": [float(times[0]), float(times[-1])], "discarded_ns": cx.discard_ns,
        "image_shift_events": cx.image_shifts, "protocol_sha256": cx.run["protocol_sha256"],
        "prmtop_sha256": cx.run["prmtop_sha256"],
        "his_states": [{k: h[k] for k in ("residue", "propka_pka", "state", "flag")} for h in cx.qc["his_states"]],
        "endpoints": endpoints, "distributions": dist, "histograms": histograms,
        "h72_environment": h72_env, "mg": mg_info,
        "image": {"min_nm": float(image_min.min()), "fraction_below_warn": float((image_min < cfg["image_warn_nm"]).mean()),
                  "warn_nm": cfg["image_warn_nm"]},
        "gbsa": gbsa_info,
        "definitions": {"salt_bridge": f"min distance between charged-group atoms <= {sb_cut} nm ({protocol['analysis']['charge_atoms_note']})",
                        "contacts": f"closest heavy atoms <= {cfg['contact_cutoff_nm']} nm every {cfg['contact_stride_ps']} ps",
                        "sasa": f"Shrake-Rupley (mdtraj), every {cfg['sasa_stride_ps']} ps; buried = SASA(Rab10+GTP+Mg) + SASA(KIF13A) - SASA(complex)",
                        "rmsd": "CA RMSD to the prepared start after superposition on Rab10 core CA (6-170 without switch I 36-44 and switch II 66-80)",
                        "dssp": "mdtraj simplified DSSP, H = helix"},
    }
    if summary["image"]["fraction_below_warn"] > 0:
        summary["quality_warning"] = "Complex within image_warn_nm of its periodic image in some frames"
    write_json(out_dir / "analysis.json", summary)
    inputs = [run_dir / "run.json", run_dir / "trajectory.dcd", system_dir / "complex.prmtop", system_dir / "residue_map.json"]
    for name in ["analysis.json", "timeseries.tsv", "intermolecular_salt_bridges.tsv", "interface_contacts.tsv"]:
        provenance(out_dir / name, inputs, params={"analysis": cfg, "allow_smoke": allow_smoke})
    return summary


def solvent_mg_cn(cx: Complex, cfg: dict) -> dict:
    import mdtraj as md
    path = Path(cx.run["solvent_dir"]) / "system_solvent.dcd"
    if not path.exists():
        return {"available": False}
    top = md.load_prmtop(str(cx.system_dir / "complex.prmtop"))
    mg_atom = [a.index for a in top.atoms if a.residue.name == "MG"][0]
    water_o = [a.index for a in top.atoms if a.residue.name in {"WAT", "HOH"} and a.name == "O"]
    sel = [mg_atom] + water_o
    counts = []
    for chunk in md.iterload(str(path), top=top, chunk=50, atom_indices=sel):
        d = md.compute_distances(chunk, np.array([[0, k] for k in range(1, len(sel))]), periodic=True)
        counts.append((d <= cfg["mg_cutoff_nm"]).sum(axis=1))
    counts = np.concatenate(counts)
    protein_gtp = len(cx.mg["pairs"])
    return {"available": True, "frames": len(counts), "water_O_within_cutoff": describe(counts),
            "values": {int(k): int(v) for k, v in zip(*np.unique(counts, return_counts=True))},
            "note": f"plus {protein_gtp} protein/GTP ligands if they stay bound (see restrained_pairs)"}


def gbsa_energies(cx: Complex, protocol: dict) -> dict:
    """GBn2 (mbondi3) single-trajectory interaction energy, frames every gbsa_stride_ps. Relative only."""
    import openmm as mm
    from openmm import app, unit
    cfg = protocol["analysis"]
    stride = max(1, int(round(cfg["gbsa_stride_ps"] / cx.frame_ps)))
    frames = cx.traj[::stride]
    salt = protocol["solvation"]["salt_molar"] * unit.molar
    temperature = protocol["simulation"]["temperature_k"] * unit.kelvin
    rows = cx.rmap["residues"]
    rec = [i for r in rows if r["prmtop_residue_index"] in set(cx.rmap["receptor_residue_indices"]) for i in range(r["first_atom"], r["end_atom"])]
    lig = [i for r in rows if r["prmtop_residue_index"] in set(cx.rmap["ligand_residue_indices"]) for i in range(r["first_atom"], r["end_atom"])]
    n = cx.rmap["n_solute_atoms"]
    contexts = {}
    for name, atoms in [("complex", list(range(n))), ("receptor", rec), ("ligand", lig)]:
        prmtop = app.AmberPrmtopFile(str(cx.system_dir / f"{name}_dry.prmtop"))
        if prmtop.topology.getNumAtoms() != len(atoms):
            raise ValueError(f"{name}_dry.prmtop atom count mismatch")
        system = prmtop.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, implicitSolvent=app.GBn2,
                                     implicitSolventSaltConc=salt, temperature=temperature)
        context = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("CPU"))
        contexts[name] = (context, atoms)
    energies = []
    for f in range(len(frames)):
        e = {}
        for name, (context, atoms) in contexts.items():
            context.setPositions(frames.xyz[f, atoms])
            e[name] = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        energies.append(e["complex"] - e["receptor"] - e["ligand"])
    energies = np.array(energies)
    return {"frames": len(energies), "stride_ps": cfg["gbsa_stride_ps"], "mean_kj_mol": float(energies.mean()),
            "sd_kj_mol": float(energies.std()), "values_kj_mol": energies.round(2).tolist(),
            "label": "relative interaction energy only (GBn2/mbondi3, ff19SB, 0.15 M, no entropy, no SA beyond the default ACE term); not a binding free energy"}


# ----------------------------------------------------------------------------- aggregation

def t_ci(values: np.ndarray, df: float | None = None) -> tuple[float, float] | tuple[None, None]:
    from scipy import stats
    n = len(values)
    if n < 2:
        return None, None
    df = n - 1 if df is None else df
    half = stats.t.ppf(0.975, df) * values.std(ddof=1) / math.sqrt(n)
    return float(values.mean() - half), float(values.mean() + half)


def welch(a: np.ndarray, b: np.ndarray) -> dict:
    from scipy import stats
    diff = float(a.mean() - b.mean())
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    se = math.sqrt(va + vb)
    if se == 0:
        return {"difference": diff, "se": 0.0, "df": None, "ci95": [diff, diff], "t": None}
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    half = stats.t.ppf(0.975, df) * se
    return {"difference": diff, "se": se, "df": df, "ci95": [diff - half, diff + half], "t": diff / se}


def permutation_exact(a: np.ndarray, b: np.ndarray) -> dict:
    pooled = np.concatenate([a, b])
    observed = abs(a.mean() - b.mean())
    diffs = []
    for idx in itertools.combinations(range(len(pooled)), len(a)):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(idx)] = True
        diffs.append(abs(pooled[mask].mean() - pooled[~mask].mean()))
    diffs = np.array(diffs)
    tol = 1e-12 * max(1.0, observed)
    p = float((diffs >= observed - tol).mean())
    p_min_data = float((diffs >= diffs.max() - tol).mean())
    total = len(diffs)
    p_min_design = (2 if len(a) == len(b) else 1) / total
    return {"p_two_sided": p, "relabelings": total, "p_min_attainable_design": p_min_design, "p_min_attainable_these_values": p_min_data}


def signflip_exact(d: np.ndarray) -> dict:
    observed = abs(d.mean())
    stats_ = [abs((d * np.array(signs)).mean()) for signs in itertools.product((-1, 1), repeat=len(d))]
    stats_ = np.array(stats_)
    tol = 1e-12 * max(1.0, observed)
    return {"p_two_sided": float((stats_ >= observed - tol).mean()), "sign_assignments": len(stats_),
            "p_min_attainable_design": 2 / len(stats_)}


def aggregate(runs_root: Path, out_dir: Path, protocol: dict, protocol_path: Path, allow_smoke: bool) -> dict:
    analyses = []
    for path in sorted(runs_root.glob("*/analysis/analysis.json")):
        item = read_json(path)
        if item["smoke"] and not allow_smoke:
            continue
        if not item["eligible_for_inference"] and not allow_smoke:
            continue
        analyses.append(item)
    if not analyses:
        raise ValueError("No analyses to aggregate")
    if len({a["protocol_sha256"] for a in analyses}) != 1:
        raise ValueError("Cannot combine runs from different protocol files")
    by_arm = defaultdict(list)
    for item in analyses:
        by_arm[item["arm"]].append(item)
    for arm, items in by_arm.items():
        if len({i["start"] for i in items}) != len(items) or len({i["seed"] for i in items}) != len(items):
            raise ValueError(f"{arm}: duplicate starts or seeds")
    endpoints = sorted({k for item in analyses for k, v in item["endpoints"].items() if isinstance(v, (int, float))})
    group_rows, diff_rows = [], []
    for endpoint in endpoints:
        for arm in protocol["arms"]:
            items = [i for i in by_arm.get(arm, []) if i["endpoints"].get(endpoint) is not None]
            if not items:
                continue
            values = np.array([i["endpoints"][endpoint] for i in items], dtype=float)
            low, high = t_ci(values)
            group_rows.append({"endpoint": endpoint, "arm": arm, "n": len(values),
                               "replicas": ";".join(f"{i['start']}={i['endpoints'][endpoint]:.6g}" for i in items),
                               "mean": float(values.mean()), "sd": float(values.std(ddof=1)) if len(values) > 1 else None,
                               "t_ci95_low": low, "t_ci95_high": high})
        for arm_a, arm_b, kind in protocol["analysis"]["comparisons"]:
            items_a = [i for i in by_arm.get(arm_a, []) if i["endpoints"].get(endpoint) is not None]
            items_b = [i for i in by_arm.get(arm_b, []) if i["endpoints"].get(endpoint) is not None]
            if len(items_a) < 2 or len(items_b) < 2:
                continue
            row = {"endpoint": endpoint, "contrast": f"{arm_a} - {arm_b}", "kind": kind}
            if kind == "welch":
                a = np.array([i["endpoints"][endpoint] for i in items_a], dtype=float)
                b = np.array([i["endpoints"][endpoint] for i in items_b], dtype=float)
                w, perm = welch(a, b), permutation_exact(a, b)
                row.update({"n_a": len(a), "n_b": len(b), "difference": w["difference"], "se": w["se"], "df": w["df"],
                            "ci95_low": w["ci95"][0], "ci95_high": w["ci95"][1], "p_exact": perm["p_two_sided"],
                            "p_min_design": perm["p_min_attainable_design"],
                            "p_min_these_values": perm["p_min_attainable_these_values"], "permutations": perm["relabelings"]})
            else:
                keyed_a = {i["pair_key"]: i["endpoints"][endpoint] for i in items_a}
                keyed_b = {i["pair_key"]: i["endpoints"][endpoint] for i in items_b}
                keys = sorted(set(keyed_a) & set(keyed_b) - {None})
                if len(keys) < 2:
                    continue
                d = np.array([keyed_a[k] - keyed_b[k] for k in keys], dtype=float)
                low, high = t_ci(d)
                flip = signflip_exact(d)
                row.update({"n_a": len(keys), "n_b": len(keys), "difference": float(d.mean()),
                            "se": float(d.std(ddof=1) / math.sqrt(len(d))), "df": len(d) - 1,
                            "ci95_low": low, "ci95_high": high, "p_exact": flip["p_two_sided"],
                            "p_min_design": flip["p_min_attainable_design"], "p_min_these_values": None,
                            "permutations": flip["sign_assignments"],
                            "pairs": ";".join(f"{k}={v:.6g}" for k, v in zip(keys, d))})
            diff_rows.append(row)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(out_dir / "aggregate_groups.tsv", group_rows)
    write_tsv(out_dir / "aggregate_differences.tsv", diff_rows)
    summary = {"utc": utc_now(), "runs": [{k: a[k] for k in ("system_id", "arm", "start", "eligible_for_inference", "smoke")} for a in analyses],
               "protocol_sha256": analyses[0]["protocol_sha256"], "groups": group_rows, "differences": diff_rows,
               "inference_note": "Replicas are independent starting models (three per arm); frames are never treated as replicates. "
                                 "With 3 vs 3 the exact two-sided permutation p cannot go below 0.10 (sign-flip with 3 pairs: 0.25); "
                                 "t intervals with 2-4 df are wide. Treat results as exploratory effect sizes."}
    write_json(out_dir / "aggregate.json", summary)
    provenance(out_dir / "aggregate.json", [protocol_path] + [runs_root / a["system_id"] / "analysis" / "analysis.json" for a in analyses],
               params={"allow_smoke": allow_smoke})
    return summary


# ----------------------------------------------------------------------------- self-test

def selftest(out_dir: Path, system_dir: Path, protocol: dict, protocol_path: Path) -> dict:
    """Synthetic run built from a prepared system: frames 0-9 = start; frames 10-19 = KIF13A pulled 1.5 nm away.

    Expected: E1130-R70 occupancy 0.5 (start distance is a salt bridge), buried SASA ~0 and no interface
    contacts in the pulled half; aggregation statistics checked against hand-computed values.
    """
    import mdtraj as md
    from mdtraj.formats import DCDTrajectoryFile
    run_dir = out_dir / "synthetic_run"
    run_dir.mkdir(parents=True, exist_ok=True)
    qc = read_json(system_dir / "qc.json")
    mg = read_json(system_dir / "mg_restraints.json")
    rmap = read_json(system_dir / "residue_map.json")
    top = md.load_prmtop(str(system_dir / "complex.prmtop"))
    full = md.load(str(system_dir / "complex.inpcrd"), top=top)
    ions = [a.index for a in top.atoms if a.residue.name in {"Na+", "Cl-"}]
    shell = sorted(i for w in mg["shell_waters_unrestrained"] for i in w.values())
    selection = sorted(set(range(rmap["n_solute_atoms"])) | set(shell) | set(ions))
    sub = full.atom_slice(selection)
    kif = [k for k, old in enumerate(selection) if old < rmap["n_solute_atoms"] and any(
        r["molecule"] == "KIF13A" and r["first_atom"] <= old < r["end_atom"] for r in rmap["residues"])]
    rab = [k for k, old in enumerate(selection) if any(r["molecule"] == "Rab10" and r["first_atom"] <= old < r["end_atom"] for r in rmap["residues"])]
    direction = sub.xyz[0, kif].mean(axis=0) - sub.xyz[0, rab].mean(axis=0)
    direction /= np.linalg.norm(direction)
    frames = np.repeat(sub.xyz, 20, axis=0)
    frames[10:, kif] += 1.5 * direction
    lengths = np.repeat(sub.unitcell_lengths * 10, 20, axis=0)
    angles = np.repeat(sub.unitcell_angles, 20, axis=0)
    with DCDTrajectoryFile(str(run_dir / "trajectory.dcd"), "w") as handle:
        handle.write(frames * 10, cell_lengths=lengths, cell_angles=angles)
    write_json(run_dir / "trajectory_atoms.json", {"atom_indices": selection})
    run = {"system_id": "SYNTHETIC_" + qc["system"]["system_id"], "arm": "WT", "start": "synthetic", "replica": 1,
           "pair_key": None, "seed": 0, "completed": True, "smoke": True, "duration_override_ns": None,
           "production_start_step": 0, "report_interval_steps": 10, "timestep_fs": 2.0,
           "completed_production_steps": 200, "target_production_steps": 200,
           "protocol_sha256": sha256(protocol_path), "prmtop_sha256": sha256(system_dir / "complex.prmtop"),
           "solvent_dir": str(out_dir / "none")}
    write_json(run_dir / "run.json", run)
    fine = json.loads(json.dumps(protocol))
    fine["analysis"].update({"contact_stride_ps": 0.02, "sasa_stride_ps": 0.02, "image_stride_ps": 0.02, "gbsa_stride_ps": 0.2})
    summary = analyze_run(run_dir, system_dir, fine, allow_smoke=True, gbsa=True, solvent=False)
    strided = list(csv.DictReader((run_dir / "analysis" / "timeseries_strided.tsv").open(), delimiter="	"))
    buried = [float(r["value"]) for r in strided if r["series"] == "buried_sasa_nm2"]
    contacts = [float(r["value"]) for r in strided if r["series"] == "contact_pairs"]
    ts = list(csv.DictReader((run_dir / "analysis" / "timeseries.tsv").open(), delimiter="\t"))
    e_first = [float(r["e1130_r70_nm"]) for r in ts[:10]]
    e_second = [float(r["e1130_r70_nm"]) for r in ts[10:]]
    checks = {
        "e1130_r70_occupancy_is_0.5": abs(summary["endpoints"]["e1130_r70_sb_occupancy"] - 0.5) < 1e-9,
        "bridge_present_at_start": max(e_first) <= 0.4,
        "bridge_absent_when_pulled": min(e_second) > 1.0,
        "contacts_only_in_bound_half": all(float(r["intermolecular_sb_count"]) == 0 for r in ts[10:]),
        "gbsa_bound_more_negative": summary["gbsa"]["values_kj_mol"][0] < summary["gbsa"]["values_kj_mol"][-1] - 100,
        "gbsa_pulled_near_zero": abs(summary["gbsa"]["values_kj_mol"][-1]) < 50,
        "buried_sasa_positive_bound": min(buried[:10]) > 5.0,
        "buried_sasa_zero_pulled": max(abs(x) for x in buried[10:]) < 0.05,
        "interface_contacts_zero_pulled": max(contacts[10:]) == 0 and min(contacts[:10]) > 10,
        "rmsd_zero_for_start_frames": float(np.max([float(r["kif_helix2_rmsd_nm"]) for r in ts[:10]])) < 1e-3,
        "helix2_rmsd_equals_pull": all(abs(float(r["kif_helix2_rmsd_nm"]) - 1.5) < 1e-3 for r in ts[10:]),
    }
    # aggregation statistics on hand-made values
    a, b = np.array([1.0, 2.0, 3.0]), np.array([4.0, 5.0, 6.0])
    perm = permutation_exact(a, b)
    w = welch(a, b)
    flip = signflip_exact(np.array([1.0, 2.0, 3.0]))
    checks.update({
        "perm_p_extreme_is_0.1": abs(perm["p_two_sided"] - 0.1) < 1e-12 and perm["relabelings"] == 20,
        "perm_p_min_design_0.1": abs(perm["p_min_attainable_design"] - 0.1) < 1e-12,
        "welch_diff_-3_df_4": abs(w["difference"] + 3) < 1e-12 and abs(w["df"] - 4) < 1e-9,
        "welch_ci_matches_t": abs(w["ci95"][1] - (-3 + 2.7764451051977987 * math.sqrt(2 / 3))) < 1e-9,
        "signflip_p_0.25": abs(flip["p_two_sided"] - 0.25) < 1e-12,
    })
    checks = {k: (None if v is None else bool(v)) for k, v in checks.items()}
    result = {"checks": checks, "all_passed": all(v for v in checks.values() if v is not None),
              "synthetic_endpoints": summary["endpoints"]}
    write_json(out_dir / "selftest.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    sub = parser.add_subparsers(dest="command", required=True)
    one = sub.add_parser("run")
    one.add_argument("--run-dir", type=Path, required=True)
    one.add_argument("--system-dir", type=Path, required=True)
    one.add_argument("--out-dir", type=Path)
    one.add_argument("--allow-smoke", action="store_true")
    one.add_argument("--gbsa", action="store_true")
    one.add_argument("--solvent", action="store_true", help="also count Mg-water O from the full-system DCD")
    one.add_argument("--discard-ns", type=float, help="override the burn-in (pilots); the run is then not eligible for inference")
    agg = sub.add_parser("aggregate")
    agg.add_argument("--runs-root", type=Path, required=True)
    agg.add_argument("--out-dir", type=Path, required=True)
    agg.add_argument("--allow-smoke", action="store_true")
    st = sub.add_parser("selftest")
    st.add_argument("--out-dir", type=Path, required=True)
    st.add_argument("--system-dir", type=Path, required=True)
    args = parser.parse_args()
    protocol = load_protocol(args.protocol)
    if args.command == "run":
        summary = analyze_run(args.run_dir, args.system_dir, protocol, allow_smoke=args.allow_smoke,
                              gbsa=args.gbsa, solvent=args.solvent, out_dir=args.out_dir, discard_ns=args.discard_ns)
        print(json.dumps({"analyzed": summary["system_id"], "frames": summary["frames"],
                          "eligible_for_inference": summary["eligible_for_inference"],
                          "endpoints": summary["endpoints"]}, indent=1))
    elif args.command == "aggregate":
        summary = aggregate(args.runs_root, args.out_dir, protocol, args.protocol, args.allow_smoke)
        print(json.dumps({"runs": len(summary["runs"]), "differences": len(summary["differences"])}))
    else:
        result = selftest(args.out_dir, args.system_dir, protocol, args.protocol)
        print(json.dumps(result, indent=1))
        if not result["all_passed"]:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
