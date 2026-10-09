"""Analyze frozen local monomer endpoints and compare independent replica means.

    python scripts/k1122e_md/analyze.py --run-root RUN --output-root RESULTS

Coordinates in the protein-only DCD must match protein.pdb and the residue map.
The alpha-torsion endpoint is an explicit geometric definition, not DSSP.
Absolute monomer/solvated energies are excluded from WT-mutant contrasts.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

from common import DEFAULT_CONFIG, load_config, provenance, read_json, sha256, write_json

CHARGE_ATOMS = {
    "LYS": (1, ["NZ"]),
    "ARG": (1, ["NH1", "NH2"]),
    "ASP": (-1, ["OD1", "OD2"]),
    "GLU": (-1, ["OE1", "OE2"]),
}


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("Refusing to write a table with no rows")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def autocorrelation_quality(values):
    """Initial-positive-sequence estimate for diagnostics within one trajectory."""
    import numpy as np
    values = np.asarray(values, dtype=float)
    centered = values - values.mean()
    if len(values) < 4 or float(np.var(values)) == 0:
        return {"integrated_autocorrelation_frames": 1.0,
                "effective_frames_diagnostic": float(len(values)),
                "constant_or_short_trace": True}
    size = 1 << (2 * len(values) - 1).bit_length()
    spectrum = np.fft.rfft(centered, n=size)
    covariance = np.fft.irfft(spectrum * np.conjugate(spectrum), n=size)[:len(values)]
    covariance /= np.arange(len(values), 0, -1)
    correlation = covariance / covariance[0]
    paired = correlation[1:len(values) - (len(values) % 2)]
    paired = paired[:len(paired) // 2 * 2].reshape(-1, 2).sum(axis=1)
    stop = next((i for i, item in enumerate(paired) if item <= 0), len(paired))
    tau = max(1.0, float(1 + 2 * paired[:stop].sum()))
    return {"integrated_autocorrelation_frames": tau,
            "effective_frames_diagnostic": float(len(values) / tau),
            "constant_or_short_trace": False}


def min_atom_distances(xyz, left: list[int], right: list[int]):
    import numpy as np
    result = np.empty(len(xyz), dtype=float)
    for offset in range(0, len(xyz), 512):
        block = xyz[offset:offset + 512]
        distances = np.linalg.norm(block[:, left, None, :] - block[:, None, right, :], axis=-1)
        result[offset:offset + len(block)] = distances.min(axis=(1, 2))
    return result


def alpha_fraction(traj, residues: dict[int, Any], positions: list[int], settings: dict):
    import mdtraj as md
    import numpy as np
    phi_atoms, psi_atoms = [], []
    for position in positions:
        atom = {r: {a.name: a.index for a in residues[r].atoms}
                for r in [position - 1, position, position + 1]}
        phi_atoms.append([atom[position - 1]["C"], atom[position]["N"],
                          atom[position]["CA"], atom[position]["C"]])
        psi_atoms.append([atom[position]["N"], atom[position]["CA"],
                          atom[position]["C"], atom[position + 1]["N"]])
    phi = np.rad2deg(md.compute_dihedrals(traj, np.asarray(phi_atoms), periodic=False))
    psi = np.rad2deg(md.compute_dihedrals(traj, np.asarray(psi_atoms), periodic=False))
    phi_low, phi_high = settings["alpha_phi_deg"]
    psi_low, psi_high = settings["alpha_psi_deg"]
    helix = (phi >= phi_low) & (phi <= phi_high) & (psi >= psi_low) & (psi <= psi_high)
    return helix.mean(axis=1), helix


def hbond_fraction(traj, residues: dict[int, Any], start: int, end: int, settings: dict):
    import mdtraj as md
    import numpy as np
    pairs, angles = [], []
    for acceptor in range(start, end - 3):
        donor = acceptor + 4
        left = {a.name: a.index for a in residues[acceptor].atoms}
        right = {a.name: a.index for a in residues[donor].atoms}
        if "H" not in right:
            continue
        pairs.append([left["O"], right["N"]])
        angles.append([right["N"], right["H"], left["O"]])
    if not pairs:
        raise ValueError("No core i,i+4 donors available")
    distance = md.compute_distances(traj, np.asarray(pairs), periodic=False)
    angle = np.rad2deg(md.compute_angles(traj, np.asarray(angles), periodic=False))
    bonded = (distance <= settings["hbond_acceptor_donor_nm"]) & (
        angle >= settings["hbond_n_h_o_angle_deg_min"])
    return bonded.mean(axis=1), len(pairs)


def periodic_self_contacts(traj, ca_indices: list[int], positions: list[int], settings: dict):
    import numpy as np
    if traj.unitcell_lengths is None:
        raise ValueError("Periodic DCD box vectors missing")
    if not np.allclose(traj.unitcell_angles, 90, atol=0.1):
        raise ValueError("Self-image diagnostic expects an orthorhombic box")
    pairs = [(i, j) for i in range(len(positions)) for j in range(i + 1, len(positions))
             if abs(positions[i] - positions[j]) >= settings["image_contact_min_sequence_separation"]]
    left = np.asarray([ca_indices[i] for i, _ in pairs])
    right = np.asarray([ca_indices[j] for _, j in pairs])
    result = np.zeros(len(traj), dtype=float)
    cutoff = settings["image_contact_ca_distance_nm"]
    for offset in range(0, len(traj), 256):
        block = traj.xyz[offset:offset + 256]
        delta = block[:, left] - block[:, right]
        box = traj.unitcell_lengths[offset:offset + len(block), None, :]
        direct = np.linalg.norm(delta, axis=-1)
        image = np.linalg.norm(delta - np.round(delta / box) * box, axis=-1)
        result[offset:offset + len(block)] = ((image < cutoff) & (direct > cutoff + 0.1)).sum(axis=1)
    return result


def analyze_task(task_dir: Path, *, allow_smoke: bool = False) -> dict:
    import mdtraj as md
    import numpy as np

    build_path, run_path = task_dir / "build.json", task_dir / "run.json"
    build, run = read_json(build_path), read_json(run_path)
    task = build["task"]
    settings = build["config"]["analysis"]
    if not run["completed"]:
        raise ValueError(f"Task incomplete: {task['task_id']}")
    if run["smoke"] and not allow_smoke:
        raise ValueError("Smoke trajectory excluded from scientific analysis")
    trajectory_path = task_dir / "trajectory.dcd"
    traj = md.load_dcd(str(trajectory_path), top=str(task_dir / "protein.pdb"))
    if not np.isfinite(traj.xyz).all():
        raise ValueError("Nonfinite trajectory coordinates")
    if traj.n_atoms != build["protein_atoms"]:
        raise ValueError("DCD atom count does not match protein-only topology")
    mapping = build["protein_residue_map"]
    if traj.topology.n_residues != len(mapping):
        raise ValueError("Residue mapping length does not match trajectory")
    residues = {}
    for entry, residue in zip(mapping, traj.topology.residues):
        if residue.name != entry["residue_name"]:
            raise ValueError("Residue-name mapping mismatch")
        if entry["mane_position"] is not None:
            residues[entry["mane_position"]] = residue
    start_step, interval = run["production_start_step"], run["report_interval_steps"]
    first_report = ((start_step // interval) + 1) * interval
    times_ns = (first_report + np.arange(len(traj)) * interval - start_step) * run["timestep_fs"] / 1_000_000
    expected_frames = ((start_step + run["completed_production_steps"]) // interval - start_step // interval)
    if len(traj) != expected_frames:
        raise ValueError("Trajectory frame count does not match completed production")
    mask = times_ns > (0 if run["smoke"] else settings["discard_production_ns"])
    traj, times_ns = traj[mask], times_ns[mask]
    if len(traj) == 0:
        raise ValueError("No post-burn-in trajectory frames")
    core_start, core_end = settings["core"]
    core = list(range(core_start, core_end + 1))
    window_start, window_end = settings["site_window"]
    traces = {}
    traces["core_alpha_torsion_fraction"], per_residue_alpha = alpha_fraction(traj, residues, core, settings)
    traces["site_alpha_torsion_fraction"], _ = alpha_fraction(
        traj, residues, list(range(window_start, window_end + 1)), settings)
    traces["core_i_i4_hbond_fraction"], hbond_donors = hbond_fraction(
        traj, residues, core_start, core_end, settings)
    heavy = {i: [a.index for a in residues[i].atoms if a.element.symbol != "H"] for i in core}
    site = settings["site"]
    traces["site_nonlocal_heavy_contact_count"] = np.zeros(len(traj))
    contact_rows = []
    for other in core:
        if abs(site - other) < settings["contact_min_sequence_separation"]:
            continue
        distance = min_atom_distances(traj.xyz, heavy[site], heavy[other])
        occupation = distance <= settings["contact_heavy_distance_nm"]
        traces["site_nonlocal_heavy_contact_count"] += occupation
        contact_rows.append({"task_id": task["task_id"], "metric": "site_heavy_contact",
                             "residue_1": site, "name_1": residues[site].name,
                             "residue_2": other, "name_2": residues[other].name,
                             "occupancy": float(occupation.mean()),
                             "minimum_distance_mean_nm": float(distance.mean())})
    charged = {}
    for i in core:
        residue = residues[i]
        if residue.name in CHARGE_ATOMS:
            sign, names = CHARGE_ATOMS[residue.name]
            atoms = {a.name: a.index for a in residue.atoms}
            charged[i] = (sign, [atoms[name] for name in names])
    traces["core_salt_bridge_count"] = np.zeros(len(traj))
    traces["site_salt_bridge_count"] = np.zeros(len(traj))
    for left, right in itertools.combinations(charged, 2):
        sign_l, atoms_l = charged[left]
        sign_r, atoms_r = charged[right]
        if sign_l == sign_r:
            continue
        distance = min_atom_distances(traj.xyz, atoms_l, atoms_r)
        occupation = distance <= settings["salt_bridge_charged_distance_nm"]
        traces["core_salt_bridge_count"] += occupation
        if site in {left, right}:
            traces["site_salt_bridge_count"] += occupation
        contact_rows.append({"task_id": task["task_id"], "metric": "salt_bridge",
                             "residue_1": left, "name_1": residues[left].name,
                             "residue_2": right, "name_2": residues[right].name,
                             "occupancy": float(occupation.mean()),
                             "minimum_distance_mean_nm": float(distance.mean())})
    ca_indices = [next(a.index for a in residues[i].atoms if a.name == "CA") for i in core]
    all_positions = sorted(residues)
    all_ca = [next(a.index for a in residues[i].atoms if a.name == "CA") for i in all_positions]
    traces["periodic_self_contact_pair_count"] = periodic_self_contacts(traj, all_ca, all_positions, settings)
    reference = md.load_pdb(str(task_dir / "protein.pdb"))
    traces["core_ca_rmsd_nm"] = md.rmsd(traj, reference, atom_indices=np.asarray(ca_indices))
    endpoint_keys = settings["primary_endpoints"] + settings["secondary_endpoints"] + settings["quality_endpoints"]
    rows = [{"production_time_ns": float(time), **{key: float(traces[key][j]) for key in endpoint_keys}}
            for j, time in enumerate(times_ns)]
    write_tsv(task_dir / "metrics_timeseries.tsv", rows)
    write_tsv(task_dir / "contact_occupancies.tsv", contact_rows)
    write_tsv(task_dir / "alpha_residue_occupancies.tsv", [
        {"task_id": task["task_id"], "mane_position": i,
         "alpha_torsion_occupancy": float(per_residue_alpha[:, j].mean())}
        for j, i in enumerate(core)])
    blocks = []
    block_indices = np.floor((times_ns - (0 if run["smoke"] else settings["discard_production_ns"])) /
                             settings["block_ns"]).astype(int)
    expected_block_frames = settings["block_ns"] * 1_000_000 / (interval * run["timestep_fs"])
    for number in np.unique(block_indices):
        selected = block_indices == number
        if selected.sum() < 0.95 * expected_block_frames:
            continue
        blocks.append({"block": int(number), "frames": int(selected.sum()),
                       **{key: float(traces[key][selected].mean()) for key in endpoint_keys}})
    eligible = (not run["smoke"] and run["duration_override_ns"] is None
                and run["completed_production_steps"] == run["target_production_steps"])
    window = settings.get("comparison_window_ns")
    comparison_mask = ((times_ns > window[0]) & (times_ns <= window[1] + 1e-9)
                       if window and eligible else np.ones(len(traj), dtype=bool))
    if not comparison_mask.any():
        raise ValueError("No frames in the matched comparison window")
    image_fraction = float((traces["periodic_self_contact_pair_count"] > 0).mean())
    summary = {
        "task": task, "eligible_for_inference": eligible,
        "frames": len(traj), "analysis_time_ns": [float(times_ns[0]), float(times_ns[-1])],
        "hbond_possible_donors": hbond_donors, "complete_blocks": blocks,
        "periodic_self_contact_frame_fraction": image_fraction,
        "endpoints": {key: {"mean": float(traces[key].mean()),
                             "std_frames_descriptive": float(traces[key].std()),
                             "first_vs_second_half_difference": float(
                                 traces[key][len(traj) // 2:].mean() - traces[key][:max(1, len(traj) // 2)].mean()),
                             **autocorrelation_quality(traces[key])} for key in endpoint_keys},
        "comparison_window_ns": window,
        "comparison_frames": int(comparison_mask.sum()),
        "comparison_endpoints": {key: {"mean": float(traces[key][comparison_mask].mean())}
                                 for key in endpoint_keys},
        "claim_limit": build["config"]["interpretation"]["claim_limit"],
        "static_energy_excluded_from_contrasts": True,
        "config_sha256": task["config_sha256"], "source_sha256": build["source_sha256"],
        "af_confidence": build.get("af_confidence"),
    }
    if image_fraction > 0:
        summary["quality_warning"] = "Periodic solute self-contact detected; inspect and expand solvent box before mechanistic interpretation."
    write_json(task_dir / "analysis.json", summary)
    for name in ["metrics_timeseries.tsv", "contact_occupancies.tsv", "alpha_residue_occupancies.tsv", "analysis.json"]:
        provenance(task_dir / name, [build_path, run_path, trajectory_path], params={"task": task, "analysis": settings})
    return summary


def compare_replicas(summaries: list[dict], config: dict) -> list[dict]:
    import numpy as np
    eligible = [item for item in summaries if item["eligible_for_inference"]]
    grouped = defaultdict(list)
    for item in eligible:
        grouped[(item["task"]["boundary"], item["task"]["variant"])].append(item)
    results = []
    settings = config["analysis"]
    rng = np.random.default_rng(settings["bootstrap_seed"])
    for boundary in config["boundaries"]:
        wildtype = grouped[(boundary, "WT")]
        for variant in config["variants"]:
            if variant == "WT":
                continue
            mutants = grouped[(boundary, variant)]
            valid = len(wildtype) == config["replicas"] and len(mutants) == config["replicas"]
            if valid:
                for group in [wildtype, mutants]:
                    if len({item["task"]["replica"] for item in group}) != config["replicas"]:
                        raise ValueError("Duplicate replica identifiers cannot be treated as independent runs")
                    if len({item["task"]["seed"] for item in group}) != config["replicas"]:
                        raise ValueError("Duplicate random seeds cannot be treated as independent runs")
            self_image_warning = any(item["periodic_self_contact_frame_fraction"] > 0
                                     for item in wildtype + mutants)
            for endpoint in settings["primary_endpoints"] + settings["secondary_endpoints"]:
                row = {"boundary": boundary, "variant": variant, "endpoint": endpoint,
                       "endpoint_role": "primary" if endpoint in settings["primary_endpoints"] else "secondary",
                       "comparison_window_ns": settings.get("comparison_window_ns"),
                       "n_wt_replicas": len(wildtype), "n_variant_replicas": len(mutants),
                       "status": ("blocked_self_image_contacts" if valid and self_image_warning else
                                  "exploratory_complete" if valid else "insufficient_replicas"),
                       "wt_mean": None, "variant_mean": None, "difference_variant_minus_wt": None,
                       "bootstrap_percentile_2_5": None, "bootstrap_percentile_97_5": None,
                       "wt_replica_range": None, "variant_replica_range": None,
                       "quality_flag": "self_image_contacts" if self_image_warning else "none"}
                if valid and not self_image_warning:
                    wt = np.asarray([item.get("comparison_endpoints", item["endpoints"])[endpoint]["mean"] for item in wildtype])
                    mut = np.asarray([item.get("comparison_endpoints", item["endpoints"])[endpoint]["mean"] for item in mutants])
                    n_boot = settings["bootstrap_iterations"]
                    boot = (mut[rng.integers(len(mut), size=(n_boot, len(mut)))].mean(axis=1) -
                            wt[rng.integers(len(wt), size=(n_boot, len(wt)))].mean(axis=1))
                    row.update({"wt_mean": float(wt.mean()), "variant_mean": float(mut.mean()),
                                "difference_variant_minus_wt": float(mut.mean() - wt.mean()),
                                "bootstrap_percentile_2_5": float(np.quantile(boot, 0.025)),
                                "bootstrap_percentile_97_5": float(np.quantile(boot, 0.975)),
                                "wt_replica_range": float(np.ptp(wt)), "variant_replica_range": float(np.ptp(mut))})
                results.append(row)
    return results


def boundary_sensitivity(comparisons: list[dict], config: dict) -> list[dict]:
    """Do not promote an inconsistent local-fragment result to a robust mechanism."""
    rows = []
    for variant in config["variants"]:
        if variant == "WT":
            continue
        for endpoint in config["analysis"]["primary_endpoints"]:
            selected = [item for item in comparisons if item["variant"] == variant and item["endpoint"] == endpoint]
            usable = all(item["status"] == "exploratory_complete" for item in selected)
            usable = usable and len(selected) == len(config["boundaries"])
            differences = [item["difference_variant_minus_wt"] for item in selected] if usable else []
            same_sign = usable and (all(value > 0 for value in differences) or all(value < 0 for value in differences))
            intervals_exclude_zero = usable and all(
                item["bootstrap_percentile_2_5"] > 0 or item["bootstrap_percentile_97_5"] < 0 for item in selected)
            rows.append({"variant": variant, "endpoint": endpoint,
                         "all_boundaries_complete_and_quality_pass": usable,
                         "same_effect_direction": same_sign,
                         "all_exploratory_intervals_exclude_zero": intervals_exclude_zero,
                         "status": "locally_consistent_exploratory_signal" if same_sign and intervals_exclude_zero else
                                   "incomplete_or_boundary_dependent_or_unresolved",
                         "clinical_interpretation": "none"})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--task-dir", type=Path)
    group.add_argument("--run-root", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--allow-smoke", action="store_true", help="Diagnostic only, excludes inference")
    args = parser.parse_args()
    if args.task_dir:
        summary = analyze_task(args.task_dir, allow_smoke=args.allow_smoke)
        print(json.dumps({"analyzed": summary["task"]["task_id"], "frames": summary["frames"],
                          "eligible_for_inference": summary["eligible_for_inference"]}))
        return
    if not args.output_root:
        parser.error("--run-root requires --output-root")
    config = load_config(args.config)
    summaries = []
    for run_path in sorted(args.run_root.glob("*/run.json")):
        state = read_json(run_path)
        if not state["completed"] or state["smoke"] or state["duration_override_ns"] is not None:
            continue
        summary = analyze_task(run_path.parent)
        if summary["config_sha256"] != sha256(args.config):
            raise ValueError("Cannot combine different frozen protocol hashes")
        summaries.append(summary)
    if not summaries:
        raise ValueError("No completed scientific trajectories available")
    if len({item["source_sha256"] for item in summaries}) != 1:
        raise ValueError("Cannot combine different source models")
    comparisons = compare_replicas(summaries, config)
    output = args.output_root / "KIF13A_K1122E_replica_comparisons_20261004.tsv"
    write_tsv(output, comparisons)
    boundary_output = args.output_root / "KIF13A_K1122E_boundary_sensitivity_20261004.tsv"
    write_tsv(boundary_output, boundary_sensitivity(comparisons, config))
    manifest_output = args.output_root / "KIF13A_K1122E_analysis_manifest_20261004.json"
    write_json(manifest_output, {
        "protocol": config, "tasks": [item["task"] for item in summaries],
        "inference_limit": config["analysis"]["inference"],
        "boundary_sensitivity": "Compare matched production 10-50 ns windows across primary and padded fragments. Full-duration endpoints remain descriptive. Discordant signs or self-image warnings block a robust local mechanism claim.",
        "energies": "No total-potential WT-mutant subtraction. No binding or folding ddG.",
    })
    provenance(output, [args.config] + [args.run_root / item["task"]["task_id"] / "analysis.json" for item in summaries],
               params={"bootstrap_iterations": config["analysis"]["bootstrap_iterations"]})
    provenance(boundary_output, [args.config, output], params={"gate": "same sign and no self-image contacts"})
    provenance(manifest_output, [args.config, output, boundary_output], params={"completed_tasks": len(summaries)})
    print(json.dumps({"comparison_table": str(output), "completed_tasks": len(summaries)}))


if __name__ == "__main__":
    main()
