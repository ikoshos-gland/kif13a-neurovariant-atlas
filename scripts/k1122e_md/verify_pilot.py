"""Require completed ff19SB/OPC CUDA pilots before the main array can run."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from common import DEFAULT_CONFIG, provenance, read_json, sha256, write_json


def verify(task_dir: Path, config_path: Path) -> dict:
    build = read_json(task_dir / "build.json")
    run = read_json(task_dir / "run.json")
    if build["task"]["config_sha256"] != sha256(config_path) or run["config_sha256"] != sha256(config_path):
        raise ValueError("Pilot and current protocol hashes differ")
    if run["platform"] != "CUDA" or run["smoke"] or not run["completed"]:
        raise ValueError("Need a completed full-equilibration CUDA pilot")
    if run["duration_override_ns"] != 0.1 or run["completed_production_steps"] != run["target_production_steps"]:
        raise ValueError("Pilot production duration or completion mismatch")
    if not build["solvent_validation"].get("all_waters_passed") or build["ff19sb_cmap_torsions"] <= 0:
        raise ValueError("OPC water or ff19SB backbone check failed")
    with (task_dir / "states.csv").open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = [[float(value) for value in row] for row in reader]
    if not rows or not all(math.isfinite(value) for row in rows for value in row):
        raise ValueError("Empty or nonfinite pilot state report")
    temperature_column = next(j for j, name in enumerate(header) if "Temperature" in name)
    volume_column = next(j for j, name in enumerate(header) if "Volume" in name)
    temperatures = [row[temperature_column] for row in rows]
    if not all(290 <= value <= 330 for value in temperatures) or not all(row[volume_column] > 0 for row in rows):
        raise ValueError("Pilot temperature or volume outside the preflight bounds")
    result = {"task_id": build["task"]["task_id"], "passed": True,
              "protocol_sha256": sha256(config_path), "forcefield": build["config"]["simulation"]["forcefield"],
              "system_particles": build["system_atoms"], "solvent": build["solvent_validation"],
              "ff19sb_cmap_torsions": build["ff19sb_cmap_torsions"],
              "temperature_range_k": [min(temperatures), max(temperatures)],
              "equilibration_ns": 1.2, "production_ns": 0.1,
              "wall_seconds": run["wall_seconds"], "scientific_inference_eligible": False}
    output = task_dir / "pilot_validation.json"
    write_json(output, result)
    provenance(output, [config_path, task_dir / "build.json", task_dir / "run.json", task_dir / "states.csv"],
               params={"temperature_bounds_k": [290, 330], "purpose": "engine preflight only"})
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    args = parser.parse_args()
    print(json.dumps(verify(args.task_dir, args.config)), flush=True)


if __name__ == "__main__":
    main()
