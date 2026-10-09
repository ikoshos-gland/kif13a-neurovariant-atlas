"""Frozen KIF13A monomer protocol, task mapping and portable provenance helpers.

No scientific dependency is needed to validate or print the 30-task manifest:
    python scripts/k1122e_md/common.py --output tasks.json
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_CONFIG = HERE / "protocol.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(64 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = read_json(path)
    if config["schema_version"] != 1:
        raise ValueError("Unsupported protocol schema")
    core_start, core_end = config["analysis"]["core"]
    if not core_start <= config["analysis"]["site"] <= core_end:
        raise ValueError("Site must be inside the analysis core")
    if config["replicas"] < 3:
        raise ValueError("Protocol requires at least three independent replicas")
    for name, (start, end) in config["boundaries"].items():
        if not 1 < start < core_start < core_end < end < config["source"]["protein_length"]:
            raise ValueError(f"Invalid padded fragment: {name}")
        for mutation in config["variants"].values():
            if mutation and not start <= mutation["position"] <= end:
                raise ValueError("Mutation lies outside a fragment")
    sim = config["simulation"]
    if sim["box_policy"] != "cubic_max_initial_atom_extent_plus_padding":
        raise ValueError("The frozen protocol requires a cubic solvent box")
    if not 0 < config["analysis"]["discard_production_ns"] < sim["production_ns"]:
        raise ValueError("Production burn-in must be shorter than production")
    window = config["analysis"].get("comparison_window_ns")
    if window:
        durations = sim.get("production_ns_by_boundary", {}).values()
        if not config["analysis"]["discard_production_ns"] <= window[0] < window[1] <= min(durations or [sim["production_ns"]]):
            raise ValueError("Comparison window must fit every boundary after burn-in")
    if sim["water_model"] == "opc" and (sim.get("solvent_template_model") != "tip4pew"
            or sim["forcefield"] != ["amber19/protein.ff19SB.xml", "amber19/opc.xml"]):
        raise ValueError("OPC main model requires ff19SB, OPC ions and a four-site coordinate template")
    return config


def task_manifest(config: dict[str, Any], config_path: Path) -> list[dict[str, Any]]:
    tasks = []
    config_hash = sha256(config_path)
    for boundary, residues in config["boundaries"].items():
        for variant, mutation in config["variants"].items():
            for replica in range(1, config["replicas"] + 1):
                index = len(tasks)
                tasks.append({
                    "task_index": index,
                    "task_id": f"{boundary}_{variant}_r{replica:02d}",
                    "boundary": boundary,
                    "residues": residues,
                    "variant": variant,
                    "mutation": mutation,
                    "replica": replica,
                    "production_ns": config["simulation"].get("production_ns_by_boundary", {}).get(
                        boundary, config["simulation"]["production_ns"]),
                    "seed": config["seed_base"] + 100 * index + 1,
                    "protocol_id": config["protocol_id"],
                    "config_sha256": config_hash,
                })
    return tasks


def task_at(config_path: Path, task_index: int) -> tuple[dict[str, Any], dict[str, Any]]:
    config = load_config(config_path)
    tasks = task_manifest(config, config_path)
    if not 0 <= task_index < len(tasks):
        raise ValueError(f"Task index must be in 0..{len(tasks) - 1}")
    return config, tasks[task_index]


def relative_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return str(path.resolve())


def provenance(path: Path, inputs: list[Path], *, params: dict[str, Any],
               source: str = "KIF13A local monomer calculation") -> None:
    """Track derived outputs without modifying the raw-download ledger.

    The raw AlphaFold inputs and their existing sidecars are preserved. Simulation
    outputs are not raw data pulls, so each receives a derivative .prov.json.
    """
    record = {
        "source": source,
        "tool": "scripts/k1122e_md",
        "utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "path": relative_path(path),
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
        "inputs": [{"path": relative_path(item), "sha256": sha256(item)} for item in inputs],
        "params": params,
    }
    write_json(Path(str(path) + ".prov.json"), record)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    tasks = task_manifest(load_config(args.config), args.config)
    if args.output:
        write_json(args.output, tasks)
        provenance(args.output, [args.config], params={"tasks": len(tasks)})
    else:
        print(json.dumps(tasks, indent=2))


if __name__ == "__main__":
    main()
